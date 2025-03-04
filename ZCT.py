#Python libraries
#import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from scipy.interpolate import CubicSpline
from scipy.interpolate import Rbf
from scipy.signal import savgol_filter

import gmpy2
gmpy2.get_context().precision = 256

#Eyringpy modules
from spline import *
from criteria import *
from plot_map import *
from constants import *
from numerical_methods import *
from tunneling import *

def divide_coordinates(x, y, xref = "none"):
	"""
    Divides reaction path into two sides (left and right).

    Parameters:
    x      : array_like
             Reaction coordinate.
    y      : array_like
             Energy profile.
    xref   : scalar
             x value that divides the x-axis into left and right.

    Returns:
    xleft  : array_like
             x values lower than xref
    yleft  : array_like
             y values corresponding to xleft
    xright : array_like
             x values higher than xref
    yright  : array_like
             y values corresponding to xright
    """
	if xref == "none":
		ymax = max(y)
		xref = x[ymax_coordiante.index(ymax)]

	xleft = [ ]
	xright = [ ]
	yleft = [ ]
	yright = [ ]

	for i in range(len(x)):
		if x[i] <= xref:
			xleft.append(x[i])
			yleft.append(y[i])
		if x[i] >= xref:
			xright.append(x[i])
			yright.append(y[i])

	divided_reaction_path = { }
	divided_reaction_path['xleft'] = xleft
	divided_reaction_path['xright'] = xright
	divided_reaction_path['yleft'] = yleft
	divided_reaction_path['yright'] = yright

	return divided_reaction_path

def get_turning_points(s,VGa,E):
	"""
    Get the classical turning points eference a reference energy value.

    Parameters:
    s   : array_like
          Reaction coordinate.
    V   : array_like
          Energy profile. Usually, it is the vibrationally
          adiabatic ground-state potential energies.
    E   : scalar
          Reference energy value.

    Returns:
    VGa : array_like
          The result of summing V + zpe. 
    """
	potential_data = VGa_critical_points(s, VGa)
	VAG_s = potential_data['VAG_s']		#Coordinate of the maximum of the VGa curve
	path_divided = divide_coordinates(s, VGa, xref = VAG_s)
	s_left = path_divided['xleft']
	s_right = path_divided['xright']
	VGa_left = path_divided['yleft']
	VGa_right = path_divided['yright']

	s1 = "none"
	s2 = "none"

	turning_points = [ ]
	dE = 0.01
	for i in range(len(VGa_left)):
		if VGa_left[i] < E:
			continue
		else:
			if abs(VGa_left[i] - E) <= dE:
				dE = abs(E - VGa_left[i])
				s1 = s_left[i]
	turning_points.append(s1)

	dE = 0.01
	for j in range(len(VGa_right)):
		if VGa_left[i] < E:
			continue
		else:
			if abs(VGa_right[j]-E) <= dE:
				dE = abs(E - VGa_right[j])
				s2 = s_right[j]

	if s2 == "none":
		s2 = s1
	
	turning_points.append(s2)

	return turning_points

###--------------------------------------------------------------------------------
#### Functions to find the transmission coefficient

def theta(s_coords,VGa_energies,mueff,E):
	#Function to get the action integral (theta):

	#Find the turning points:
	turning_points = get_turning_points(s_coords,VGa_energies,E)

	s1 = turning_points[0]	#Minimum turning point
	s2 = turning_points[1]	#Maximum turning point

	#Calculate the action integral (theta):
	all_theta_integrand = [ ]
	for i in range(len(VGa_energies)):
		theta_integrand = np.sqrt(2*mueff[i]*abs(E-VGa_energies[i]))		
		all_theta_integrand.append(theta_integrand)

	theta_integral = cumulative_trapz(s_coords,all_theta_integrand,xmin=s1,xmax=s2)

	return theta_integral

def tunneling_probability(s_coords,VGa_energies,mueff,E,E0,VAG):
	#Function to calculate the tunneling probability:
	if E < E0:
		probability = 0
	elif E0 <= E <= VAG:
		theta_int = 2*theta(s_coords,VGa_energies,mueff,E)
		if theta_int < -700 or theta_int > 700:
			exponential = np.exp(np.float128(theta_int))
		else:
			exponential = np.exp(theta_int)
		if exponential == "inf": 
			probability = 0
		else:
			probability = 1/(1+exponential)
	elif E > (2*VAG-E0):
		probability=1
	else:
		theta_int = theta(s_coords,VGa_energies,mueff,2*VAG-E)
		exponential = gmpy2.exp(gmpy2.mpfr(2*theta_int))
		probability = 1 - (1/(1 + exponential))

	return probability


def kappa(s_coords,VGa_energies,E0, VAG, rxvts=0.0, T="none", mueff = [ ]):
	mu = 1/(9.1093837015E-031/1.66053906660E-027)
	BETA = 1/(BOLTZMANN_CONSTANT_AU * T)

	if len(mueff) == 0:
		mu = 1/(9.1093837015E-031/1.66053906660E-027)
		mueff = [mu] * len(VGa_energies)

	if T == "none":
		T = 298.15 #Kelvin

	mueff_all = [ ]
	for i in mueff:
		if i < mu:
			mueff_all.append(i)
		else:
			mueff_all.append(mu)

	# Interpolate reaction coordinate s
	s_coords_n = np.linspace(min(s_coords),max(s_coords),1000)
	s_coords_new = np.append(s_coords_n, rxvts)
	s_coords_new.sort()

	# Interpolate VGa(s)	
	VGa_energies_new = spline(s_coords,VGa_energies,s_coords_new,interp_order=3)

	# Interpolate mueff(s)
	mueff_new = spline(s_coords,mueff_all,s_coords_new,interp_order=1)

	# Calculate mueff(s)/mu
	mueff_mu = [imueff/mu for imueff in mueff_new]

	# ************************** Tunneling ************************** #
	#Etun_list, tun_weights = gauquad_pointsweights(30,VGa_energies[0],VAG)
	Etun_list, tun_weights = gauquad_pointsweights(30,E0,VAG)
	all_Etun_prob = [ ]
	Etun_cont = [ ]

	kappa_tun = 0.0
	for i in range(len(Etun_list)):
		tun_prob = tunneling_probability(s_coords_new,VGa_energies_new,mueff_new,Etun_list[i],E0,VAG)
		tun_integrand = tun_prob * gmpy2.exp(gmpy2.mpfr(-BETA*(Etun_list[i]-VAG))) * BETA
		all_Etun_prob.append(tun_prob)

		Etun_cont.append(tun_weights[i] * tun_integrand)
		kappa_tun += tun_weights[i] * tun_integrand

	all_Etun_cont = [ ]
	for ituncont in Etun_cont:
		all_Etun_cont.append((ituncont/kappa_tun)*100)

	# ************************** Non-classical reflection ************************** #
	Eref_list, ref_weights = gauquad_pointsweights(30,VAG,(2*VAG-E0))
	
	kappa_ref = 0.0
	for i in range(len(Eref_list)):
		ref_prob = tunneling_probability(s_coords_new,VGa_energies_new,mueff_new,Eref_list[i],E0,VAG)
		ref_integrand2 = ref_prob * gmpy2.exp(gmpy2.mpfr(-BETA*(Eref_list[i]-VAG))) * BETA
		kappa_ref += ref_weights[i] * ref_integrand2

	# kappa correction for VTST
	position_rxvts = np.where(s_coords_new == rxvts)
	VAG_vts = VGa_energies_new[position_rxvts[0][0]]
	kcvt_cag =  np.exp(BETA*(VAG_vts - VAG))

	# kappa = kappa from tunneling + kappa from non-classical reflection
	kappa_total = kappa_tun + kappa_ref

	# Convert VGa and Etun from HA to kcal mol-1
	VGa_energies_new_kcal = [iVG * HARTREES_TO_KCAL_MOL for iVG in VGa_energies_new]
	Etun_list_kcal = [jVG * HARTREES_TO_KCAL_MOL for jVG in Etun_list]

	kappa_data = { }
	kappa_data = {'scoord':s_coords_new,'VGa':VGa_energies_new_kcal, 'mueff_mu':mueff_mu,
	'Etun':Etun_list_kcal, 'Etunprob':all_Etun_prob, 'Etuncont':all_Etun_cont,
	'kappatun':kappa_tun, 'kapparef':kappa_ref, 'kappa':kappa_total,'kcvt_cag':kcvt_cag,}

	return kappa_data

