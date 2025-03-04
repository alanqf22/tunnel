# -*- coding: utf-8 -*-


## /********************************************************/
## /********************************************************/
##                         CONSTANTS       
## Function: Contains the constants used by Eyringpy
## according to the NIST (National Institute of Stan-
## dars and Technology, http:///www.physics.nist.gov/
## cuu/Constants/).
##
## References:
##
## 1. Mohr, P. J.; Newell, D. B.; Taylor, B. N. Rev.
##    Mod. Phys. 2016, 88 (3), 1-73.
## /********************************************************/
## /********************************************************/

## Constants
ATOMIC_MASS_CONSTANT = 1.66053906660E-027    # kg
AVOGADRO_CONSTANT = 6.022140857E23           # mol-1
BOLTZMANN_CONSTANT = 1.38064852E-23          # J K-1
IDEAL_GAS_CONSTANT = 8.3144598               # J mol-1 K-1
PLANCK_CONSTANT = 6.626070040E-34            # J s-1
PLANCK_CONSTANT_OVER_2PI = 1.054571800E-34   # J s-1
SPEED_OF_LIGHT = 299792458.0                 # m s-1
STANDARD_PRESSURE = 101325.0                 # Pa
BOHR_RADIUS = 0.52917721067E-10              # m
ELECTRON_MASS = 9.1093837015E-31             # kg
PI = 3.141592653589793					 	
TWOPI = 2*PI

## Conversions

ANGSTROMS_TO_METERS = 1E-10
CALORIES_TO_JOULES = 4.184    # cal to J
CALORIES_MOL_TO_JOULES = CALORIES_TO_JOULES / AVOGADRO_CONSTANT    							# cal mol-1 to J
CUBICMETER_TO_LITERS = 1000.0    															# m3 to L
HARTREES_TO_JOULES = 4.359744650E-18    													# Hartrees to J
HARTREES_TO_CALORIES_MOL = (HARTREES_TO_JOULES * AVOGADRO_CONSTANT) / CALORIES_TO_JOULES    # Hartrees to cal mol-1
HARTREES_TO_KCAL_MOL = HARTREES_TO_CALORIES_MOL/1000    									# Hartrees to cal mol-1
HARTREES_TO_ELECTRONVOLT = 27.211396132 													# Hartrees to eV
IDEAL_GAS_CONSTANT_CALORIES = IDEAL_GAS_CONSTANT / CALORIES_TO_JOULES    					# cal mol-1 K-1
LITER_MOL_TO_CUBICCENTIMETER_MOLECULE = 1000.0 / AVOGADRO_CONSTANT    						# L mol-1 to cm3 molecule-1
BOLTZMANN_CONSTANT_KCAL = (BOLTZMANN_CONSTANT * AVOGADRO_CONSTANT) / 4184.0   				# kcal mol-1
BOHR_RADIUS_CM = 100*0.52917721067E-10              												# cm

## a.u.
ATOMIC_MASS_UNIT = ELECTRON_MASS/ATOMIC_MASS_CONSTANT										
MU = 1/ATOMIC_MASS_UNIT
PLANCK_CONSTANT_OVER_2PI_AU = 1              												# a.u.
BOLTZMANN_CONSTANT_AU = BOLTZMANN_CONSTANT/HARTREES_TO_JOULES								# a.u.
SPEED_OF_LIGHT_AU = SPEED_OF_LIGHT / (BOHR_RADIUS /(1.0 / 4.1341373337e+16)) 				# a.u.										# a.u.
