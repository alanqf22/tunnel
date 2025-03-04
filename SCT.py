# -*- coding: utf-8 -*-

#Python libraries
import numpy as np
import matplotlib.pyplot as plt

#Eyring modules
from funcs import *
from hessian import *
from constants import *
from numerical_methods import *

## /********************************************************/
## /********************************************************/
##                          SCT 
## Function: calculates the terms to solve the SCT model.
## /********************************************************/
## /********************************************************/

##### -------------------------------------------------------------------- #####
def kappa_corrected(curvatura, umbral = 1, max_kappa = 0.5):
    """
    Corrige los valores de curvatura (kappa) grandes para que se mantengan dentro de un umbral.
    
    Parámetros:
    - curvatura: lista de valores de curvatura (kappa)
    - umbral: valor a partir del cual kappa es considerado grande
    - max_kappa: valor máximo permitido para la curvatura después de la corrección
    
    Retorna:
    - curvatura_corregida: lista de kappa corregida
    """
    curvatura_array = np.array(curvatura)

    # Crear una copia para almacenar la curvatura corregida
    curvatura_corregida = np.copy(curvatura_array)

    # Encontrar los índices donde la curvatura es mayor que el umbral
    indices_grande_kappa = np.where(curvatura_array > umbral)
    
    # Aplicar la corrección en los índices identificados
    curvatura_corregida[indices_grande_kappa] = max_kappa + np.log(curvatura_array[indices_grande_kappa] / umbral)

    return curvatura_corregida


def page_mciver(gms1, Fms):
    #Function: Calculates curvature vector by Page and McIver method
    if len(gms1) == 0 or len(Fms) == 0: return None
    #REF: Page, M., & McIver Jr, J. W. (1988).
    #On evaluating the reaction path Hamiltonian.
    #The Journal of chemical physics, 88(2), 922-935

    gms = [ ]
    for i in gms1:
        gms.append(i)

    #Calculate v0:
    #Equation:n v0 = gms_i/c ; c = norm(gms)
    c = np.linalg.norm(gms)
    v0 = [-g_i/c for g_i in gms]

    #Calculate v1:
    #Equation: v1 = [Fms * v0 - (v0^T * Fms * v0) * v0]/c
    minuend = np.array(Fms * np.matrix(v0).transpose()).transpose()[0]

    subtrahend = float(np.matrix(v0) * Fms * np.matrix(v0).transpose())
    subtrahend2 = [subtrahend*iv for iv in v0]
    v1 = (minuend - subtrahend2) / c 
    
    return v1

##### -------------------------------------------------------------------- #####

def get_BmF(evecs,v1,coord):
    #Function: Calculates the reaction-path curvature coupling elements
    #Equation: BmF = -[sign(s)]*Lm^T(s)*F(s)*v1(s)
    BmF = [-sign_of_number(coord)*np.dot(Lm,v1) for Lm in evecs]
    return BmF

##### -------------------------------------------------------------------- #####

def get_turn_points(freqs):
    #Function: Calculates turning point of a harmonic potential of frequency
    #Equation: tm(s) = [h_bar/(mu*angfreq(s)]^1/2
    #mu = 1/ATOMIC_MASS_UNIT
    mu = 1/(9.1093837015E-031/1.66053906660E-027)
    all_turn_points = [ ]
    for angfreq in freqs:
        if angfreq < 0.0:
            turn_point = 1e10
        else:
            turn_point = np.sqrt(PLANCK_CONSTANT_OVER_2PI_AU / angfreq / mu)
        all_turn_points.append(turn_point)
    return all_turn_points

##### -------------------------------------------------------------------- #####

def get_kappa(BmF):
    #Function: Calculates the reaction path curvature (kappa)
    #Equation: k(s) = (sum(BmF(s)^2))^1/2
    kappa = np.sqrt(sum([ibmf**2 for ibmf in BmF]))
    return kappa

##### -------------------------------------------------------------------- #####

def get_tbar(BmF, tpoints, kappa):
    #Function: Calculates the associated turning point for zero-point motion in the harmonic potential (tbar)
    tbar = sum([(ibmf**2)/(itm**4) for (ibmf,itm) in zip(BmF,tpoints)]) ** (-0.25) * np.sqrt(kappa)
    return tbar

##### -------------------------------------------------------------------- #####

def mu_eff(all_Gcc, all_Fcc, all_ccords, all_masses, all_rxcoord, degrees):
    #Function: Calculates the effective reduced mass for Small-Curvature Tunneling Model
    all_k = [ ]         #Reaction path curvature
    all_tbar = [ ]      #Associated turning point for zero-point motion in the harmonic potential
    all_zpe = [ ]       #Zero-Point Energy

    for idata in range(len(all_Gcc)):
        rxcoord = all_rxcoord[idata]

        if rxcoord != 0.0 or rxcoord == 0.0:
            Gcc=all_Gcc[idata]
            Fcc=all_Fcc[idata]
            ccords = all_ccords[idata]
            masses = all_masses[idata]

            if rxcoord == 0.0:
                ccords = shift2com(ccords,masses)

            masses = [im/ATOMIC_MASS_UNIT for im in masses]
            N =len(masses)

            ##----------------- ENERGY GRADIENT -----------------##
            #Convert gradient from cartesian to mass-scaled
            mass_grad = cart2mass_grad(Gcc,masses)
            #Put all elements from Gradient into a single list:
            Gms = [ ]
            for i in mass_grad:
                Gms.extend(i)
            Gms = np.array(Gms).transpose()

            ##----------------- FORCE CONSTANTS MATRIX -----------------##
            #Create a square matrix (Hessian) of size 3Nx3N
            hessian = low2matrix (3*N,Fcc)
            # Convert Hessian matrix from cartesian to mass-scaled
            Fms = cart2mass_force(hessian,masses)

            #Diaginalize of the Hessian matrix with translations and rotations
            evals1, evecs1, freqs1 = diagonalize(Fms)

            #Frequencies (with translations and rotations) in cm^-1
            all_freq_cm1_NO = [ ]
            for ifr in freqs1:
                freq_cm1_NO = ifr/TWOPI/SPEED_OF_LIGHT_AU/(BOHR_RADIUS_CM)
                all_freq_cm1_NO.append(freq_cm1_NO)

            if rxcoord == 0.0:
                #Create projection matrix in the rotating and translating frame for the Hessian matrix
                proj_matrix = get_projectionmatrix(ccords,masses,None)

            else:
                #Create projection matrix in the rotating and translating frame for the Hessian matrix
                proj_matrix = get_projectionmatrix(ccords,masses,Gms)

            #Remove the rotating and translating frame from the Hessian matrix
            Fic = project_hessian(Fms, N, proj_matrix)
            Fic = np.array(Fic)

            #Diaginalize of the Hessian matrix without translations and rotations
            evals, evecs, freqs = diagonalize(Fic)

            if rxcoord == 0.0:
                evals, evecs, freqs = remove_extra_freqs(evals, evecs, freqs, Gms, degrees)

            else:
                evals, evecs, freqs = remove_extra_freqs(evals, evecs, freqs, Gms, degrees)
     
            #Frequencies (without translations and rotations) in cm^-1
            all_freq_cm1 = [ ]
            for ifre in freqs:
                freq_cm1 = ifre/TWOPI/SPEED_OF_LIGHT_AU/(BOHR_RADIUS_CM)
                all_freq_cm1.append(freq_cm1)

            #Calculate zero point energy in HA
            zpe = 0.0
            for ifr in range(len(freqs)):
                if freqs[ifr] >= 0.0:
                    zpe += 0.5*freqs[ifr]
            all_zpe.append(zpe)

            #Calculate curvature vector:
            v1 = page_mciver(Gms, Fms)
        
            #Calculate reaction-path curvature coupling elements:
            BmF = get_BmF(evecs,v1,rxcoord)

            #Calculate turning point of a harmonic potential of frequency:
            tpoints = get_turn_points(freqs)

            #Calculate reaction path curvature:
            k = get_kappa(BmF)
            all_k.append(k)

            #Calculate t_bar:
            tbar = get_tbar(BmF, tpoints, k)
            all_tbar.append(tbar)
            
        else:
            all_k.append(None)
            all_tbar.append(None)

    #Interpolates kappa and tbar for TS:
    ts_pos = None
    if 0.0 in all_rxcoord:
        ts_pos = all_rxcoord.index(0.0)     #Position TS
    if ts_pos != None:
        all_k[ts_pos] = (all_k[ts_pos+1] + all_k[ts_pos-1])/2
        all_tbar[ts_pos] = (all_tbar[ts_pos+1] + all_tbar[ts_pos-1])/2

    #Derivates tbar:
    all_dtbar = finite_difference_mixed(all_rxcoord, all_tbar)
    
    #Calculate effective reduced mass:
    arg = [-2*ik*it - (ik*it)**2 + (idt)**2 for (ik,it,idt) in zip(all_k,all_tbar,all_dtbar)]
    all_mueff_mu = [np.exp(iarg) for iarg in arg]
    if 0.0 in all_rxcoord:
        all_mueff_mu[ts_pos] = (all_mueff_mu[ts_pos+1] + all_mueff_mu[ts_pos-1])/2
    all_mueff = [im * (1.0/ATOMIC_MASS_UNIT) for im in all_mueff_mu]

    """
    plt.figure(figsize=(8, 6), dpi=300) #Size
    plt.plot(all_rxcoord, all_k, 'ko', label=r"$\kappa$", markersize=0.8)
    plt.xlabel(r"$s$ (amu$^{1/2}$ Bohr)", fontsize=9)
    plt.ylabel(r"$\kappa$", fontsize=9)
    plt.axvline(x=xalfa, color='blue', linestyle='-', label=f'Force minimum, x = {xalfa}', linewidth=0.7)
    plt.axvline(x=xgamma, color='red', linestyle='-', label=f'Force maximum, x = {xgamma}', linewidth=0.7)
    #plt.title("Reaction path curvature", fontsize=16)
    #plt.legend(loc="best", fontsize=6)
    plt.tight_layout()
    plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.15)
    plt.show()
    """

    mueff_data = { }
    mueff_data = {'zpe':all_zpe, 'kappa':all_k, 'tbar':all_tbar, 'dtbar':all_dtbar,
                    'mueff':all_mueff, 'mueff_mu':all_mueff_mu}

    return mueff_data



    