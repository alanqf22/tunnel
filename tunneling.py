#!/Library/Frameworks/Python.framework/Versions/3.6/bin/python3
# -*- coding: utf-8 -*-

##Python modules
import os
import sys
import time
import datetime
from os import path
from decimal import Decimal, getcontext
from scipy.interpolate import UnivariateSpline

##Eyringpy modules
from SCT import *
from ZCT import *
from funcs import *
from spline import *
from constants import *
from inputreader import *
from sort_by_rxcoord import *
from gaussian_parser import *

## /********************************************************/
## /********************************************************/
##                         TUNNELING 
## Function: Calculates multidimensional tunneling
##           transmission coefficient.
## /********************************************************/
## /********************************************************/

def VGa(V,zpe):
    """
    Calculate vibrationally adiabatic ground-state potential energy curve (VGa).

    Parameters:
    V   : array_like
          Potential energies.
    zpe : array_like
          Zero-point energies.

    Returns:
    VGa : array_like
          The result of summing V + zpe. 
    """
    VGa = [ ]
    for i in range(len(V)):
        VGa.append(V[i]+zpe[i])
    return VGa

def VGa_critical_points(s, VGa):
    """
    Get the critical points of the VGa curve.

    Parameters:
    s : array_like
        Reaction coordinates.
    VGa : array_like
        Vibrationally adiabatic ground-state potential energies.

    Returns:
    E0   : Maximum between the VGa of reactants or the VGa of products.
    E0s  : Reaction coordinate of E0.
    VAG  : Maximum value of the VGa curve.
    VAGs : Reaction coordinate of VGA.
    """
    E0 = max(VGa[0], VGa[-1])
    E0s = s[VGa.index(E0)]
    VAG = max(VGa)
    VAGs = s[VGa.index(VAG)]

    VGa_critical_points = { }
    VGa_critical_points['E0'] = E0
    VGa_critical_points['E0_s'] = E0s
    VGa_critical_points['VAG'] = VAG
    VGa_critical_points['VAG_s'] = VAGs
    return VGa_critical_points

if __name__ == "__main__":
    start_time = time.time()
    basename = os.path.basename(sys.argv[0])
    tunneling = os.path.splitext(basename)[0]
    username = os.path.expanduser('~').split(os.sep)[-1]
    dirname = os.path.basename(os.getcwd())
    now = datetime.datetime.now()
    date = now.strftime("%Y-%m-%d %H:%M:%S")
    inputfile = dirname + ".eif"                    #Input file

    #Data from input reader:
    inpdata = inputreader(inputfile)                #Data from *.eif
    reacts_list = inpdata['reacts']                 #Reactant(s)
    prods_list = inpdata['prods']                   #Product(s)
    irc_list = inpdata['irc_list']                  #IRC Gaussian output(s)
    pathpoints_dirname = inpdata['pathpoints']      #IRC path points: *.fchk files
    tunnel_method = inpdata['tunnel']               #Tunneling model
    temperatures = inpdata['temperature']['tlist']  #Temperature(s) in Kelvin
    rxcoords = inpdata['rxcoords']
    rxvtst = inpdata['rxvtst']
    range_remove = inpdata['rx2remove']
    degrees = inpdata['degrees']
    
    inverse = "false"

    #Verifying input data
    #VTS reaction coordinates and temperatures
    if len(rxvtst) == 0:
        print("Erorr: Reaction coordinates of Variational Transition States must be specified.")

    if len(temperatures) == 0:
        print("Erorr: Temperatures must be specified.")

    if len(rxvtst) != len(temperatures):
        print("Erorr: Number of temperatures must be equal to the number")
        print("       of coordinates of the variational transition states.")
        sys.exit(1)

    #Saving *.log/.out files for ZCT or *.fchk files for SCT
    pathpoints_list = []
    if tunnel_method == "zct":
        if pathpoints_dirname:
            pathpoints_files = glob.glob(f"{pathpoints_dirname}/*.out") or glob.glob(f"{pathpoints_dirname}/*.log")
            pathpoints_list.extend(pathpoints_files)
        else:
            print("Error: PATHPOINTS keyword must be specified to perform ZCT model.")
            sys.exit(1)

    if tunnel_method == "sct":
        if pathpoints_dirname:
            pathpoints_files = glob.glob(f"{pathpoints_dirname}/*.fchk")
            pathpoints_list.extend(pathpoints_files)
        else:
            print("Error: PATHPOINTS keyword must be specified to perform SCT model.")
            sys.exit(1)

    if len(pathpoints_list) == 0:
        print("Error: No files in the " + str(pathpoints_dirname) + " folder!")
        sys.exit(1)
    else:
        pathpoints_list.sort()

    #Reaction Coordinates
    all_rxcoord = [ ]
    if len(rxcoords) > 0:
        all_rxcoord = rxcoords
    else:
        if len(irc_list) > 0: 
            irc_data = sort_by_rxcoord(irc_list)
            coords = irc_data['rxcoord']
            if len(pathpoints_list) == coords:
                all_rxcoord = coords
            else:
                names_without_extension = [os.path.splitext(os.path.basename(file))[0] for file in pathpoints_list]
                rx_positions = [irx+1 for irx in range(len(coords))]
                for icoord in rx_positions:
                    for ifile in names_without_extension:
                        if int(icoord) == int(ifile):
                            all_rxcoord.append(coords[icoord-1])
        else:
            print("Erorr: Reaction coordinates must be specified!")
            sys.exit(1)

    if inverse == "true":
        pathpoints_list.sort(reverse=True)
        all_rxcoord = [-x for x in reversed(all_rxcoord)]

    if tunnel_method == "sct":
        #Input data to calculate effective reduced mass
        all_Gcc = get_cartesian_gradient(pathpoints_list)           #Energy gradient in cartesian coordinates
        all_Fcc = get_cartesian_force_constants(pathpoints_list)    #Force constants matrix in cartesian coordinates
        all_ccords = get_cartesian_coords(pathpoints_list)          #Cartesian coordinates
        all_masses = get_atomic_weights(pathpoints_list)            #Atomic weights

        #Calculating effective reduced mass
        mueff_data = mu_eff(all_Gcc, all_Fcc, all_ccords, all_masses, all_rxcoord, degrees)
        all_k = mueff_data['kappa']                                 #Reaction path curvature
        all_mu = mueff_data['mueff']                                #Effective reduced mass
        
        #Input data to calculate Vibratinally-adiabatic ground-state potential energy
        all_V = get_electronic_energy_fchk(pathpoints_list)         #Electronic energy in HA
        all_zpe = mueff_data['zpe']                                 #Zero-point energy in HA

    if tunnel_method == "zct":
        all_V = get_electronic_energy(pathpoints_list)              #Electronic energy in HA
        all_zpe = get_zpe(pathpoints_list)                          #Zero-point energy in HA
        all_mu = [MU] * len(all_V)

    #Vibrationally adiabatic ground-state potential energy
    all_Ein = all_V
    all_Ein_kcal = [iE * HARTREES_TO_KCAL_MOL for iE in all_Ein]

    all_V = [iV - max(all_V[0], all_V[-1]) for iV in all_V]                         
    all_VGa = VGa(all_V,all_zpe)                                    #VGa(s) = V(s) + zpe(s) [in HA]

    all_V_kcal = [iV * HARTREES_TO_KCAL_MOL for iV in all_V]        #Electronic energy in HA
    all_zpe_kcal = [iZ * HARTREES_TO_KCAL_MOL for iZ in all_zpe]    #Zero-point energy in HA
    all_VGa_kcal = [i*HARTREES_TO_KCAL_MOL for i in all_VGa]        #VGa(s) in kcal mol-1

    VGadata =  VGa_critical_points(all_rxcoord, all_VGa)
    E0 = VGadata['E0']
    VAG = VGadata['VAG']
    sAG = VGadata['VAG_s']

    print("VAG:")
    print(sAG)
    print("")

    # Print VGa and mueff/mu for initial data
    eof = "VGa_mu_initial.eof"
    with open(eof, 'w') as eof_file:
        eof_file.write(f"{'s':<15}{'E':<18}{'V(s)-Eref':<18}{'zpe(s)':<18}{'VaG(s)':<18}{'mueff/mu':<18}\n")
        eof_file.write(f"{'':<15}{'[kcal mol⁻¹]':<18}{'[kcal mol⁻¹]':<18}{'[kcal mol⁻¹]':<18}{'[kcal mol⁻¹]':<18}{'':<18}\n")
        eof_file.write("-" * 85 + "\n")
        for i in range(len(all_VGa_kcal)):
            eof_file.write(f"{all_rxcoord[i]:<15.5f}{all_Ein_kcal[i]:<18.2f}{all_V_kcal[i]:<18.2f}{all_zpe_kcal[i]:<18.2f}{all_VGa_kcal[i]:<18.2f}{(all_mu[i]/MU):<18.5E}\n")
        eof_file.write("                                    \n")

    all_rxcoord = all_rxcoord[:-1]
    all_VGa_kcal = all_VGa_kcal[:-1]
    all_V = all_V[:-1]
    all_zpe = all_zpe[:-1]
    all_VGa = all_VGa[:-1]
    all_mu = all_mu[:-1]   
    
    # Plot VaG(s) profile
    # Force minima and maxima for various regions
    regions = {
    "heavy_R1": (-1.7995, 1.6995),
    "heavy_R2": (-1.4996, 2.2994),
    "heavy_R3": (-1.7996, 1.7996),
    "heavy_R4": (-1.7996, 2.8994),
    "heavy_R5": (-1.5996, 2.6994),
    "heavy_R6": (-1.5996, 2.7994),
    "heavy_R7": (-1.5997, 1.7996),
    "heavy_R8": (-1.1997,2.0995),
    "heavy_R9": (-1.5996,2.0995),
    "heavy_R10": (-1.0998, 1.0997),
    "heavy_R11": (-0.99974, 1.99953),
    "heavy_R12": (-1.3996,2.2994),
    "heavy_R13": (-1.29969, 3.09930)
    }

    xalfa, xgamma = regions["heavy_R12"]  # Set active region
    plt.figure(figsize=(8, 6), dpi=300) #Size
    plt.plot(all_rxcoord, all_VGa_kcal, 'ko', label=r"$V_a^G$", markersize=0.8)
    plt.xlabel(r"$s$ (amu$^{1/2}$ Bohr)", fontsize=9)
    plt.ylabel(r"$V_a^G$ (kcal mol$^{-1}$)", fontsize=9)
    plt.axvline(x=xalfa, color='blue', linestyle='-', label=f'Force minimum, x = {xalfa}', linewidth=0.7)
    plt.axvline(x=xgamma, color='red', linestyle='-', label=f'Force maximum, x = {xgamma}', linewidth=0.7)
    #plt.title("Perfil de Energía - VGa", fontsize=16)
    #plt.legend(loc="best", fontsize=6)
    plt.tight_layout()
    plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.15)
    plt.show()

    column_widths = {'T': 12,
    'Tunneling': 15,
    '%Tunneling': 10,
    'Reflection': 15,
    '%Reflection': 10,
    'Kappa': 15,
    'Kappa - CAG': 10}

    ########## Small-Curvature tunneling data and output file ##########
    with open("kappa_sct.eof", "w") as file1:
        # Write the table header with aligned columns
        file1.write(f"{'T [K]':<{column_widths['T']}}{'Tunneling':<{column_widths['Tunneling']}}"
                   f"{'%Tun':<{column_widths['%Tunneling']}}{'Reflection':<{column_widths['Reflection']}}"
                   f"{'%Ref':<{column_widths['%Reflection']}}{'Kappa-SCT':<{column_widths['Kappa']}}"
                   f"{'Kappa-CAG':<{column_widths['Kappa - CAG']}}\n")
        file1.write("-" * 90 + "\n")  # Optional: adds a line to separate header from data

        for iT in range(len(temperatures)):
            #Calculate SCT kappa 
            sct_data = kappa(all_rxcoord, all_VGa, E0, VAG, rxvtst[iT], temperatures[iT], all_mu)

            # Convert relevant values from numpy.float128 to Decimal safely
            sct_ktun = Decimal(str(sct_data['kappatun']))
            sct_kref = Decimal(str(sct_data['kapparef']))
            sct_kappa = Decimal(str(sct_data['kappa']))
            sct_kcvt_cag = sct_data['kcvt_cag']      
            sct_Etun_kcal = sct_data['Etun']
            sct_Etunprob = sct_data['Etunprob']
            sct_Etuncont = sct_data['Etuncont']     

            # Calculate the percentage of SCT Tunneling and Reflection
            sct_pct_tun = (sct_ktun / sct_kappa) * Decimal(100)
            sct_pct_ref = (sct_kref / sct_kappa) * Decimal(100)

            # Write SCT kappa data in table format
            file1.write(f"{temperatures[iT]:<{column_widths['T']}.1f}"
                        f"{sct_ktun:<{column_widths['Tunneling']}.5E}"
                        f"{sct_pct_tun:<{column_widths['%Tunneling']}.1f}"
                        f"{sct_kref:<{column_widths['Reflection']}.5E}"
                        f"{sct_pct_ref:<{column_widths['%Reflection']}.1f}"
                        f"{sct_kappa:<{column_widths['Kappa']}.5E}"
                        f"{sct_kcvt_cag:<{column_widths['Kappa - CAG']}.1f}\n")
            file1.write("                                    \n")

            # Write the contributions of each energy to SCT kappa
            with open("tunprobs_" + str(temperatures[iT]) + "_sct.eof", "w") as file2:
                # Write the header for the columns, adding units for energy
                file2.write("{:<20} {:<15} {:<15}\n".format("Etun [kcal mol⁻¹]", "TunProb", "Cont [%]"))
                file2.write("-" * 70 + "\n") 
                for i in range(len(sct_Etun_kcal)):
                    file2.write("{:<20.4f} {:<15.2E} {:<15.2E}\n".format(sct_Etun_kcal[i], sct_Etunprob[i], sct_Etuncont[i]))
                file2.write("                                    \n")

    ########## Zero-Curvature tunneling data and output file ##########
    with open("kappa_zct.eof", "w") as file:
        # Write the table header with aligned columns
        file.write(f"{'T [K]':<{column_widths['T']}}{'Tunneling':<{column_widths['Tunneling']}}"
                   f"{'%Tun':<{column_widths['%Tunneling']}}{'Reflection':<{column_widths['Reflection']}}"
                   f"{'%Ref':<{column_widths['%Reflection']}}{'Kappa-ZCT':<{column_widths['Kappa']}}"
                   f"{'Kappa-CAG':<{column_widths['Kappa - CAG']}}\n")
        file.write("-" * 90 + "\n")  # Optional: adds a line to separate header from data

        for iT in range(len(temperatures)):
            zct_data = kappa(all_rxcoord, all_VGa, E0, VAG, rxvtst[iT], temperatures[iT])

            # Convert relevant values from numpy.float128 to Decimal safely
            zct_ktun = Decimal(str(zct_data['kappatun']))
            zct_kref = Decimal(str(zct_data['kapparef']))
            zct_kappa = Decimal(str(zct_data['kappa']))
            zct_kcvt_cag = zct_data['kcvt_cag']      
            zct_Etun_kcal = zct_data['Etun']
            zct_Etunprob = zct_data['Etunprob']
            zct_Etuncont = zct_data['Etuncont']                 

            # Calculate the percentage of Tunneling and Reflection using Decimal
            zct_pct_tun = (zct_ktun / zct_kappa) * Decimal(100)
            zct_pct_ref = (zct_kref / zct_kappa) * Decimal(100)

            # Write the values in table format with aligned columns
            file.write(f"{temperatures[iT]:<{column_widths['T']}.2f}"
                        f"{zct_ktun:<{column_widths['Tunneling']}.5E}"
                        f"{zct_pct_tun:<{column_widths['%Tunneling']}.1f}"
                        f"{zct_kref:<{column_widths['Reflection']}.5E}"
                        f"{zct_pct_ref:<{column_widths['%Reflection']}.1f}"
                        f"{zct_kappa:<{column_widths['Kappa']}.5E}"
                        f"{zct_kcvt_cag:<{column_widths['Kappa - CAG']}.1f}\n")
            file.write("                                    \n")

            # Write the contributions of each energy to SCT kappa
            with open("tunprobs_" + str(temperatures[iT]) + "_zct.eof", "w") as file2:
                # Write the header for the columns, adding units for energy
                file2.write("{:<20} {:<15} {:<15}\n".format("Etun [kcal mol⁻¹]", "TunProb", "Cont [%]"))
                file2.write("-" * 70 + "\n") 
                for i in range(len(zct_Etun_kcal)):
                    file2.write("{:<20.4f} {:<15.2E} {:<15.2E}\n".format(zct_Etun_kcal[i], zct_Etunprob[i], zct_Etuncont[i]))
                file2.write("                                    \n")

    all_rxcoord_final = sct_data['scoord']
    all_VGa_final_kcal = sct_data['VGa']
    mueff_mu_final = sct_data['mueff_mu']
    
    # Print interpolated VGa with mueff(s)/mu
    eof = "VGa_mu_interp.eof"
    with open(eof, 'w') as eof_file:
        eof_file.write(f"{'coordinate':<20}{'VGa':<18}{'mueff/mu':<18}\n")
        eof_file.write(f"{'':<20}{'[kcal mol⁻¹]':<18}{'':<18}\n")
        eof_file.write("-" * 55 + "\n")
        for i in range(len(all_rxcoord_final)):
            eof_file.write(f"{all_rxcoord_final[i]:<20.5f}{all_VGa_final_kcal[i]:<18.2f}{mueff_mu_final[i]:<18.5E}\n")
        eof_file.write("                                    \n")

    # Plot mueff(s)/mu profile
    plt.figure(figsize=(8, 6), dpi=300) #Size
    plt.plot(all_rxcoord_final, mueff_mu_final, 'k-', label=r"$mu^{eff}/mu$", markersize=0.8)
    plt.xlabel(r"$s$ (amu$^{1/2}$ Bohr)", fontsize=9)
    plt.ylabel(r"$\mu^{eff}/\mu$", fontsize=9)
    plt.axvline(x=xalfa, color='blue', linestyle='-', label=f'Force minimum, x = {xalfa}', linewidth=0.7)
    plt.axvline(x=xgamma, color='red', linestyle='-', label=f'Force maximum, x = {xgamma}', linewidth=0.7)
    plt.tight_layout()
    plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.15)
    plt.show()

    #Plot Tunneling Probabilities
    E0toVAG_energies = [ ]
    E0toVAG_scoords = [ ]
    for irx,iE in zip(all_rxcoord_final,all_VGa_final_kcal):
        if iE >= E0 * HARTREES_TO_KCAL_MOL:
            E0toVAG_scoords.append(irx)
            E0toVAG_energies.append(iE)

    plot_color_map(E0toVAG_scoords,E0toVAG_energies,sct_Etunprob)
    plot_color_map(E0toVAG_scoords,E0toVAG_energies,zct_Etunprob)

    