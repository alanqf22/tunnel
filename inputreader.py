# -*- coding: utf-8 -*-

## /********************************************************/
## /********************************************************/
##                      INPUTREADER 
## Function: Subroutine that reads the eyringpy input
## file and extracts the information to carry out a
## calculation (method, tunneling model, output file
## names of reactant(s), reactive complex, transition
## state, product complex, product(s), vertical pro-
## ducts; methods for reactions in solution (diffu-
## sion, electron transfer, acid-base), phase of the
## reaction, temperature, symmetry numbers, reaction
## symmetry, pH, pKa, solvent, solvent viscosity, 
## reaction distance and radii of the solutes)    
## /********************************************************/
## /********************************************************/


## Python modules
import os
import sys
from os import path

## Eyringpy modules
from check_gauss_initparam import *

def inputreader(inputfile):

    method = "tst"        # default value
    tunnel = "none"       # default value
    diff = "false"        # default value
    et = "none"           # default value
    ab = "none"           # default value
    phase = "gas"         # default value
    solvent = "none"
    tunn_available = ["wig", "eck", "zct", "sct"]
    diff_available = ["true", "false"]
    et_available = ["marcus"]
    phase_available = ["gas", "sol"]
    freq_files = [ ]
    sp_files = [ ]
    nsym = [ ]
    reacts = [ ]          
    reacts_sp = [ ]       
    rc = [ ]              
    rc_sp = [ ]          
    ts = [ ]              
    ts_sp = [ ]          
    pc = [ ]              
    pc_sp = [ ]           
    prods = [ ]          
    prods_sp = [ ]        
    prodsv = [ ]          
    nsym_reacts = [ ]                                                     
    nsym_rc = [ ]                        
    nsym_ts = [ ]                        
    nsym_pc = [ ]                         
    nsym_prods = [ ]                                                
    rxsym = 1.0           # default value                
    ph = [ ]                             
    visc =  [ ]                          
    rxd = [ ]                             
    rad = [ ]
    pka_r1 = [ ]
    pka_r2 = [ ]
    pka = [ ]                             
    tlist = [ ]                           
    phlist = [ ]                           
    temp = { }
    ph = { }
    temp = {'tlist' : tlist, 'tmin' : "none", 'tmax' : "none", 'tstep' : "none"}    
    ph = {'phlist' : phlist, 'phmin' : "none", 'phmax' : "none", 'phstep' : "none"}

    #IRC-analysis
    pathorientation_available = ["forward", "reverse"]
    rforce_available = ["true", "false"]
    interp_available = ["spline"]
    discrim_available = ["true", "false"]
    xyz_available = ["true", "false"]
    homolumo_available = ["true", "false"]
    dipole_available = ["true", "false"]
    plot_available = ["rforce", "angle", "bondlength", "wiberg", "charge", "homolumo", "dipole", "asm"]
    asm_available = ["true", "false"]
    eda_available = ["true", "false"]

    rforce = "false"                # default value
    interp = "spline"               # default value
    order = int("4")                # default value
    smooth = float("0.0")           # default value
    npoints = int("0")              # default value
    discrim = "false"               # default value
    homolumo = "false"              # default value
    dipole = "false"                # default value

    xyz = "false"                   # default value

    pathpoints = None               # default value
    frag1points = None              # default value
    frag2points = None              # default value
    
    memory = "none"                 # default value
    functional = "none"             # default value
    basis = "none"                  # default value
    pop = "none"                    # default value

    asm = "false"                   # default value
    eda = "false"                   # default value

    pathorientation = "forward"     # default value

    irc_list = [ ]
    scan_list = [ ]
    angle1 = [ ]
    angle2 = [ ]
    angle3 = [ ]
    bondlen1 = [ ]
    bondlen2 = [ ]
    wiberg1 = [ ]
    wiberg2 = [ ]
    charge = [ ]
    plot = [ ]

    frags = [ ]

    rxcoords = [ ]
    rxvtst = [ ]
    rx2remove = [ ]
    degrees = None

    with open(inputfile, 'r') as eif:
        for line in eif:
            lin = line.strip()    
            if len(lin) != 0:   
                if not lin.startswith("#"):    
                    readline = lin.split()    
                    key = readline[0].lower()
                    
                    # Kinetic method
                    if key == "method":    
                        if len(readline) == 2:
                            if str(readline[1]).isalpha():    
                                if str(readline[1]).lower() == "tst":    
                                    method = readline[1].lower()
                                    
                    # Tunneling
                    if key == "tunnel":    
                        if len(readline) == 2:
                            if str(readline[1]).lower() not in tunn_available:
                                print("Error: the available tunneling methods are Wigner (wig), Eckart (eck), zero curvature (zct), and small curvature (sct)!")
                                sys.exit(1)
                            else:
                                tunnel = readline[1].lower()
                                
                    # Diffusion
                    if key == "diff":    
                        if len(readline) == 2:
                            if str(readline[1]).lower() not in diff_available:
                                print("Error: diff can only be true or false!")
                                sys.exit(1)
                            else:
                                diff = readline[1].lower()
                                
                    # Electron transfer
                    if key == "et":    
                        if len(readline) == 2:
                            if str(readline[1]).lower() not in et_available:
                                print("Error: the available electron transfer method is the Marcus theory (marcus)!")
                                sys.exit(1)
                            else:
                                et = readline[1].lower()
                                
                    # Reaction phase
                    if key == "phase":   
                        if len(readline) == 2:
                            if str(readline[1]).lower() not in phase_available:
                                print("Error: the available phases are gas (gas) and solution (sol)!")
                                sys.exit(1)
                            else:
                                phase = readline[1].lower()
                                
                    # Temperature in Kelvin
                    if key == "temp":    
                        if len(readline) <= 3 or len(readline) > 4:
                            # Temperature list
                            for itemp in range(1, len(readline)):
                                if is_number(readline[itemp]):    
                                    if float(readline[itemp]) >= 10.0:    
                                        if float(readline[itemp]) not in tlist:
                                            tlist.append(float(readline[itemp]))    
                                    else:
                                        print("Error: the lowest allowed temperature is 10 Kelvin!")
                                        sys.exit(1)                                  
                                else:
                                    print("Error: the temperature must be a number!")
                                    sys.exit(1)
                        if len(readline) == 4:
                            if is_number(readline[1]) and is_number(readline[2]) and is_number(readline[3]):
                                #Temperature range
                                if float(readline[1]) < float(readline[2]) and float(readline[2]) > float(readline[3]):    
                                    if float(readline[1]) >= 10.0:
                                        tmin = float(readline[1])    # Maximum temperature
                                        tmax = float(readline[2])    # Minimum temperature
                                        tstep = float(readline[3])   # Interval
                                        temp = {'tlist' : tlist, 'tmin' : tmin, 'tmax' : tmax, 'tstep' : tstep}  
                                    else:
                                        print("Error: the lowest allowed temperature is 10 Kelvin!")
                                        sys.exit(1)  
                                # Temperature list
                                else:    
                                    for itemp in range(1, len(readline)):
                                        if float(readline[itemp]) >= 10.0:
                                            if float(readline[itemp]) not in tlist:
                                                tlist.append(float(readline[itemp]))
                                        else:
                                            print("Error: the lowest allowed temperature is 10 Kelvin!")
                                            sys.exit(1)  
                            else:
                                print("Error: the temperature must be a number!")
                                sys.exit(1)  
                                
                    # Specific temperature(s)
                    if key == "stemp":    
                        if len(readline) > 1:     
                            for itemp in range(1, len(readline)):
                                if is_number(readline[itemp]):    
                                    if float(readline[itemp]) >= 10.0:
                                        if float(readline[itemp]) not in tlist:
                                            tlist.append(float(readline[itemp]))    
                                    else:
                                        print("Error: the lowest allowed temperature is 10 Kelvin!")
                                        sys.exit(1)                                 
                                else:
                                    print("Error: the temperature must be a number!")
                                    sys.exit(1)  
                                    
                    # Solvent
                    if key == "solvent":    
                        if len(readline) == 2:
                            if str(readline[1]).isalpha():
                                solvent = readline[1].lower()
                                
                    # Reactant(s) output file
                    if key == "react1":    
                        if len(readline) == 2:
                            reacts.insert(0, readline[1])
                        if len(readline) == 3:
                            reacts.insert(0, readline[1])
                            reacts_sp.insert(0, readline[2])
                    if key == "react2":    
                        if len(readline) == 2:
                            reacts.insert(1, readline[1])
                        if len(readline) == 3:
                            reacts.insert(1, readline[1])
                            reacts_sp.insert(1, readline[2])
                            
                    # Pre-reactive complex output file
                    if key == "rc":    
                        if len(readline) == 2:
                            rc.append(readline[1])
                        if len(readline) == 3:
                            rc.append(readline[1])
                            rc_sp.append(readline[2])
                            
                    # Transition state output file
                    if key == "ts":    
                        if len(readline) == 2:
                            ts.append(readline[1])
                        if len(readline) == 3:
                            ts.append(readline[1])
                            ts_sp.append(readline[2])
                            
                    # Product complex output file
                    if key == "pc":    
                        if len(readline) == 2:
                            pc.append(readline[1])
                        if len(readline) == 3:
                            pc.append(readline[1])
                            pc_sp.append(readline[2])
                            
                    # Product(s) output file
                    if key == "prod1":   
                        if len(readline) == 2:
                            prods.insert(0, readline[1])
                        if len(readline) == 3:
                            prods.insert(0, readline[1])
                            prods_sp.insert(0, readline[2])
                    if key == "prod2":    
                        if len(readline) == 2:
                            prods.insert(1, readline[1])
                        if len(readline) == 3:
                            prods.insert(1, readline[1])
                            prods_sp.insert(1, readline[2])
                            
                    # Vertical product(s) output file
                    if key == "prodv1":    
                         if len(readline) == 2:
                            prodsv.insert(0, readline[1])
                    if key == "prodv2":    
                         if len(readline) == 2:
                            prodsv.insert(1, readline[1])
                            
                    # Reactant(s) symmetry number
                    if key == "nsym_r1":    
                        if len(readline) == 2:
                            if is_number(readline[1]):
                                if float(readline[1]) >= 1.0:
                                    nsym_reacts.insert(0, float(readline[1]))
                                else:
                                    print("Error: the lowest symmetry number is 1!")
                                    sys.exit(1)  
                            else:
                                print("Error: the symmetry number must be a number!")
                                sys.exit(1)  
                    if key == "nsym_r2":    
                         if len(readline) == 2:
                             if is_number(readline[1]):
                                 if float(readline[1]) >= 1.0:
                                     nsym_reacts.insert(1, float(readline[1]))
                                 else:
                                     print("Error: the lowest symmetry number is 1!")
                                     sys.exit(1)  
                             else:
                                 print("Error: the symmetry number must be a number!")
                                 sys.exit(1)  
                                 
                    # Pre-reactive complex symmetry number
                    if key == "nsym_rc":    
                         if len(readline) == 2:
                             if is_number(readline[1]):
                                 if float(readline[1]) >= 1.0:
                                     nsym_rc.append(float(readline[1]))
                                 else:
                                     print("Error: the lowest symmetry number is 1!")
                                     sys.exit(1)  
                             else:
                                 print("Error: the symmetry number must be a number!")
                                 sys.exit(1)  
                            
                    # Transition state symmetry number
                    if key == "nsym_ts":    
                         if len(readline) == 2:
                             if is_number(readline[1]):
                                 if float(readline[1]) >= 1.0:
                                     nsym_ts.append(float(readline[1]))
                                 else:
                                     print("Error: the lowest symmetry number is 1!")
                                     sys.exit(1)  
                             else:
                                 print("Error: the symmetry number must be a number!")
                                 sys.exit(1)  
                                 
                    # Product complex symmetry number
                    if key == "nsym_pc":    
                         if len(readline) == 2:
                             if is_number(readline[1]):
                                 if float(readline[1]) >= 1.0:
                                     nsym_pc.append(float(readline[1]))
                                 else:
                                     print("Error: the lowest symmetry number is 1!")
                                     sys.exit(1)  
                             else:
                                 print("Error: the symmetry number must be a number!")
                                 sys.exit(1)  
                                 
                    # Product(s) symmetry number
                    if key == "nsym_p1":    
                        if len(readline) == 2:
                            if is_number(readline[1]):
                                if float(readline[1]) >= 1.0:
                                    nsym_prods.insert(0, float(readline[1]))
                                else:
                                    print("Error: the lowest symmetry number is 1!")
                                    sys.exit(1)  
                            else:
                                print("Error: the symmetry number must be a number!")
                                sys.exit(1)  
                    if key == "nsym_p2":    
                        if len(readline) == 2:
                            if is_number(readline[1]):
                                if float(readline[1]) >= 1.0:
                                    nsym_prods.insert(1, float(readline[1]))
                                else:
                                    print("Error: the lowest symmetry number is 1!")
                                    sys.exit(1)  
                            else:
                                print("Error: the symmetry number must be a number!")
                                sys.exit(1)  
                                
                    # Reaction symmetry
                    if key == "rxsym":    
                        if len(readline) == 2:
                            if is_number(readline[1]):
                                if float(readline[1]) >= 1.0:
                                    rxsym = float(readline[1])
                                else:
                                    print("Error: the lowest reaction symmetry number is 1!")
                                    sys.exit(1)  
                            else:
                                print("Error: the reaction symmetry must be a number!")
                                sys.exit(1)  
                                
                    # Solvent viscosity in Pa s
                    if key == "visc":    
                        if len(readline) >= 2:
                            for ivisc in range(1, len(readline)):
                                if is_number(readline[ivisc]):
                                    if float(readline[ivisc]) not in visc:
                                        visc.append(float(readline[ivisc]))
                                else:
                                    print("Error: the viscosity must be a number!")
                                    sys.exit(1)  
                                    
                    # Reaction distance in Angstrom
                    if key == "rxd":    
                        if len(readline) == 2:
                            if is_number(readline[1]):
                                if float(readline[1]) > 0:
                                    rxd.append(float(readline[1]))
                                else:
                                    print("Error: the reaction distance must be a number greater than zero!")
                                    sys.exit(1)  
                            else:
                                print("Error: the reaction distance must be a number!")
                                sys.exit(1)  
                                
                    # Radius of the solute(s) in Angstrom
                    if key == "rad1":    
                        if len(readline) == 2:
                            if is_number(readline[1]):
                                if float(readline[1]) > 0:
                                    rad.insert(0, float(readline[1]))
                                else:
                                    print("Error: the radious of the solute must be a number greater than zero!")
                                    sys.exit(1)  
                            else:
                                print("Error: the radius of the solute must be a number!")
                                sys.exit(1)  
                    if key == "rad2":   
                        if len(readline) == 2:
                            if is_number(readline[1]):
                                if float(readline[1]) > 0:
                                    rad.insert(1, float(readline[1]))
                                else:
                                    print("Error: the radious of the solute must be a number greater than zero!")
                                    sys.exit(1)  
                            else:
                                print("Error: the radius of the solute must be a number!")
                                sys.exit(1)  
                                
                    # pKa of reactant(s)
                    if key == "pka_r1":    
                        if len(readline) >=2:
                            for ipka in range(1, len(readline)):
                                if is_number(readline[ipka]):
                                    if float(readline[ipka]) not in pka_r1:
                                        pka_r1.append(float(readline[ipka]))
                                else:
                                    print("Error: the pka must be a number!")
                                    sys.exit(1)  
                        if len(readline) > 4:
                            print("Error: the keyword pka supports up to 3 numerical values!")
                            sys.exit(1)  
                    if key == "pka_r2":    
                        if len(readline) >=2:
                            for ipka in range(1, len(readline)):
                                if is_number(readline[ipka]):

                                    if float(readline[ipka]) not in pka_r2:
                                        pka_r2.append(float(readline[ipka]))
                                else:
                                    print("Error: the pka must be a number!")
                                    sys.exit(1)
                        if len(readline) > 4:
                            print("Error: the keyword pka supports up to 3 numerical values!")
                            sys.exit(1)  
                            
                    # Type of acid-base specie
                    if key == "ab":    
                        if len(readline) == 2:
                            if is_number(readline[1]):
                                if 0 <= int(readline[1]) <= 3:
                                    ab = int(readline[1])
                                else:
                                    print("Error: the keyword ab supports only integers between 0 and 3!")
                                    sys.exit(1)  
                            else:
                                print("Error: the keyword ab supports only one integer value!")
                                sys.exit(1)  
                                
                    # pH of the environment
                    if key == "ph":    
                        if len(readline) <= 3 or len(readline) > 4:
                            # pH list
                            for iph in range(1, len(readline)):
                                if is_number(readline[iph]):
                                    if 0 <= float(readline[iph]) <= 14.0:
                                        if float(readline[iph]) not in phlist:
                                            phlist.append(float(readline[iph]))
                                    else:
                                        print("Error: the pH must be between 0 and 14!")
                                        sys.exit(1)  
                        if len(readline) == 4:
                            if is_number(readline[1]) and is_number(readline[2]) and is_number(readline[3]):
                                if float(readline[1]) < float(readline[2]) and float(readline[2]) > float(readline[3]):
                                    #pH range
                                    if float(readline[1]) >= 0.0:
                                        if float(readline[2]) <= 14.0:
                                            phmin = float(readline[1])    # Maximum pH
                                            phmax = float(readline[2])   # Minimum pH
                                            phstep = float(readline[3])    # Interval
                                            ph = {'phlist' : phlist, 'phmin' : phmin, 'phmax' : phmax, 'phstep' : phstep}
                                        else:
                                            print("Error: the highest allowed pH is 14!")
                                            sys.exit(1)  
                                    else:
                                        print("Error: the lowest allowed pH is 0!")
                                        sys.exit(1)       
                                else:
                                    # pH list
                                    for iph in range(1, len(readline)):
                                        if 0 <= float(readline[iph]) <= 14.0:
                                            if float(readline[iph]) not in phlist:
                                                phlist.append(float(readline[iph]))
                                        else:
                                            print("Error: the pH must be between 0 and 14!")
                                            sys.exit(1)  
                            else:
                                print("Error: the pH must be a number!")
                                sys.exit(1)  

                    # Specific pH
                    if key == "sph":    
                        if len(readline) > 1:     
                            for iph in range(1, len(readline)):
                                if is_number(readline[iph]):
                                    if 0 <= float(readline[iph]) <= 14.0:
                                        if float(readline[iph]) not in phlist:
                                            phlist.append(float(readline[iph]))
                                    else:
                                        print("Error: the pH must be between 0 and 14!")
                                        sys.exit(1)  
                                else:
                                    print("Error: the pH must be a number!")
                                    sys.exit(1)

                    ######################################## IRC-Analysis ########################################
                    # IRC(s) output file
                    if key == "irc":
                        if len(readline) == 2:
                            irc_list.insert(0, readline[1])
                        if len(readline) == 3:
                            irc_list.insert(0, readline[1])
                            irc_list.insert(1, readline[2])
                        if len(readline) == 1:
                            print("IRC output(s) must be declared!")
                            sys.exit(1)
                        if len(readline) > 3:
                            print("IRC output(s) must be one or two!")
                            sys.exit(1)
                    
                    # Path orientation
                    if key == "pathorientation":
                        if len(readline) == 2:
                            if str(readline[1]).lower() == "forward":
                                pathorientation = "forward"
                            if str(readline[1]).lower() == "reverse":
                                pathorientation = "reverse"
                            if str(readline[1]).lower() not in pathorientation_available:
                                print("Error: the only available options for PATHORIENTATION are forward or reverse")
                                sys.exit(1)
                        if len(readline) == 1:
                            pathorientation = "forward"
                        if len(readline) > 2:
                            print("Error: PATHORIENTATION must have only one value!")
                            sys.exit(1)

                    ## ------------------- Keywords related to the reaction force analysis -------------------- ##
                    # Reaction Force Analysis
                    if key == "rforce":
                        if len(readline) == 2:
                            if str(readline[1]).lower() == "true":
                                rforce = "true"
                            if str(readline[1]).lower() == "false":
                                rforce == "false"
                            if str(readline[1]).lower() not in rforce_available:
                                print("Error: the only available options for RFORCE are true or false")
                                sys.exit(1)
                                
                        if len(readline) == 1:
                            rforce = "false"
                            
                        if len(readline) > 2:
                            print("Error: RFORCE must have only one value!")
                            sys.exit(1)
                                    
                    # Interpolation Method
                    if key == "interp":    
                        if len(readline) == 2:
                            if str(readline[1]).isalpha():
                                interp = readline[1].lower()
                                if str(readline[1]).lower() not in interp_available:
                                    print("Error: the only interpolation method available is spline!")
                                    sys.exit(1)
                                else:
                                    interp = readline[1]
                        if len(readline) == 1:
                            interp = "spline"
                        if len(readline) >= 3:
                            print("Error: the INTERP entry must be just one!")
                            sys.exit(1)
                                    
                    # Order
                    if key == "order":
                        if len(readline) == 2:
                            if is_number(readline[1]):
                                if 3 <= int(readline[1]) <= 5:
                                    order = int(readline[1])
                                else:
                                    print("Error: the interpolation order only can take integer values between 3 and 5!")
                                    sys.exit(1)
                            else:
                                print("Error: the interpolation order must be a number!")
                                sys.exit(1)  
                    
                    # Smooth
                    if key == "smooth":    
                        if len(readline) == 2:
                            if is_number(readline[1]):
                                if float(readline[1]) >= 0.0:
                                    smooth = float(readline[1])
                                else:
                                    print("Error: the lowest smooth factor is 0.0!")
                                    sys.exit(1)  
                            else:
                                print("Error: smooth factor must be a number!")
                                sys.exit(1)
                    
                    # New points
                    if key == "npoints":    
                        if len(readline) == 2:
                            if is_number(readline[1]):
                                if int(readline[1]) >= 0:
                                    npoints = int(readline[1])
                                else:
                                    print("Error: Number of new points can't be negative or a float number!")
                                    sys.exit(1)  
                            else:
                                print("Error: NPOINTS must be a number!")
                                sys.exit(1)

                    # Discrimination
                    if key == "discrim":    
                        if len(readline) == 2:
                            if str(readline[1]).isalpha():
                                discrim = readline[1].lower()
                                if str(readline[1]).lower() not in discrim_available:
                                    print("Error: the only available options for DISCRIM are true or false!")
                                    sys.exit(1)
                                else:
                                    discrim = readline[1]
                        if len(readline) == 1:
                            discrim = "false"
                        if len(readline) >= 3:
                            print("Error: the DISCRIM entry must be just one!")
                            sys.exit(1)

                    ## --------------------- Keywords related to the properties analysis ---------------------- ##
                    # Homo lumo
                    if key == "homolumo":
                        if len(readline) == 2:
                            if str(readline[1]).lower() == "true":
                                homolumo = "true"
                            if str(readline[1]).lower() == "false":
                                homolumo == "false"
                            if str(readline[1]).lower() not in homolumo_available:
                                print("Error: the only available options for HOMOLUMO are true or false!")
                                sys.exit(1)
                        if len(readline) == 1:
                            homolumo = "false"
                        if len(readline) > 2:
                            print("Error: HOMOLUMO must have only one value!")
                            sys.exit(1)
                    
                    # Dipole moment
                    if key == "dipole":
                        if len(readline) == 2:
                            if str(readline[1]).lower() == "true":
                                dipole = "true"
                            if str(readline[1]).lower() == "false":
                                dipole == "false"
                            if str(readline[1]).lower() not in dipole_available:
                                print("Error: the only available options for DIPOLE are true or false!")
                                sys.exit(1)
                        if len(readline) == 1:
                            dipole = "false"
                        if len(readline) > 2:
                            print("Error: DIPOLE must have only one value!")
                            sys.exit(1)
                    
                    # Charge
                    if key == "charge":
                        if len(readline) >= 2:
                            for i in range(1, len(readline)):
                                if is_number(readline[i]):
                                    charge.append(int(readline[i]))
                        if len(readline) < 2:
                            print("CHARGE argument must be at least one!")
                            sys.exit(1)
                    
                    # Bond length
                    if key == "bondlength":
                        bondlen = [ ]
                        for i in range(1, len(readline)):
                            bondlen.append(readline[i].split('-'))
                        for m in bondlen:
                            if is_number(m[0]):
                                bondlen1.append(int(m[0]))
                            if is_number(m[1]):
                                bondlen2.append(int(m[1]))
                            if len(bondlen1) != len(bondlen2):
                                print("Values in BONDLEN must be provided in pairs!")
                                sys.exit(1)
                    
                    # Wiberg
                    if key == "wiberg":
                        wiberg = [ ]
                        for i in range(1, len(readline)):
                            wiberg.append(readline[i].split('-'))
                        for m in wiberg:
                            if is_number(m[0]):
                                wiberg1.append(int(m[0]))
                            if is_number(m[1]):
                                wiberg2.append(int(m[1]))
                            if len(wiberg1) != len(wiberg2):
                                print("Values in WIBERG must be provided in pairs!")
                                sys.exit(1)
                    
                    # Angle
                    if key == "angle":
                        angle = [ ]
                        for i in range(1, len(readline)):
                            angle.append(readline[i].split('-'))
                        for m in angle:
                            if is_number(m[0]):
                                angle1.append(int(m[0]))
                            if is_number(m[1]):
                                angle2.append(int(m[1]))
                            if is_number(m[1]):
                                angle3.append(int(m[2]))
                            if len(angle1) != len(angle2)or len(angle1) != len(angle3) or len(angle2) != len(angle3):
                                print("Values in ANGLE must be provided in trios!")
                                sys.exit(1)
                    
                    # xyz file
                    if key == "xyz":
                        if len(readline) == 2:
                            if str(readline[1]) == "true":
                                xyz = "true"
                            if str(readline[1]) == "false":
                                xyz == "false"
                            if str(readline[1]) not in xyz_available:
                                print("Error: the only available options for XYZ are true or false")
                        if len(readline) == 1:
                            xyz = "false"
                        if len(readline) > 2:
                            print("Error: XYZ must have only one value!")
                            sys.exit(1)                    
                                
                    # Gaussian's NBO input files specifications
                    if key == "nboinput":
                        if len(readline) == 5:
                            if str(readline[1]):
                                memory = str(readline[1])
                            if str(readline[2]):
                                functional = str(readline[2])
                            if str(readline[3]):
                                basis = str(readline[3])
                            if str(readline[4]):
                                pop = str(readline[4])
                                
                        if len(readline) < 5:
                            print("Error: NBOINPUT must have four entry values: memory, functional, basis and POP keyword specification!")
                            sys.exit(1)
                    
                    # Plots
                    if key == "plot":
                        if len(readline) >= 2:
                            for i in range(1, len(readline)):
                                if str(readline[i]).isalpha():
                                    if str(readline[i]).lower() in plot_available:
                                        plot.append(str(readline[i].lower()))
                                    else:
                                        print(readline[i])
                                        print("is not a valid argument in the PLOT keyword")
                                        sys.exit(1)
                                        
                        if len(readline) < 2:
                            print("PLOT argument must be at least one!")
                            sys.exit(1)


                    ## ---------------- Keywords related to the Distortion/Interaction - Activation Strain Model ----------------- ##
                    # Activation Strain Model
                    if key == "asm":
                        if len(readline) == 2:
                            if str(readline[1]).lower() == "true":
                                asm = "true"
                            if str(readline[1]).lower() == "false":
                                asm == "false"
                            if str(readline[1]).lower() not in asm_available:
                                print("Error: the only available options for ASM are true or false")
                                sys.exit(1)
                        if len(readline) == 1:
                            asm = "false"
                        if len(readline) > 2:
                            print("Error: ASM must have only one value!")
                            sys.exit(1)

                    # Energy Decomposition Analysis
                    if key == "eda":
                        if len(readline) == 2:
                            if str(readline[1]).lower() == "true":
                                eda = "true"
                            if str(readline[1]).lower() == "false":
                                eda == "false"
                            if str(readline[1]).lower() not in eda_available:
                                print("Error: the only available options for EDA are true or false")
                                sys.exit(1)
                        if len(readline) == 1:
                            eda = "false"
                        if len(readline) > 2:
                            print("Error: EDA must have only one value!")
                            sys.exit(1)

                    # Scan output file
                    if key == "scan":
                        if len(readline) == 2:
                            scan_list.insert(0, readline[1])
                        if len(readline) == 1:
                            print("Scan file must be declared!")
                            sys.exit(1)
                        if len(readline) > 3:
                            print("Scan file must be one!")
                            sys.exit(1)

                    # Directories names
                    if key == "pathpoints":
                        if len(readline) == 2:
                            if str(readline[1]):
                                pathpoints = str(readline[1])
                        else:
                            print("Error: PATHPOINTS must have only one value with the folder name that contains the single point or NBO files!")
                            sys.exit(1)

                    if key == "frag1points":
                        if len(readline) == 2:
                            if str(readline[1]):
                                frag1points = str(readline[1])
                        else:
                            print("Error: FRAG1POINTS must have only one value with the folder name that contains the single point files!")
                            sys.exit(1)

                    if key == "frag2points":
                        if len(readline) == 2:
                            if str(readline[1]):
                                frag2points = str(readline[1])
                        else:
                            print("Error: FRAG2POINTS must have only one value with the folder name that contains the single point files!")
                            sys.exit(1)


                    # Reaction Coordinates
                    if key == "rxcoords":    
                        if len(readline) >= 5:
                            for irx in range(1, len(readline)):
                                if is_number(readline[irx]):    
                                    rxcoords.append(float(readline[irx]))                                     
                                else:
                                    print("Error: " + readline[irx] + " is not a number!")
                                    sys.exit(1)
                        else:
                            print("Error: at least 5 coordinates are needed to reproduce the potential energy curve.")
                            sys.exit(1)

                    # Reaction Coordinates of the Variational Transition States 
                    if key == "rxvtst":    
                        if len(readline) >= 2:
                            for irx in range(1, len(readline)):
                                if is_number(readline[irx]):    
                                    rxvtst.append(float(readline[irx]))                                     
                                else:
                                    print("Error: " + readline[irx] + " is not a number!")
                                    sys.exit(1)

                    # Reaction Coordinates to remove 
                    if key == "rx2remove":    
                        if len(readline) >= 2:
                            for irx in range(1, len(readline)):
                                if is_number(readline[irx]):    
                                    rx2remove.append(float(readline[irx]))                                     
                                else:
                                    print("Error: " + readline[irx] + " is not a number!")
                                    sys.exit(1)

                    # Degrees of freedom
                    if key == "degrees":
                        if len(readline) == 2:
                            if is_number(readline[1]):
                                if int(readline[1]) > 1:
                                    degrees = int(readline[1])
                                else:
                                    print("Error: degrees of freedom must be a number greater than or equal to 1.")
                                    sys.exit(1)
                            else:
                                print("Error: degrees of freedom must be a number!")
                                sys.exit(1)
                    
        if len(reacts) != 0 and len(prods) != 0:
            for ireact in reacts:
                freq_files.append(ireact)
            for iprod in prods:
                freq_files.append(iprod)
            if len(ts) == 1:
                freq_files.append(ts[0])
            if len(rc) == 1:
                freq_files.append(rc[0])
            if len(pc) == 1:
                freq_files.append(pc[0])

            if len(reacts_sp) > 0:
                for ireact_sp in reacts_sp:
                    sp_files.append(ireact_sp)
            if len(prods_sp) > 0:
                for iprod_sp in prods_sp:
                    sp_files.append(iprod_sp)
            if len(ts_sp) == 1:
                sp_files.append(ts_sp[0])
            if len(rc_sp) == 1:
                sp_files.append(rc_sp[0])
            if len(pc_sp) == 1:
                sp_files.append(pc_sp[0])

            if len(nsym_reacts) > 0 and len(nsym_prods) > 0:
                for isym_r in nsym_reacts:
                    nsym.append(isym_r)
                for isym_p in nsym_prods:
                    nsym.append(isym_p)
                if len(ts) == 1:
                    if len(nsym_ts) == 1:
                        nsym.append(nsym_ts[0])
                if len(rc) == 1:
                    if len(nsym_rc) == 1:
                        nsym.append(nsym_rc[0])
                if len(pc) == 1:
                    if len(nsym_pc) == 1:
                        nsym.append(nsym_pc[0])

        if len(pka_r1) > 0:
            pka.append(sorted(pka_r1, key=float, reverse=True))
        if len(pka_r2) > 0:
            pka.append(sorted(pka_r2, key=float, reverse=True))

                                                   
    inputdata = { }
    inputdata = {'method' : method, 'tunnel' : tunnel, 'diff' : diff,
                 'et' : et, 'ab' : ab, 'phase' : phase, 'solvent' : solvent,
                 'freq_files' : freq_files, 'sp_files' : sp_files,
                 'nsym' : nsym, 'temperature' : temp, 'reacts' : reacts,
                 'reacts_sp' : reacts_sp, 'rc' : rc, 'rc_sp' : rc_sp,
                 'ts' : ts, 'ts_sp' : ts_sp, 'pc' : pc, 'pc_sp' : pc_sp,
                 'prods' : prods, 'prods_sp' : prods_sp, 'prodsv' : prodsv,
                 'nsym_reacts' : nsym_reacts, 'nsym_rc' : nsym_rc,
                 'nsym_ts' : nsym_ts, 'nsym_pc' : nsym_pc,
                 'nsym_prods' : nsym_prods, 'rxsym' : rxsym, 'visc' : visc,
                 'rxd' : rxd, 'rad' : rad, 'ph' : ph, 'pka' : pka, 'rforce' : rforce,
                 'interp' : interp, 'order' : order, 'smooth' : smooth,
                 'npoints' : npoints, 'irc_list' : irc_list, 'xyz' : xyz,
                 'angle1' : angle1, 'angle2' : angle2, 'angle3' : angle3, 
                 'bondlen1' : bondlen1, 'bondlen2' : bondlen2, 'wiberg1' : wiberg1,
                 'wiberg2': wiberg2, 'charge' : charge, 'discrim' : discrim,
                 'homolumo' : homolumo, 'dipole' : dipole, 'plot' : plot,
                 'memory' : memory, 'functional' : functional, 'basis' : basis,
                 'pop' : pop, 'pathpoints' : pathpoints, 'frag1points' : frag1points,
                 'frag2points' : frag2points, 'asm' : asm, 'eda' : eda,'pathorientation' : pathorientation,
                 'scan_list' : scan_list, 'rxcoords' : rxcoords, 'rxvtst' : rxvtst,
                 'rx2remove' : rx2remove, 'degrees' : degrees}
    
    return inputdata
