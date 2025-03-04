# -*- coding: utf-8 -*-

## /********************************************************/
## /********************************************************/
##             CHECK GAUSSIAN INITIAL PARAMETERS
## Function: Subroutine that extracts initial parameters  
## (program, funtional, basis set, if a file corresponds 
## to a frequency calculation, the termination, etc.)
## /********************************************************/
## /********************************************************/


## Python modules

import re
import sys
import os
from os import path

## Eyring modules

from level import *

def check_outputfile(filename):

    # Check initial paramaters from Gaussian output files
    
    all_program = [ ]
    all_termination = [ ]
    all_thermo = [ ]
    all_functional = [ ]
    all_basis = [ ]
    all_level = [ ]
    all_mpn = [ ]
    all_cc = [ ]
    
    for ifile in filename:
        if not path.isfile(ifile):
            print("The file", ifile, "does not exist!")
            sys.exit(1)
        else:
            program = "unknown"
            termination = "not normal"
            thermo = "false"
            functional = "unknown"
            basis = "unknown"
            mpn = "unknown"
            level = "unknown"
            cc = "unknown"
            with open(ifile,'r') as f:
                for line in f:
                    if len(line) != 0:    
                        if "Gaussian 09" in line or "Gaussian 16" in line:    
                            program = "Gaussian"
                        if "Normal termination of Gaussian 09" in line or "Normal termination of Gaussian 16":    
                            termination = "normal"
                        if "Thermochemistry" in line:    
                            thermo = "true"                            
                        line = line.strip()
                        if not line.startswith("# OF"):
                            if not line.startswith("#N"):
                                if line.startswith("#"):
                                    lin = re.split("[# / \s]+", line)                                    
                                    
                for istr in lin:
                    istr.replace('/', '')
                    for ifunc in functionals:
                        if ifunc.lower() == istr.lower():
                            functional = istr   # Functional
                        if ifunc.lower() == istr.lower().lstrip("r"):
                            functional = istr   # Functional
                        if ifunc.lower() == istr.lower().lstrip("u"):
                            functional = istr   # Functional
                        if ifunc.lower() == istr.lower().replace('/gen', ''):
                            functional = istr   # Functional
                            
                    for impn in MPn:
                        if impn.lower() == istr.lower():
                            mpn = istr   # MPn method
                        if impn.lower() == istr.lower().lstrip("r"):
                            mpn = istr   # MPn method
                        if impn.lower() == istr.lower().lstrip("u"):
                            mpn = istr   # MPn method
                    for icc in CC:
                        if icc.lower() == istr.lower():
                            cc = istr   # CC method
                    for ibasis in basis_sets:
                        if ibasis.lower() == istr.lower():
                            basis = istr    # Basis set
                
                if functional != "unknown":
                    if functional.lower().startswith("u"):
                        level = functional.lower().replace("u", "") + "/" + basis.lower()
                    elif functional.lower().startswith("r"):
                        level = functional.lower().replace("r", "") + "/" + basis.lower()  
                    else:
                        level = functional.lower() + "/" + basis.lower()
                elif mpn != "unknown":
                    if mpn.lower().startswith("u"):
                        level = mpn.lower().replace("u", "") + "/" + basis.lower()
                    elif mpn.lower().startswith("r"):
                        level = mpn.lower().replace("r", "") + "/" + basis.lower()  
                    else:
                        level = mpn.lower() + "/" + basis.lower()
                elif cc != "unknown":
                    level = cc + "/" + basis.lower()
                else:
                    level = "unknown" 
                    
        all_program.append(program)
        all_termination.append(termination)
        all_thermo.append(thermo)
        all_functional.append(functional)
        all_mpn.append(mpn)
        all_cc.append(cc)
        all_basis.append(basis)
        all_level.append(level)
        
    
        
                            
    outputfiledata = { }
    outputfiledata = {'program' : all_program,
                      'termination' : all_termination,
                      'thermo' : all_thermo,
                      'functional' : all_functional,
                      'mpn' : all_mpn,
                      'cc' : all_cc,
                      'basis' : all_basis,
                      'e_method': all_level}
    return outputfiledata

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


