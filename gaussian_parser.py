#Python libraries
import re
import sys
from os import path
from math import isnan
#import numpy as np

#Eyring modules
from level import *
from atomic import *
from utils import *
from check_gauss_initparam import *

## /********************************************************/
## /********************************************************/
##                     GAUSSIAN_PARSER 
## Function: parsers Gaussian *.out, *.log and *.fchk files.
## /********************************************************/
## /********************************************************/

##### -------------------------------------------------------------------- #####
#####                        IRC GAUSSIAN FILE PARSER                      #####
##### -------------------------------------------------------------------- #####

def get_irc_gaussian(filename):
    outputfiledata = check_outputfile(filename)
    program = outputfiledata['program']
    termination = outputfiledata['termination']
    all_rxcoord = [ ]
    all_energy = [ ]
    all_charge = [ ]
    all_multiplicity = [ ]
    file_for = [ ]
    file_rev = [ ]
    file_both = [ ]

    for ifile in range(len(filename)):
        energy = [ ]
        rxcoord = [ ]
        if program[ifile] != "Gaussian":
            print (filename[ifile], "is not a Gaussian output file!")
            sys.exit(1)
        else:
            if termination[ifile] != "normal":
                print (filename[ifile], "calculation did not terminate normally!")
                sys.exit(1)
            else:
                readlines = 0
                with open(filename[ifile],'r') as f:
                    for line in f:
                        line = line.strip()
                        lin = line.split()
                        if len(lin) > 0:
                            if "Energy" in line and "RxCoord" in line:
                                readlines = 1
                            if "----------------------------" in line:
                                readlines = 0
                            if readlines == 1:
                                if lin[0].isdigit():
                                    li = line.replace("-", " -").split()
                                    energy.append(float(li[1]))
                                    rxcoord.append(float(li[2]))
                            if "Charge" in line and "Multiplicity" in line:
                                charge = int(lin[2])
                                multiplicity = int(lin[5])

        #IRC Forward
        if rxcoord[0] == 0.0 and rxcoord[-1] > 0.0:
            file_for.append(filename[ifile])
            all_rxcoord.append(rxcoord)
            all_energy.append(energy)

        #IRC Reverse
        if rxcoord[0] < 0.0 and rxcoord[-1] == 0.0:
            rxcoord_rev = [ ]
            energy_rev = [ ]
            for i,j in zip(rxcoord, energy):
                rxcoord_rev.insert(0, i)
                energy_rev.insert(0, j)
            file_rev.append(filename[ifile])
            all_rxcoord.append(rxcoord_rev)
            all_energy.append(energy_rev)

        #IRC (Forward and Reverse together)
        if rxcoord[0] < 0.0 and rxcoord[-1] > 0.0:
            rxcoord_for = [ ]
            rxcoord_rev = [ ]
            energy_for = [ ]
            energy_rev = [ ]
            for i, j in zip(rxcoord,energy):
                if i >= 0.0:  #IRC Forward  
                    rxcoord_for.append(i)
                    energy_for.append(j)
                if i < 0.0:   #IRC Reverse
                    rxcoord_rev.insert(0, i)
                    energy_rev.insert(0, j)
            file_both.append(filename[ifile])
            all_rxcoord.append(rxcoord_for + rxcoord_rev)
            all_energy.append(energy_for + energy_rev)

    all_charge.append(charge)
    all_multiplicity.append(multiplicity)

    irc = { }
    irc['rxcoord'] = all_rxcoord
    irc['energy'] = all_energy
    irc['charge'] = all_charge
    irc['multiplicity'] = all_multiplicity
    irc['for'] = file_for 
    irc['rev'] = file_rev
    irc['both'] = file_both

    return irc

##### -------------------------------------------------------------------- #####
#####                     FREQUENCY GAUSSIAN FILE PARSER                   #####
##### -------------------------------------------------------------------- #####

def get_frequency(filename):
    all_freq_p = [ ]
    all_freq_n = [ ]
    all_nnfreq = [ ]
    for i in range(len(filename)):
        freq_p = [ ]
        freq_n = [ ]
        with open(filename[i],'r') as f:
                for line in f:
                    lin = line.split()
                    if len(lin) > 0:
                        # Frequencies (cm-1)
                        if lin[0] == "Frequencies" and lin[1] == "--":
                            for i in range(2,len(lin)):
                                if float(lin[i]) > 0.0:
                                    freq_p.append(float(lin[i]))     # Positive frequencies  
                                else:
                                    freq_n.append(abs(float(lin[i])))    # Negative frequencies

        nnfreq = int(len(freq_n))    # Number of negative frequencies
        all_freq_p.append(freq_p)
        all_freq_n.append(freq_n)
        all_nnfreq.append(nnfreq)
        
    frequency = {'positive' : all_freq_p,
                 'negative' : all_freq_n,
                 'negative_frequency' : all_nnfreq}
    
    return frequency

##### -------------------------------------------------------------------- #####

def get_zpe(filename):
    all_zpe = [ ]
    for i in range(len(filename)):
        zpe = 0.0   
        with open(filename[i],'r') as f:
                for line in f:
                    lin = line.split()
                    if len(lin) > 0:
                        if "Zero-point correction" in line:
                            if isnan(float(lin[2])) is False:
                                zpe = float(lin[2])
        all_zpe.append(zpe)
    return all_zpe

##### -------------------------------------------------------------------- #####

def get_electronic_energy(filename):
    all_eelec = [ ]
    for i in range(len(filename)):
        readline = 0
        electronic_energy = 0.0
        with open(filename[i],'r') as f:
                for line in f:
                    lin = line.split()
                    if len(lin) > 0:
                        # Electronic energy (Hartree particle-1)
                        if "SCF Done" in line:    # DFT
                            linex = line.replace('\x00', '')
                            linx = linex.split()
                            
                            if len(linx) < 10:
                                electronic_energy = float(linx[4])
                        if "EUMP2" in line:    # MP2
                            electronic_energy = float(lin[5].replace('D', 'E'))
                        if "T5(CCSD)" in line:  
                            readline = 1
                        if "Discarding MO integrals." in line:
                            readline = 0
                        if readline == 1:
                            if "CCSD(T)" in line:    # CCSD(T)
                                electronic_energy = float(lin[1].replace('D', 'E'))
        all_eelec.append(electronic_energy)
    return all_eelec

##### -------------------------------------------------------------------- #####
#####                            FCHK FILE PARSER                          #####
##### -------------------------------------------------------------------- #####

def get_cartesian_gradient(filename):
    #Function to parser Cartesian Gradient from fchk Gaussian files
    all_CG = [ ]
    for i in range(len(filename)):
        readline = 0
        CG = [ ]
        with open(filename[i],'r') as f:
                for line in f:
                    lin = line.split()
                    if len(lin) > 0:
                        if "Cartesian Gradient" in line:
                            readline = 1
                        elif "Dipole Moment" in line:
                            readline = 0
                        elif "Nonadiabatic coupling" in line:
                            readline = 0
                        elif "Cartesian Force Constants" in line:
                            readline = 0
                        if readline == 1:
                            if is_number(lin[0])== True:
                                for i in range(len(lin)):
                                    CG.append(float(lin[i]))
        all_CG.append(CG)
    return all_CG

##### -------------------------------------------------------------------- #####

def get_cartesian_force_constants(filename):
    #Function to parser Cratesian Force Constants from fchk Gaussian files
    all_CFC = [ ]
    for i in range(len(filename)):
        readline = 0
        CFC = [ ]
        with open(filename[i],'r') as f:
                for line in f:
                    lin = line.split()
                    if len(lin) > 0:
                        if "Cartesian Force Constants" in line: readline = 1
                        elif "Dipole Moment" in line: readline = 0
                        elif "Nonadiabatic coupling" in line: readline = 0
                        if readline == 1:
                            if is_number(lin[0])== True:
                                for i in range(len(lin)):
                                    CFC.append(float(lin[i]))
        all_CFC.append(CFC)
    return all_CFC

##### -------------------------------------------------------------------- #####

def get_atomic_weights(filename):
    #Function to parser Real Atomic Weights from fchk Gaussian files
    all_AW = [ ]
    for i in range(len(filename)):
        readline = 0
        AW = [ ]
        with open(filename[i],'r') as f:
                for line in f:
                    lin = line.split()
                    if len(lin) > 0:
                        if "Real atomic weights" in line: readline = 1
                        if "Atom fragment info" in line: readline = 0
                        if readline == 1:
                            if is_number(lin[0])== True:
                                for i in range(len(lin)):
                                    AW.append(float(lin[i]))
        all_AW.append(AW)
    return all_AW

##### -------------------------------------------------------------------- #####

def get_cartesian_coords(filename):
    #Function to parser Current Cartesian Coordinates from fchk Gaussian files
    all_coords = [ ]
    for i in range(len(filename)):
        readline = 0
        coords = [ ]
        with open(filename[i],'r') as f:
                for line in f:
                    lin = line.split()
                    if len(lin) > 0:
                        if "Current cartesian coordinates" in line: readline = 1
                        if "Number of symbols" in line: readline = 0
                        if readline == 1:
                            if is_number(lin[0])== True:
                                for i in range(len(lin)):
                                    coords.append(float(lin[i]))
        all_coords.append(coords)
    return all_coords

def get_electronic_energy_fchk(filename):
    #Function to parser electronic energy from fchk Gaussian files
    all_eelec = [ ]
    for i in range(len(filename)):
        electronic_energy = 0.0
        with open(filename[i],'r') as f:
                for line in f:
                    lin = line.split()
                    if len(lin) > 0:
                        if "SCF Energy" in line:
                            electronic_energy = float(lin[3])
        all_eelec.append(electronic_energy)
    return all_eelec
