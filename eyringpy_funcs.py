# -*- coding: utf-8 -*-

## /********************************************************/
## /********************************************************/
##                       UTILS
## Function: Subroutine that performs general tasks. 
## /********************************************************/
## /********************************************************/

## Python modules
import glob
import os
from os import path
from os.path import basename

##---------------------------------------------------------------------##
## Getting Gaussian output files
def get_out_file():
    
    out = [ ]
    for i in glob.glob("*.out"):
        out.append(i)

    return out
##---------------------------------------------------------------------##
## Getting xyz files
def get_xyz_file():
    
    xyz = [ ]
    for i in glob.glob("*.xyz"):
        xyz.append(i)
        
    return xyz
##---------------------------------------------------------------------##
## Calculating energy difference (DE or D(E+ZPE))
def energy_difference(energy):

    de = [ ]
    for i in range(len(energy)):
        ide = energy[i] - energy[0]
        de.append(ide)

    return de
##---------------------------------------------------------------------##
def get_output(path):
    
    out = [ ]
    for i in glob.glob(str(path)):
        out.append(i)
        
    return out
##---------------------------------------------------------------------##
def sort_files_by_int(file_list):

    files_sorted = sorted(file_list, key=lambda x:int(basename(os.path.splitext(x)[0])))

    return files_sorted
##---------------------------------------------------------------------##   
