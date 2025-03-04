# -*- coding: utf-8 -*-

## /********************************************************/
## /********************************************************/
##                     GET XYZ GAUSSIAN
## Function: Subroutine that extracts the xyz coordinates
## from Gaussian (09 and 16) output files.
## /********************************************************/
## /********************************************************/

## Python modules
import sys

## Eyringpy modules
from atomic import *

##----------------------------------------------------------
def get_xyz_gaussian(filename):
    all_atomic = [ ]
    all_atom = [ ]
    all_x = [ ]
    all_y = [ ]
    all_z = [ ]
    for i in range(len(filename)):
        atomics = [ ]
        atoms = [ ]
        xs = [ ]
        ys = [ ]
        zs = [ ]
        with open(filename[i],'r') as f:
            for line in f:
                lin = line.split()
                if len(lin) > 0:
                    if "Input orientation:" in line or "Standard orientation:" in line:
                        atomic = [ ]
                        x = [ ]
                        y = [ ]
                        z = [ ]
                        for i in range(5):
                            line = next(f)
                        while not line.startswith("--------"):
                            lin = line.split()
                            if len(lin) == 6 and lin[0].isdigit() and lin[1].isdigit() and lin[2].isdigit():
                                atomic.append(int(lin[1]))
                                x.append(float(lin[3]))
                                y.append(float(lin[4]))
                                z.append(float(lin[5]))
                            else:
                                break
                            line = next(f)
                            
                        atom = get_chemical_symbol(atomic)
                        atomics.append(atomic)
                        atoms.append(atom)
                        xs.append(x)
                        ys.append(y)
                        zs.append(z)
        all_atomic.append(atomics)
        all_atom.append(atoms)
        all_x.append(xs)
        all_y.append(ys)
        all_z.append(zs)

        

    xyz = { }
    xyz['atomic'] = all_atomic
    xyz['atom'] = all_atom
    xyz['x'] = all_x 
    xyz['y'] = all_y
    xyz['z'] = all_z

    return xyz
