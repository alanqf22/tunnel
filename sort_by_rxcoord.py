# -*- coding: utf-8 -*-

## /********************************************************/
## /********************************************************/
##                      SORT BY RXCOORD
## Function: Subroutine that sorts information from
## Intrinsic Reaction Coordinate (IRC) computations
## (performed with Gaussian 09 and 16) according to 
## the reaction coordinate in ascending order. The
## IRC information can be from the two-file version
## (IRC Forward and IRC Reverse) or the one-file 
## version (IRC complete).
## /********************************************************/
## /********************************************************/

## Eyringpy modules
from get_xyz_gaussian import *
from gaussian_parser import *

##----------------------------------------------------------
def sort_by_rxcoord(filename):

    irc = get_irc_gaussian(filename)
    rxcoord = irc['rxcoord']
    energy = irc['energy']
    charge = irc['charge'][0]
    multiplicity = irc['multiplicity'][0]
    reverse = irc['rev']
    forward = irc['for']
    forward_reverse = irc['both']
    xyz = get_xyz_gaussian(filename)
    atomic_coord = xyz['atomic']
    atom_coord = xyz['atom']
    x_coord = xyz['x']
    y_coord = xyz['y']
    z_coord = xyz['z']
    
    #Removing the last element from xyz
    atomic = [ ]
    atom = [ ]
    x = [ ]
    y = [ ]
    z = [ ]
    for i in range(len(x_coord)):
        atomic.append(atomic_coord[i][:-1])
        atom.append(atom_coord[i][:-1])
        x.append(x_coord[i][:-1])
        y.append(y_coord[i][:-1])
        z.append(z_coord[i][:-1])
        
    #Putting together the lists
    rxcoord_merge = [ ]
    energy_merge = [ ]
    atomic_merge = [ ]
    atom_merge = [ ]
    x_merge = [ ]
    y_merge = [ ]
    z_merge = [ ]
    for j in range(len(rxcoord)):
        rxcoord_merge = rxcoord_merge + rxcoord[j]
        energy_merge = energy_merge + energy[j]
        atomic_merge = atomic_merge + atomic[j]
        atom_merge = atom_merge + atom[j]
        x_merge = x_merge + x[j]
        y_merge = y_merge + y[j]
        z_merge = z_merge + z[j]

    #Sorting out by rxcoord (from negatives to positives)
    data_zipped = zip(rxcoord_merge, energy_merge, atomic_merge, atom_merge, x_merge, y_merge, z_merge)
    data_sorted = sorted(data_zipped, key=lambda x: float(x[0]))

    all_rxcoord = [ ]
    all_energy = [ ]
    all_atomic = [ ]
    all_atom = [ ]
    all_x = [ ]
    all_y = [ ]
    all_z = [ ]
    for k in range(len(data_sorted)):
        all_rxcoord.append(data_sorted[k][0])
        all_energy.append(data_sorted[k][1])
        all_atomic.append(data_sorted[k][2])
        all_atom.append(data_sorted[k][3])
        all_x.append(data_sorted[k][4])
        all_y.append(data_sorted[k][5])
        all_z.append(data_sorted[k][6])

    #Removing repeated values
    all_rxcoord_final = [ ]
    all_energy_final = [ ]
    all_atomic_final = [ ]
    all_atom_final = [ ]
    all_x_final = [ ]
    all_y_final = [ ]
    all_z_final = [ ]
    for m in range(len(all_rxcoord)):
        if all_rxcoord[m] not in all_rxcoord_final:
            all_rxcoord_final.append(all_rxcoord[m])
            all_energy_final.append(all_energy[m])
            all_atomic_final.append(all_atomic[m])
            all_atom_final.append(all_atom[m])
            all_x_final.append(all_x[m])
            all_y_final.append(all_y[m])
            all_z_final.append(all_z[m])
            
    #Reactant energy as the zero energy value
    min_energy = all_energy_final[0]
    all_energy_final = [(i - min_energy) for i in all_energy_final]

    irc_sorted = { }
    irc_sorted['rxcoord'] = all_rxcoord_final
    irc_sorted['energy'] = all_energy_final
    irc_sorted['atomic'] = all_atomic_final
    irc_sorted['atom'] = all_atom_final
    irc_sorted['x'] = all_x_final
    irc_sorted['y'] = all_y_final
    irc_sorted['z'] = all_z_final
    irc_sorted['charge'] = charge
    irc_sorted['multiplicity'] = multiplicity
    irc_sorted['for'] = forward
    irc_sorted['rev'] = reverse
    irc_sorted['both'] = forward_reverse
                             
    return irc_sorted



