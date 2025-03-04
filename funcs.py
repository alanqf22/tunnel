# -*- coding: utf-8 -*-

#Python libraries
import numpy as np

#Eyring modules
from constants import *

def x(data,num):
    return data[3*num+0]

##### -------------------------------------------------------------------- #####

def y(data,num):
    return data[3*num+1]

##### -------------------------------------------------------------------- #####

def z(data,num):
    return data[3*num+2]

##### -------------------------------------------------------------------- #####

def sign_of_number(x):
    #Function: Determines the sign of a number
    if x >= 0.0: return 1
    else: return -1

##### -------------------------------------------------------------------- #####

def eval2angfreq(eigenvalue):
    #Function: Converts eigenvalues to angular frequencies
    mu = 1/ATOMIC_MASS_UNIT
    return sign_of_number(eigenvalue) * (abs(eigenvalue)/mu)**0.5

##### -------------------------------------------------------------------- #####

def cart2mass_grad(cartesian_gradient,masses):
    #Function: Converts gradient from cartesian to mass-scaled
    mu = 1/ATOMIC_MASS_UNIT 
    #number of atoms:
    atom_num = len(masses)
    #corrected masses:
    mass_amu = [(masses[i]/mu)**0.5 for i in range(atom_num)]
    #mass-scaled gradient:
    msg = [ [x(cartesian_gradient,i)/mass_amu[i],y(cartesian_gradient,i)/mass_amu[i],z(cartesian_gradient,i)/mass_amu[i]] for i in range(atom_num)]
    return msg

##### -------------------------------------------------------------------- #####

def cart2mass_force(Fms,masses):
    #Function: Converts reaction force constant matrix from cartesian to mass-scaled
    mu = 1/ATOMIC_MASS_UNIT
    if Fms is None or len(Fms) == 0: return Fms
    #number of atoms:
    nat = len(masses)
    #mass-scaled reaction force constant matrix:
    Fcc = [ [ 0.0 for at1 in range(3*nat) ] for at2 in range(3*nat)]
    for i in range(3*nat):
        mi = masses[int(i/3)]
        for j in range(3*nat):
            mj = masses[int(j/3)]
            f = mu/((mi*mj)**0.5)
            Fcc[i][j] = Fms[i][j] * f
    return Fcc

##### -------------------------------------------------------------------- #####

def low2matrix (n,elements):
    #Function: Converts a lower triangular matrix to a square matrix
    #Square matrix size: nxn
    #Elements of the lower triangular matrix: elements

    # Create a lower triangular matrix initialized with zeros
    lower_triangular_matrix = np.zeros((n, n), dtype=float)

    # Fill the lower triangular matrix with the given elements
    count = 0
    for i in range(n):
        for j in range(i + 1):
            lower_triangular_matrix[i][j] = float(elements[count])
            count += 1

    # Create a symmetric matrix initially filled with zeros
    matrix = np.zeros((n, n))

    # Copy the values from the lower triangular array to the upper triangular part of new matrix (symmetric matrix)
    for i in range(n):
        for j in range(i, n):
            matrix[i][j] = lower_triangular_matrix[j][i]
            matrix[j][i] = lower_triangular_matrix[j][i]

    return matrix

##### -------------------------------------------------------------------- #####

def zero_point_energy(frequency):    # cal mol-1
    
    all_zpe = [ ]
    for ifreql in range(len(frequency)):
        zpe = 0.0
        for ifreq in frequency[ifreql]:
            zpe = zpe + ((0.5 * ifreq * PLANCK_CONSTANT * SPEED_OF_LIGHT * 100 * AVOGADRO_CONSTANT) / CALORIES_TO_JOULES)
        all_zpe.append(zpe)

    return all_zpe

##### -------------------------------------------------------------------- #####

def get_com(xcc,masses,indices=None):
    '''
    Returns the centre of mass of the selected atoms (indices)
    If no indices given, all atoms are considered
    '''
    if indices is None: indices = range(len(masses))
    tmass = sum([           masses[idx] for idx in indices])
    com_x = sum([x(xcc,idx)*masses[idx] for idx in indices])/tmass
    com_y = sum([y(xcc,idx)*masses[idx] for idx in indices])/tmass
    com_z = sum([z(xcc,idx)*masses[idx] for idx in indices])/tmass
    return [com_x,com_y,com_z]

##### -------------------------------------------------------------------- #####

def set_origin(xcc,x0):
    '''
    returns xcc with x0 as origin
    '''
    nat = howmanyatoms(xcc)
    return [xi-xj for xi,xj in zip(xcc,nat*x0)]

##### -------------------------------------------------------------------- #####

def howmanyatoms(xcc): return len(xcc)//3

##### -------------------------------------------------------------------- #####

def shift2com(xcc,masses):
    '''
    Function to shift to center of mass
    '''
    com = get_com(xcc,masses)
    xcc = set_origin(xcc,com)
    return xcc

##### -------------------------------------------------------------------- #####

#all_xcc = [-1.99108435e-16, 1.59473799, 0.0, -1.29360257, -1.26513574, 0.0, 2.35592723, 1.29132356, 0.0, 1.8502025, -1.25157884, 0.0]
#all_masses = [1.40030740E+01, 3.19720718E+01, 1.59949146E+01, 1.00782504E+00]

#a = shift2com(all_xcc,all_masses)
#print(a)
