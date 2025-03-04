# -*- coding: utf-8 -*-

#Python libraries
import numpy as np

#Eyring modules
from funcs import *

## /********************************************************/
## /********************************************************/
##                          HESSIAN 
## Function: treats the data from the force constant matrix
## to convert it into the Hessian matrix, diagonalize it and
## obtain its eigen values, eigen vectors and frequencies.
## /********************************************************/
## /********************************************************/

##### -------------------------------------------------------------------- #####

def get_projectionmatrix(xcc,masses,v0=None):
    PROJECT_ROT = True
    PROJECT_TRA = True
    EPS_NORM = 1e-7

    '''
    Generates matrix to project translation and rotation coordinates (mass scaled/weighted)
    Other coordinate can be projected by introducing it using v0 (in mass-scaled)
    '''
    vecs = []
    nat  = len(masses)
    # PROJECT TRA IN HESS FOR FREQS
    if PROJECT_TRA:
        # translation
        sqrtmasses = [np.sqrt(mass) for mass in masses]
        b1 = [term if ii==0 else 0.0 for term in sqrtmasses for ii in range(3)]
        b2 = [term if ii==1 else 0.0 for term in sqrtmasses for ii in range(3)]
        b3 = [term if ii==2 else 0.0 for term in sqrtmasses for ii in range(3)]
        norm1 = np.linalg.norm(b1)
        norm2 = np.linalg.norm(b2)
        norm3 = np.linalg.norm(b3)
        b1 /= norm1
        b2 /= norm2
        b3 /= norm3
        vecs += [b1,b2,b3]
    # PROJECT ROT IN HESS FOR FREQS
    if PROJECT_ROT:
       # rotation
        b4 = np.zeros(len(xcc))
        b5 = np.zeros(len(xcc))
        b6 = np.zeros(len(xcc))
        for i in range(nat):
            b4[3*i + 1] =   np.sqrt(masses[i]) * z(xcc,i)
            b4[3*i + 2] = - np.sqrt(masses[i]) * y(xcc,i)
            b5[3*i + 0] = - np.sqrt(masses[i]) * z(xcc,i)
            b5[3*i + 2] =   np.sqrt(masses[i]) * x(xcc,i)
            b6[3*i + 0] =   np.sqrt(masses[i]) * y(xcc,i)
            b6[3*i + 1] = - np.sqrt(masses[i]) * x(xcc,i)
        norm4 = np.linalg.norm(b4)
        norm5 = np.linalg.norm(b5)
        norm6 = np.linalg.norm(b6)
        if norm4 > EPS_NORM: b4 /= norm4; vecs.append(b4)
        if norm5 > EPS_NORM: b5 /= norm5; vecs.append(b5)
        if norm6 > EPS_NORM: b6 /= norm6; vecs.append(b6)
    # Gram Schmidt
    if len(vecs) != 0:
        X = np.matrix(vecs).transpose()
        X_gs, R = np.linalg.qr(X)
        projmatrix = X_gs * X_gs.H
    else:
        projmatrix = np.zeros( (3*nat,3*nat) )
    
    # PROJECT GRADIENT
    #v0 = None
    if v0 is not None:
        normv0 = np.linalg.norm(v0)
        if normv0 > 1e-7:
            v0 = np.matrix( v0 ) / normv0
            projmatrix += v0.transpose() * v0
    
    return projmatrix

##### -------------------------------------------------------------------- #####

def project_hessian(Fms,natoms,proj_matrix):
    I = np.identity(3*natoms)
    Fms_proj = (I - proj_matrix) * Fms * (I - proj_matrix)
    return Fms_proj


##### -------------------------------------------------------------------- #####

def diagonalize(matrix):
    #Function: Diagonalizes a square array to compute its eigenvalues, eigenvectors and frequencies
    eigen_vals, eigen_vecs = np.linalg.eigh(matrix)
    #Convert eigenvalues to angular frequencies:
    freqs = [eval2angfreq(ieval) for ieval in eigen_vals]
    eigen_vecs  = eigen_vecs.transpose().tolist()
    return eigen_vals, eigen_vecs, freqs

##### -------------------------------------------------------------------- #####

def remove_extra_freqs(evals, evecs, freqs, v0, nvdof):
    if v0 is None: nvdof += 1
    #Function: Remove extra frequencies
    idxs = sorted([(abs(fq),fq,idx) for idx,fq in enumerate(freqs)])
    idxs.reverse()
    idxs = idxs[:nvdof]
    idxs = sorted([(fq,idx) for absfq,fq,idx in idxs])
    idxs = [idx for fq,idx in idxs]
    freqs = [freqs[idx] for idx in idxs]
    evals  = [ evals[idx] for idx in idxs]
    evecs  = [ evecs[idx] for idx in idxs]
    return evals, evecs, freqs

##### -------------------------------------------------------------------- #####
