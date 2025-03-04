# -*- coding: utf-8 -*-

## /********************************************************/
## /********************************************************/
##                           Spline
## Function: Subroutines that use the spline interpolation
## to interpolate, derive and integrate energy profiles.
## /********************************************************/
## /********************************************************/

## Python modules
import numpy as np
from math import sqrt
from scipy.interpolate import UnivariateSpline

## Eyringpy modules
from constants import *

def spline(x, y, xnew = None, interp_order = 3, smooth_factor = 0.0, npoints = 0):
    if xnew is None:
        xnew = np.linspace(min(x), max(x), npoints)
    spl = UnivariateSpline(x, y, k = interp_order, s = smooth_factor)
    ynew = list(spl(xnew))
    return ynew
    
def spline_derivative(x, y, xnew = None, npoints = 0, interp_order = 3, smooth_factor = 0.0, derivative_order = 1):
    if xnew is None:
        xnew = np.linspace(min(x), max(x), npoints)
    spl = UnivariateSpline(x, y, k = interp_order, s = smooth_factor)
    der = spl.derivative(derivative_order)
    ynew = list(der(xnew))
    return ynew

def spline_integral(x, y, xmin, xmax, interp_order = 3, smooth_factor = 0.0):
    spl = UnivariateSpline(x, y, k = interp_order, s = smooth_factor)
    integration = spl.integral(xmin, xmax)
    return integration

def mse(ytrue, ypred):
    difference_array = np.subtract(ytrue, ypred)
    squared_array = np.square(difference_array)
    mse = squared_array.mean()
    return(mse)
