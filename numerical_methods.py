#Python libraries
import numpy as np
from scipy.integrate import cumtrapz

def finite_difference_mixed(x, y):
    """
    Calculates the derivatives using centered finite differences for the intermediate points,
    forward finite differences for the first point,
    and backward finite differences for the last point.
    
    Parameters:
    x : array_like
        List of independent variable values.
    y : array_like
        List of dependent variable values.

    Returns:
    dydx : array_like
           Approximate derivatives using finite differences.
    """
    x = np.array(x)
    y = np.array(y)
    
    dydx = np.zeros_like(y)
    
    # Forward finite difference for the first point
    dydx[0] = (y[1] - y[0]) / (x[1] - x[0])
    
    # Centered finite differences for the intermediate points
    dydx[1:-1] = (y[2:] - y[:-2]) / (x[2:] - x[:-2])
    
    # Backward finite difference for the last point
    dydx[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])
    
    return dydx

def gauquad_pointsweights(ntp,x0=-1.0,xn=1.0):
    points, weights = np.polynomial.legendre.leggauss(ntp)
    suma    = (xn+x0)/2.0
    resta   = (xn-x0)/2.0
    points  = [resta*xi+suma for xi in points]
    weights = [wi*resta for wi in weights]
    return points, weights

def cumulative_trapz(x,y,xmin="none",xmax="none"):
    """
    Computes the cumulative integrated value of y using the compound trapezoidal rule
    of the scipy.integrate.cumtrapz method. The integration range is bounded by the given
    minimum and maximum values of x.

    Parameters:
    x    : array_like
           List of independent variable values.
    y    : array_like
           List of dependent variable values.

    xmin : scalar
           Lower integration value on the x-axis.

    xmax : scalar
           Upper integration value on the x-axis.

    Returns:
    dydx : array_like
           The result of cumulative integration of y along axis. 
    """
    all_x = [ ]
    all_y = [ ]
    if xmin == "none":
        all_x = x
        all_y = y
    else:
        for i in range(len(x)):
            if x[i] >= xmin and x[i] <= xmax:
                all_x.append(x[i])
                all_y.append(y[i])

    cumtrapz_val = cumtrapz(all_y, all_x)

    if len(cumtrapz_val) > 0:
        last_cum_trapezoidal = cumtrapz_val[-1]
    else:
        last_cum_trapezoidal = 0

    return last_cum_trapezoidal

