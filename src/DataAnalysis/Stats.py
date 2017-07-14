"""!
@file DataAnalysis/Stats.py
@package DataAnalysis

@defgroup Stats Stats

@brief Statistics tools to perform error analysis on data.

This module supplements existing python capabilities and addresses gaps in
existing libraries.

@author James Bevins

@date 19Feb17
"""

import numpy as np

from sympy import symbols, diff
from math import sqrt

#------------------------------------------------------------------------------#
def red_chisq(yData, yMod, standDev=[], freeParams=2):
    """!
    @ingroup Stats

    Returns the reduced chi-square error statistic for an arbitrary model,
    \f$\frac{\chi^2}{\nu}\f$, where \f$\nu\f$ is the number of degrees of
    freedom. If individual standard deviations (array sd) are supplied,
    then the chi-square error statistic is computed as the sum of squared
    errors divided by the standard deviations.
    See http://physics.ucsc.edu/~drip/133/ch4.pdf for reference.

    @param yData: <em> Numpy array of integers or floats </em> \n
        Experimental data. \n
    @param yMod: <em> Numpy array of integers or floats </em> \n
        Model data. \n
    @param standDev: <em> Numpy array of integers or floats </em> \n
        Experimental data 1\f$\sigma\f$ standard deviation. \n
    @param freeParams: \e integer \n
        The number of free parameters in the model. \n

    @return \e float: The \f$\frac{\chi^2}{\nu}\f$ statistic \n
    """

    # Enable the use of lists as input
    if type(yData) == list:
        yData = np.asarray(yData)
    if type(yMod) == list:
        yMod = np.asarray(yMod)
    if type(standDev) == list:
        standDev = np.asarray(standDev)
        
    # Chi-square statistic
    if len(standDev) != len(yData):
        chisq = np.sum((yData-yMod)**2)
    else:
        chisq = np.sum(((yData-yMod)/standDev)**2)

    # Number of degrees of freedom
    nu = yData.size-freeParams

    return chisq/nu

#------------------------------------------------------------------------------#
def curve_fit_error4(func, args, indyErr, covar):
    """!
    @ingroup Stats

    Propogates error from a 4 parameter curve fit to a point estimate on the
    curve according to:
    http://ipnpr.jpl.nasa.gov/progress_report/42-122/122E.pdf
    http://www.itl.nist.gov/div898/handbook/mpc/section3/mpc3671.htm

    This requires that the function is partially differentiable.

    This was developed for the 4 factor or exponential fit to germanium detecors
    from:

    http://www.ezag.com/fileadmin/ezag/user-uploads/isotopes/pdf/
    Behavior_of_Several_Germanium_Detector_Full_Energy_Peak.pdf

    but would be applicable to any function with four free parameters.

    @param func: \e function \n
        The curve fit function used.  This must be partially differentiable wrt
        all variables.
    @param args: <em> list or array of floats or intergers </em> \n
        Arguments to the curve fit function.  There must be one independent
        variable and four free parameter variables
    @param indyErr: <em> float or interger </em> \n
        The error on the independent variable. \n
    @param covar: <em> array of floats or intergers </em> \n
        Covariance of the free parameters used for the fit. Must be 4 x 4
        array \n

    @return \e float: The standard deviation of the curve fit at the point
                      specified\n
    """
    assert len(args) == 5, "The number of arguments is not the 5 needed for \
                            curve_fit_error4."
    assert len(covar) == 4 and len(covar[0]) == 4, "The size of the covariance \
                        array is less than the number of free parameters (4)"
    assert hasattr(func, '__call__'), 'Invalid function handle'

    # Assign symbols to each variable
    x, a, b, c, d = symbols('x a b c d')

    # Take partial derivatives wrt each variable to get variance
    err = np.sqrt(np.diag(covar))
    variance = 0
    variance += diff(func(x, a, b, c, d), x).subs({x:args[0], a:args[1], \
                           b:args[2], c:args[3], d:args[4]})**2*indyErr**2
    variance += diff(func(x, a, b, c, d), a).subs({x:args[0], a:args[1], \
                           b:args[2], c:args[3], d:args[4]})**2*err[0]**2
    variance += diff(func(x, a, b, c, d), b).subs({x:args[0], a:args[1], \
                           b:args[2], c:args[3], d:args[4]})**2*err[1]**2
    variance += diff(func(x, a, b, c, d), c).subs({x:args[0], a:args[1], \
                           b:args[2], c:args[3], d:args[4]})**2*err[2]**2
    variance += diff(func(x, a, b, c, d), d).subs({x:args[0], a:args[1], \
                           b:args[2], c:args[3], d:args[4]})**2*err[3]**2

    # Add covariance terms
    variance += + 2*diff(func(x, a, b, c, d), a).subs({x:args[0], a:args[1], \
                           b:args[2], c:args[3], d:args[4]}) \
                   *diff(func(x, a, b, c, d), b).subs({x:args[0], a:args[1], \
                           b:args[2], c:args[3], d:args[4]})*covar[0, 1]
    variance += + 2*diff(func(x, a, b, c, d), a).subs({x:args[0], a:args[1], \
                           b:args[2], c:args[3], d:args[4]}) \
                   *diff(func(x, a, b, c, d), c).subs({x:args[0], a:args[1], \
                           b:args[2], c:args[3], d:args[4]})*covar[0, 2]
    variance += + 2*diff(func(x, a, b, c, d), a).subs({x:args[0], a:args[1], \
                           b:args[2], c:args[3], d:args[4]}) \
                   *diff(func(x, a, b, c, d), d).subs({x:args[0], a:args[1], \
                           b:args[2], c:args[3], d:args[4]})*covar[0, 3]
    variance += + 2*diff(func(x, a, b, c, d), b).subs({x:args[0], a:args[1], \
                           b:args[2], c:args[3], d:args[4]}) \
                   *diff(func(x, a, b, c, d), c).subs({x:args[0], a:args[1], \
                           b:args[2], c:args[3], d:args[4]})*covar[1, 2]
    variance += + 2*diff(func(x, a, b, c, d), b).subs({x:args[0], a:args[1], \
                           b:args[2], c:args[3], d:args[4]}) \
                   *diff(func(x, a, b, c, d), d).subs({x:args[0], a:args[1], \
                           b:args[2], c:args[3], d:args[4]})*covar[1, 3]
    variance += + 2*diff(func(x, a, b, c, d), c).subs({x:args[0], a:args[1], \
                           b:args[2], c:args[3], d:args[4]}) \
                   *diff(func(x, a, b, c, d), d).subs({x:args[0], a:args[1], \
                           b:args[2], c:args[3], d:args[4]})*covar[2, 3]

    return sqrt(variance)
