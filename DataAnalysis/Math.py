"""!
@file DataAnalysis/Math.py
@package DataAnalysis

@defgroup Math Math

@brief Math support functions and methods.

This module supplements existing python capabilities and addresses gaps in
existing libraries.

@author James Bevins

@date 21Feb17
"""

import numpy as np

#------------------------------------------------------------------------------#
def gauss(x, amplitude, centroid, width):
    """!
    @ingroup DataAnalysis
    Evaluates a gaussian at a given point provided the amplitude, centroid,
    and width.

    @param x: <em> integer or float </em>  \n
        Point at which to evaluate \n
    @param amplitude: <em> integer or float </em>  \n
        Amplitude of the peak \n
    @param centroid: <em> integer or float </em>  \n
        Location of the centroid in the same units as x \n
    @param width: <em> integer or float </em>  \n
        Width of the distribution in the same units as x \n

    @return \e float: evaluated value at point x \n
    """
    z = (x-centroid)/float(width)
    return amplitude*np.exp(-0.5*z**2)

#------------------------------------------------------------------------------#
def smeared_step(x, centroid, width, amplitude):
    """!
    @ingroup DataAnalysis
    Evaluates a step function convoluted with a gaussian at a given point
    provided the amplitude, centroid, and width.

    @param x: <em> integer or float </em>  \n
        Point at which to evaluate \n
    @param centroid: <em> integer or float </em>  \n
        Location of the centroid in the same units as x \n
    @param width: <em> integer or float </em>  \n
        Width of the distribution in the same units as x \n
    @param amplitude: <em> integer or float </em>  \n
        Amplitude of the peak \n

    @return \e float: evaluated value at point x \n
    """
    z = (x-centroid)/float(width)
    return amplitude/(1+np.exp(z))**2

#------------------------------------------------------------------------------#
def skew_gauss(x, centroid, width, amplitude, rng):
    """!
    @ingroup DataAnalysis
    Evaluates a skewed gaussian at a given point provided the amplitude,
    centroid, width, and range of the distribution.

    @param x: <em> integer or float </em>  \n
        Point at which to evaluate \n
    @param amplitude: <em> integer or float </em>  \n
        Amplitude of the peak \n
    @param centroid: <em> integer or float </em>  \n
        Location of the centroid in the same units as x \n
    @param width: <em> integer or float </em>  \n
        Width of the distribution in the same units as x \n
    @param rng: <em> integer or float </em>  \n
        Range of the distribution in the same units as x \n

    @return \e float: evaluated value at point x \n
    """
    z = (x-centroid)/float(width)
    return amplitude*(np.exp(rng*z))/(1+np.exp(z))**4

#------------------------------------------------------------------------------#
def quadratic(x, quad, linear, offset):
    """!
    @ingroup DataAnalysis
    Quadratic function.

    @param x: <em> integer or float </em>  \n
        Point at which to evaluate \n
    @param quad: <em> integer or float </em>  \n
        Quadratic term of the function \n
    @param linear: <em> integer or float </em>  \n
        Linear term of the function \n
    @param offset: <em> integer or float </em>  \n
        Offset term of the function \n

    @return \e float: evaluated value at point x \n
    """
    return quad*x**2 + linear*x + offset
