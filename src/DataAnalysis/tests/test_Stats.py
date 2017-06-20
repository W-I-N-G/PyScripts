"""!
This test suite evaluates all of the corner and edge cases for the functions
and classes in Stats.

@author James Bevins

@date 19Feb17
"""

import numpy as np

from Stats import red_chisq, curve_fit_error4

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in

#------------------------------------------------------------------------------#
def test_red_chisq():
    """!
    1) Test known cases
    """

    #1
    mass = np.array([91.161, 91.174, 91.186, 91.188])
    std = np.array([0.013, 0.011, 0.013, 0.013])
    avgMass = np.average(mass)
    assert_almost_equal(red_chisq(mass, avgMass, standDev=std, freeParams=1), \
                        0.928876, places=6)
    assert_almost_equal(red_chisq(mass, 91.17725, standDev=std, freeParams=1), \
                        0.928876, places=6)

    nMeas = np.array([6, 6, 7, 13, 8, 11, 2, 7, 6])
    nExp = np.array([6.5, 5.5, 7.6, 9.3, 9.9, 9.1, 7.3, 5.0, 5.8])
    std = np.sqrt(nExp)
    assert_almost_equal(red_chisq(nMeas, nExp, standDev=std, freeParams=3), \
                        1.169920, places=6)

#------------------------------------------------------------------------------#
def test_curve_fit_error4():
    """!
    1) Test known case
    2) Test inputs
    3) Test exceptions
    """

    #1
    def test(y, a, b, c, d):
        return -b + (b**2 - 4*c*(a-y))**(1/2.)/(2*c) + d
    var = [1, -1.83971535e-05, 1.00102485e-01, 7.03187162e-06, 0.00e+00]
    cov = np.array([[6.0046261e-10, -1.0757028e-10, 4.0166052e-12, 0.0000e+00],
                  [-1.07570279e-10, 2.34010372e-11, -9.49970127e-13, 0.0e+00],
                  [4.01660518e-12, -9.49970127e-13, 4.05090007e-14, 0.000e+00],
                  [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.0000e+00]])
    assert_almost_equal(curve_fit_error4(test, var, 0.0000376353, cov), \
                        204.062627963, places=6)

    #2
    var = np.array([1, -1.83971535e-05, 1.00102485e-01, 7.03187162e-06,
                    0.00e+00])
    assert_almost_equal(curve_fit_error4(test, var, 0.0000376353, cov), \
                        204.062627963, places=6)
    assert_almost_equal(curve_fit_error4(test, var, 1, cov), \
                        204.3063189, places=6)

    #3
    cov = [[6.0046261e-10, -1.0757028e-10, 4.0166052e-12, 0.0000e+00],
                  [-1.07570279e-10, 2.34010372e-11, -9.49970127e-13, 0.0e+00],
                  [4.01660518e-12, -9.49970127e-13, 4.05090007e-14, 0.000e+00],
                  [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.0000e+00]]
    assert_raises(TypeError, curve_fit_error4, test, var, 0.0000376353, cov)
    cov = np.array([[6.0046261e-10, -1.0757028e-10, 4.0166052e-12, 0.0000e+00],
                  [-1.07570279e-10, 2.34010372e-11, -9.49970127e-13, 0.0e+00],
                  [4.01660518e-12, -9.49970127e-13, 4.05090007e-14, 0.000e+00],
                  [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00e+00]])
    assert_raises(AssertionError, curve_fit_error4, var, var, 0.00376353, cov)
    cov = np.array([[6.0046261e-10, -1.0757028e-10, 4.0166052e-12, 0.0000e+00],
                  [-1.07570279e-10, 2.34010372e-11, -9.49970127e-13, 0.0e+00],
                  [4.01660518e-12, -9.49970127e-13, 4.05090007e-14, 0.00e+00]])
    assert_raises(AssertionError, curve_fit_error4, test, var, 0.00376353, cov)
    var = np.array([1, -1.83971535e-05, 1.00102485e-01, 7.03187162e-06])
    assert_raises(AssertionError, curve_fit_error4, test, var, 0.00376353, cov)
