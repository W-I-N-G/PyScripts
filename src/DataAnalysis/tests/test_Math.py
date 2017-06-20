"""!
This test suite evaluates all of the corner and edge cases for the functions
and classes in Math.

@author James Bevins

@date 21Feb17
"""

import os

import pandas as pd
import numpy as np

from datetime import datetime

from Math import gauss, smeared_step, skew_gauss, quadratic

import nose
from nose.plugins.skip import SkipTest
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in

#------------------------------------------------------------------------------#
def test_gauss():
    """!
    1) Test values computed by hand.
    2) Test functionality if inputs are negative.
    3) Test exceptions

    General: Test that both int and float values work for all inputs
    """

    #1
    assert_almost_equal(gauss(99, 1000, 100, 2), 882.49690, places=4)
    assert_almost_equal(gauss(99., 1000., 100., 2.), 882.49690, places=4)
    assert_equal(gauss(100., 1000., 100., 2.), 1000)
    assert_almost_equal(gauss(50., 1000., 100., 2.), 1.9919E-133, places=4)

    #2
    assert_equal(gauss(-50., 1000., 100., 2.), 0)
    assert_equal(gauss(90., 1000., -100., 2.), 0)
    assert_almost_equal(gauss(99, 1000, 100, -2), 882.49690, places=4)
    assert_almost_equal(gauss(99, -1000, 100, 2), -882.49690, places=4)

    #3
    assert_raises(TypeError, gauss, 99, 1000, "ten", 2)
    assert_raises(ValueError, gauss, 99, 1000, 100, "two")

#------------------------------------------------------------------------------#
def test_smeared_step():
    """!
    1) Test values computed by hand.
    2) Test functionality if inputs are negative.
    3) Test exceptions

    General: Test that both int and float values work for all inputs
    """

    #1
    assert_almost_equal(smeared_step(99, 100, 2, 1000), 387.4556, places=4)
    assert_almost_equal(smeared_step(99., 100., 2., 1000.), 387.4556, places=4)
    assert_equal(smeared_step(100, 100, 2, 1000), 250)
    assert_almost_equal(smeared_step(50, 100, 2, 1000), 1000, places=4)

    #2
    assert_equal(smeared_step(-50, 100, 2, 1000), 1000)
    assert_almost_equal(smeared_step(90, 100, 2, 1000), 986.65909, places=4)
    assert_almost_equal(smeared_step(99, 100, -2, 1000), 142.5370, places=4)
    assert_almost_equal(smeared_step(99, 100, 2, -1000), -387.4556, places=4)

    #3
    assert_raises(TypeError, smeared_step, 99, "ten", 1000, 2)
    assert_raises(TypeError, smeared_step, 99, 1000, 100, "two")

#------------------------------------------------------------------------------#
def test_skew_gauss():
    """!
    1) Test values computed by hand.
    2) Test functionality if inputs are negative.
    3) Test exceptions

    General: Test that both int and float values work for all inputs
    """

    #1
    assert_almost_equal(skew_gauss(99, 100, 2, 1000, 5), 12.3228, places=4)
    assert_almost_equal(skew_gauss(99., 100., 2., 1000., 5), 12.3228, places=4)
    assert_equal(skew_gauss(100, 100, 2, 1000, 5), 62.5)
    assert_almost_equal(skew_gauss(50, 100, 2, 1000, 5), 5.1664E-50, places=4)

    #2
    assert_almost_equal(skew_gauss(-50, 100, 2, 1000, 5), 1.38E-160, places=2)
    assert_almost_equal(skew_gauss(99, 100, -2, 1000, 5), 247.5091, places=4)
    assert_almost_equal(skew_gauss(99, 100, 2, -1000, 5), -12.3228, places=4)

    #3
    assert_raises(TypeError, skew_gauss, 99, "ten", 1000, 2, 5)
    assert_raises(TypeError, skew_gauss, 99, 1000, 100, "two", 5)

#------------------------------------------------------------------------------#
def test_quadratic():
    """!
    1) Test values computed by hand.
    2) Test functionality if inputs are negative.
    3) Test exceptions

    General: Test that both int and float values work for all inputs
    """

    #1
    assert_equal(quadratic(10, 5, 20, 100), 800)
    assert_equal(quadratic(10., 5., 20., 100.), 800)
    assert_equal(quadratic(10, 0, 0, 100), 100)
    assert_equal(quadratic(10, .5, 20, 100), 350)

    #2
    assert_equal(quadratic(-10, 5, 20, 100), 400)
    assert_equal(quadratic(10, -5, 20, 100), -200)
    assert_equal(quadratic(10, 5, -20, 100), 400)
    assert_equal(quadratic(10, 5, 20, -100), 600)

    #3
    assert_raises(TypeError, quadratic, "ten", 5, 20, 100)
    assert_raises(TypeError, quadratic, 10, "five", 20, -100)
