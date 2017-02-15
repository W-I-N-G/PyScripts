# This test suite evaluates all of the corner and edge cases for the functions and classes in Counting.    
#
# @author James Bevins
#
# @date 13Feb17

from BasicNuclearCalcs import *

import nose
from nose.plugins.skip import SkipTest
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in
    
#-------------------------------------------------------------------------------------------------------------#
def test_production_decay():    
    """!
    1) Test given irradiation time << half life
    2) Test given low production and short half life
    3) Test post-irradiation decay
    4) Test initial population with decay
    5) Test exceptions
    
    General: Test that both int and float values work for all inputs
    """
    
    #1
    assert_almost_equal(production_decay(1E10, 0, 100, 1E-3, 1E6, 1), 1E5, places=0)
    
    #2
    assert_almost_equal(production_decay(1E-3, 0, 100, 1E-3, 1E3, 1), 0.001442695, places=6)
    
    #3
    assert_almost_equal(production_decay(1E10, 0, 100, 1E-3, 1E6, 1, 1E10), 0.5E5, places=0)
    
    #4
    assert_almost_equal(production_decay(100, 1000, 100, 1E-3, 1E3, 1), 572.1347520, places=6)
    
    #5
    assert_raises(TypeError,production_decay,'one', 1000, 100, 1E-3, 1E3, 1)
    assert_raises(AssertionError,production_decay,-100, 1000, 100, 1E-3, 1E3, 1)
    assert_raises(AssertionError,production_decay,100, -1000, 100, 1E-3, 1E3, 1)
    assert_raises(AssertionError,production_decay,100, 1000, -100, 1E-3, 1E3, 1)
    assert_raises(AssertionError,production_decay,100, 1000, 100, -1E-3, 1E3, 1)
    assert_raises(AssertionError,production_decay,100, 1000, 100, 1E-3, -1E3, 1)
    assert_raises(AssertionError,production_decay,100, 1000, 100, 1E-3, 1E3, -1)
    
#-------------------------------------------------------------------------------------------------------------#
def test_decay():    
    """!
    1) Test known cases
    2) Test wrong units
    3) Test Exceptions
    
    General: Test that both int and float values work for all inputs
    """
    
    #1
    assert_equal(decay(100.0, 1000, 100, units='atoms'), 500)
    assert_equal(decay(100, 1000.0, 100), 18500000.0)
    assert_almost_equal(decay(1E10, 1000, 1.0, units='Bq'), 1000, places=5)
    
    #2
    assert_almost_equal(decay(1E10, 1000, 1.0, units='Frogs'), 1000, places=5)
    assert_almost_equal(decay(1E10, 1000, 1.0, units=5), 1000, places=5)
    
    #3
    assert_raises(TypeError,decay,1E10,1000,'one')
    assert_raises(AssertionError,decay,-1E10,1000,1)
    assert_raises(AssertionError,decay,1E10,-1000,1)
    assert_raises(AssertionError,decay,1E10,1000,-1)
    
#-------------------------------------------------------------------------------------------------------------#
def test_get_decay_const():    
    """!
    1) Test known cases
    2) Exceptions
    
    General: Test that both int and float values work for all inputs
    """
    
    #1
    assert_almost_equal(get_decay_const(100), 6.9314728E-3, places=5)
    assert_almost_equal(get_decay_const(50.0), 0.01386294, places=5)
    assert_almost_equal(get_decay_const(1E8), 6.9314718E-9, places=5)
    
    #2
    assert_raises(AssertionError,get_decay_const,-10)
    
#-------------------------------------------------------------------------------------------------------------#
def test_get_halflife():    
    """!
    1) Test known cases
    2) Exceptions
    
    General: Test that both int and float values work for all inputs
    """
    
    #1
    assert_almost_equal(get_halflife(100), 6.9314728E-3, places=5)
    assert_almost_equal(get_halflife(50.0), 0.01386294, places=5)
    assert_almost_equal(get_halflife(1E-8), 6.9314718E7, places=0)
    
    #2
    assert_raises(AssertionError,get_halflife,-10)
    
#-------------------------------------------------------------------------------------------------------------#
def test_activity():    
    """!
    1) Test known cases
    2) Exceptions
    
    General: Test that both int and float values work for all inputs
    """
    
    #1
    assert_almost_equal(activity(100,1000), 6.9314728, places=5)
    assert_almost_equal(activity(100,1000,100), 3.4657359, places=5)
    
    #2
    assert_raises(AssertionError,activity,-100,1000,100)
    assert_raises(AssertionError,activity,100,-1000,100)
    assert_raises(AssertionError,activity,100,1000,-100)
    
#-------------------------------------------------------------------------------------------------------------#
def test_solid_angle():    
    """!
    1) Test known cases
    2) Exceptions
    
    General: Test that both int and float values work for all inputs
    """
    
    #1
    assert_almost_equal(solid_angle(100.0,0), 6.283185307, places=5)
    assert_almost_equal(solid_angle(5,100.0), 7.839286E-3, places=5)
    assert_equal(solid_angle(0,10), 0)
    
    #2
    assert_raises(AssertionError,solid_angle,-100,10)
    assert_raises(AssertionError,solid_angle,100,-10)
    
#-------------------------------------------------------------------------------------------------------------#
def test_fractional_solid_angle():    
    """!
    1) Test known cases
    2) Exceptions
    
    General: Test that both int and float values work for all inputs
    """
    
    #1
    assert_almost_equal(fractional_solid_angle(100,0), 0.5, places=5)
    assert_almost_equal(fractional_solid_angle(5.0,100), 6.2343056E-4, places=5)
    assert_equal(fractional_solid_angle(0,10.0), 0)
    
    #2
    assert_raises(AssertionError,fractional_solid_angle,-100,10)
    assert_raises(AssertionError,fractional_solid_angle,100,-10)
    