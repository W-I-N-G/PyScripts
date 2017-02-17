# This test suite evaluates all of the corner and edge cases for the functions and classes in Counting.    
#
# @author James Bevins
#
# @date 13Feb17

import os

import pandas as pd

from datetime import datetime

from Counting import *
from BasicNuclearCalcs import fractional_solid_angle

import nose
from nose.plugins.skip import SkipTest
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in
    
#-------------------------------------------------------------------------------------------------------------#
def test_volume_solid_angle():    
    """!
    1) Test values computed by this function against those reported in the literature.  
    2) Test functionality if radius of detector is zero.
    3) Test functionality if radius of the source is zero.
    4) Test that converges to point source when radius of detector >> radius of source.
    5) Test that converges to point source when distance between detector and source is large.    
    6) Test functionality if distance between the detector and source is zero.
    7) Test that user errors in input type are caught
    
    General: Test that both int and float values work for all inputs
    """
    
    #1
    assert_almost_equal(volume_solid_angle(1, 0.5, 1),0.0343, places=4)
    assert_almost_equal(volume_solid_angle(1, 4, 1),0.3761, places=4)
    assert_almost_equal(volume_solid_angle(0.3, 2.54, 20.0),0.0040, places=4)
    assert_almost_equal(volume_solid_angle(2.0, 2.54, 5),0.0501, places=4)
    
    #2
    assert_equal(volume_solid_angle(2.0, 0.0, 3),0)
    
    #3
    assert_almost_equal(volume_solid_angle(0.0, 2.0, 3.0),fractional_solid_angle(2.0, 3.0), places=5)
    
    #4
    assert_almost_equal(volume_solid_angle(0.1, 10, 3.0),fractional_solid_angle(10.0, 3.0), places=5)
    
    #5
    assert_almost_equal(volume_solid_angle(2.54, 2.54, 300.0),fractional_solid_angle(2.54, 300), places=5)
    
    #6   
    assert_raises(AssertionError,volume_solid_angle, 2.54, 2.54, 0)
    
    #7
    assert_raises(TypeError,volume_solid_angle,"two", 2.54, 1)
    assert_raises(TypeError,volume_solid_angle,2.0, "two", 1)
    assert_raises(TypeError,volume_solid_angle,2.0, 2, "one")

#-------------------------------------------------------------------------------------------------------------#
def test_germanium_eff():
    """!
    1) Test that ouput equals hand calculated values
    2) Test that each input can be left off the input param list
    3) Test 0 for each parameter
    4) Test negative for each parameter
    5) Test that user input error exceptions work 
    """
    
    #1
    assert_almost_equal(germanium_eff(100),0.1114059, places=4)
    assert_almost_equal(germanium_eff(1000),0.0244010, places=4)
    assert_almost_equal(germanium_eff(1500),0.0148815, places=4)
    assert_almost_equal(germanium_eff(2000),0.00872368, places=4)
    
    #2
    assert_almost_equal(germanium_eff(250,a=0.05),0.2381597, places=4)
    assert_almost_equal(germanium_eff(500,a=0.05,b=0.03),-0.1997285, places=4)
    assert_almost_equal(germanium_eff(750,a=0.05,b=0.03,c=0.60),0.1337239, places=4)
    assert_almost_equal(germanium_eff(1000,a=0.05,b=0.03,c=0.60,d=0.003),0.13997, places=4)
    
    #3
    assert_almost_equal(germanium_eff(250,a=0),-0.2618402, places=4)
    assert_almost_equal(germanium_eff(500,a=0.05,b=0),0.6099624, places=4)
    assert_almost_equal(germanium_eff(750,a=0.05,b=0.03,c=0),-0.362234, places=4)
    assert_almost_equal(germanium_eff(1000,a=0.05,b=0.03,c=0.60,d=0),0.13997, places=4) 
    
    #4
    assert_almost_equal(germanium_eff(250,a=-0.05),-0.7618402, places=4)
    assert_almost_equal(germanium_eff(500,a=0.05,b=-0.03),1.4196534, places=4)
    assert_almost_equal(germanium_eff(750,a=0.05,b=0.03,c=-0.60),-0.8581933, places=4)
    assert_almost_equal(germanium_eff(1000,a=0.05,b=0.03,c=0.60,d=-0.003),0.140030, places=4)
    
    #5
    assert_raises(TypeError,germanium_eff, "five")
    assert_raises(TypeError,germanium_eff, 5, "five")
    assert_raises(TypeError,germanium_eff, 5, 5, "five")
    assert_raises(TypeError,germanium_eff, 5, 5, 5, "five")
    assert_raises(TypeError,germanium_eff, 5, 5, 5, 5, "five")
    

#-------------------------------------------------------------------------------------------------------------#
def test_parse_spe():
    """!
    1) Test that a sample spe can be read in and output expected results.  
    2) Test exceptions
    """
    
    #1
    (ct,dt,a,b,c,df)=parse_spe("/home/pyne-user/Dropbox/UCB/Computational_Tools/Scripts/Python/GeneralNuclear/tests/testFiles/test_parse.Spe")
    assert_equal(ct,120)
    assert_equal(dt,datetime.strptime('2/8/2017 18:39:13','%m/%d/%Y %H:%M:%S'))
    assert_equal(a,3.468863E-008)
    assert_equal(b,3.838918E-001)
    assert_equal(c,1.874105E-001)
    assert_true(df.equals(pd.DataFrame(['       1', '      46', '     303', '     417', '     443', '     412', '     415', '     435', '     452', '     442', '     436', '     436', '     450', '     458', '     436'],columns=['counts'])))
    
    #2
    assert_equal(parse_spe("test2.Spe"),None)
    
#-------------------------------------------------------------------------------------------------------------#
def test_peak_counts():
    """!
    1) Test output given a known result from sample spectrum
    2) Test output at different width given known results from sample spectrum; test float input parameters
    3) Test use of lists instead of arrays
    4) Test if passed an areas without a peak
    """
    
    (ct,dt,a,b,c,df)=parse_spe("/home/pyne-user/Dropbox/UCB/Computational_Tools/Scripts/Python/GeneralNuclear/tests/testFiles/test_peak_counts.Spe")    
    
    #1
    assert_almost_equal(peak_counts(np.asarray(df.index), np.asarray(df['counts']).astype(float), 1723)[0], 65901, places=0)
    assert_almost_equal(peak_counts(np.asarray(df.index), np.asarray(df['counts']).astype(float), 1723)[1], 258, places=0)
    
    #2
    assert_almost_equal(peak_counts(np.asarray(df.index), np.asarray(df['counts']).astype(float), 1723.5, width=100.5)[0], 66008, places=0)
    assert_almost_equal(peak_counts(np.asarray(df.index), np.asarray(df['counts']).astype(float), 1723.5, width=100.5)[1], 258, places=0)
    
    #3
    assert_almost_equal(peak_counts(df.index, df.counts.astype(float), 1723, width=100)[0], 66007, places=0)
    assert_almost_equal(peak_counts(df.index, df.counts.astype(float), 1723, width=100)[1], 258, places=0)
    
    #4
    assert_almost_equal(peak_counts(df.index, df.counts.astype(float), 1376, width=100)[0], 0, places=0)
    assert_almost_equal(peak_counts(df.index, df.counts.astype(float), 1376, width=100)[1], 0, places=0)
    
#-------------------------------------------------------------------------------------------------------------#
def test_foil_count_time():
    """!
    1) Test output given a known results
    2) Test exceptions
    """  
    
    #1
    assert_almost_equal(foil_count_time(0.01, 54000, 548.104260, 0.0151888013272, background=0.01, units='Bq')[0], 1254.519433, places=6)
    assert_almost_equal(foil_count_time(0.01, 16200, 1714.110718, 0.0499603363655, background=0.01, units='Bq')[0], 118.3467643, places=6)
    assert_almost_equal(foil_count_time(0.01, 128160, 46.425931, 0.0150494914458, background=0.01, units='Bq')[0], 17054.945721, places=6)
    assert_almost_equal(foil_count_time(0.01, 128160, 46.425931, 0.0150494914458, background=0.01, units='Bq')[1], 2072.133205, places=6)
    assert_almost_equal(foil_count_time(0.01, 128160, 46.425931, 0.0150494914458, units='Bq')[1], 599.099768, places=6)
    
    #2
    assert_raises(AssertionError,foil_count_time,2, 128160, 46.425931, 0.0150494914458, background=0.01, units='Bq')
    assert_raises(AssertionError,foil_count_time,0.01, -128160, 46.425931, 0.0150494914458, background=0.01, units='Bq')
    assert_raises(AssertionError,foil_count_time,0.01, 128160, -46.425931, 0.0150494914458, background=0.01, units='Bq')
    assert_raises(AssertionError,foil_count_time,0.01, 128160, 46.425931, 1.0150494914458, background=0.01, units='Bq')
    assert_raises(AssertionError,foil_count_time,0.01, 128160, 46.425931, 0.0150494914458, background=-0.01, units='Bq')
    assert_raises(AssertionError,foil_count_time,"One", 128160, 46.425931, 0.0150494914458, background=0.01, units='Bq')
