# This test suite evaluates all of the corner and edge cases for the functions and
# classes in NuclearStructure.    
#
# @author James Bevins
#
# @date 22Sep17

from GeneralNuclear.NuclearStructure import *

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in
    
#----------------------------------------------------------------------------------------#
def test_atom():    
    """!
    Tests the atom class instantiation.
    
    1) Test that a and z map to appropriate attributes
    2) Test that key work specifiers work, test floating conversion
    3) Test that invalid a and z values are precluded 
    """
    
    #1
    testAtom = Atom(1, 2)
    assert_equal(testAtom.z, 1)
    assert_equal(testAtom.a, 2)
    
    #2
    testAtom = Atom(z=2.0, a=4.0, bindingEnergy=8.0, mass=4.0026, bindingEnergyPerA=2.0,
                    massPerA=1.007)
    assert_equal(testAtom.z, 2)
    assert_equal(testAtom.a, 4)
    assert_equal(testAtom.bindingEnergy, 8.0)
    assert_equal(testAtom.mass, 4.0026)
    assert_equal(testAtom.bindingEnergyPerA, 2.0)
    assert_equal(testAtom.massPerA, 1.007)
    
    #3
    assert_raises(AssertionError,Atom, -1, 1)
    assert_raises(AssertionError,Atom, 1, -1)
    assert_raises(AssertionError,Atom, 6, 5)
    assert_raises(AssertionError,Atom, 121, 240)
    assert_raises(AssertionError,Atom, 120, 440)
    
def test_calcBE():    
    """!
    Tests the binding energy calculator from the atom class.
    
    1) Test edge with Z=A=1
    2) Test known interior
    """
    
    #1
    testAtom = Atom(1, 1)
    assert_equal(testAtom.bindingEnergy, 0)
    assert_equal(testAtom.bindingEnergyPerA, 0)
    
    #2
    testAtom = Atom(52, 125)
    assert_almost_equal(testAtom.bindingEnergy, 1054.468, 3)
    assert_almost_equal(testAtom.bindingEnergyPerA, 8.4357, 4)
    
def test_calcSEMF():    
    """!
    Tests the SEMF calculator from the atom class.
    
    1) Test edge with Z=A=1
    2) Test known interior
    """
    
    #1
    testAtom = Atom(1, 1)
    assert_equal(testAtom.mass, 1.007825)
    assert_equal(testAtom.massPerA, 1.007825)
    
    #2
    testAtom = Atom(52, 125)
    assert_almost_equal(testAtom.mass, 124.9074, 4)
    assert_almost_equal(testAtom.massPerA, 0.999259, 6)
    
def test_calcMinZ():    
    """!
    Tests the MinZ calculator.
    
    1) Test edge with A=1
    2) Test known interior
    """
    
    #1
    assert_equal(calcMinZ(1), 1)
    
    #2
    assert_equal(calcMinZ(125), 53)