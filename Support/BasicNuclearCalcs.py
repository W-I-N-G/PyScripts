#######################################################################################################
#
# Module : BasicNuclearCalcs.py
#
# Contains : Routines and functions to perform basic nuclear calculations 
#
# Author : James Bevins
#
# Last Modified: 26Jan17
#
#######################################################################################################

from math import log, pi, sqrt, exp

#-------------------------------------------------------------------------------------------------------------#
def N_0(halfLife, n, t, rate, src, vol, tt=0.0):
    """
    Calculates the initial population of atoms post-irradiation.  Accounts for production and decay during the 
    irradiation period. It optionally can calculate the number of atoms after a cool down period.
   
    Parameters
    ==========
    halfLife : float
        The half life of the decaying isotope in seconds
    n: int
        Initial number of parent atoms from previous irradiations
    t: int
        Irradiation time in seconds
    rate : float
        The reaction rate in rx per sec per src per cm3
    src : float
        The source strength in n per sec    
    vol : float
        The volume of the foil in cm3  

    Optional
    ========
    tt: int
        Post-irradiation transfer time in sec
        
    Returns
    =======
    N_0 float
        The number of atoms after irradiation time t 
    """ 
    
    lam=GetDecayConst(halfLife)
    # Production and decay during irradiation
    N_0=rate*vol*src/lam*(1-exp(-lam*t))+n*exp(-lam*t)
    
    # Decay post irradiation
    N_0=N_0*exp(-lam**tt)
    
    return N_0

#-------------------------------------------------------------------------------------------------------------#
def GetDecayConst(halfLife):
    """
    Calculates the decay constant given a half life.
   
    Parameters
    ==========
    halfLife : float or int
        The half life of the decaying isotope in seconds
        
    Returns
    =======
    decayConst : float
        The decay constant 
    """ 
    
    decayConst=log(2)/float(halfLife)
    
    return decayConst

#-------------------------------------------------------------------------------------------------------------#
def GetHalfLife(decayConst):
    """
    Calculates the half life given a decay constant.
   
    Parameters
    ==========
    decayConst : float or int
        The decay constant 
        
    Returns
    =======
    halfLife : float
        The half life of the decaying isotope in seconds
    """ 
    
    halfLife=log(2)/float(decayConst)
    
    return halfLife

#-------------------------------------------------------------------------------------------------------------#
def Activity(halfLife, n, t=0):
    """
    Calculates the activity of a given isotope at time t after production.  
   
    Parameters
    ==========
    halfLife : float
        The half life of the decaying isotope in seconds
    n: int
        Initial number of parent atoms

    Optional
    ========
    t: int
        Time post irradiation
        
    Returns
    =======
    act : float
        The activity at time t in decays/s
    """ 
    
    lam=GetDecayConst(halfLife)
    act=lam*n*exp(-lam*t)
    
    return act

#-------------------------------------------------------------------------------------------------------------#
def SolidAngle(a, d):
    """
    Calculates the solid angle assuming a point source (Knoll 4.21).    
   
    Parameters
    ==========
    a : float
        Radius of the detector in cm
    d: float
        Distance from the detector to src in cm 
        
    Returns
    =======
    omega : float
        The fractional solid angle for the given configuration
    """ 
    
    omega=2.0*pi*(1.0-d/sqrt(d**2+a**2))
    
    return omega

#-------------------------------------------------------------------------------------------------------------#
def FractionalSolidAngle(a, d):
    """
    Calculates the fractional solid angle assuming a point source.    
   
    Parameters
    ==========
    a : float
        Radius of the detector in cm
    d: float
        Distance from the detector to src in cm 
        
    Returns
    =======
    omega : float
        The fractional solid angle for the given configuration
    """ 
    
    omega=SolidAngle(a,d)/4.0/pi
    
    return omega