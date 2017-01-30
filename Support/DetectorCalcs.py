#######################################################################################################
#
# Module : BasicNuclearCalcs.py
#
# Contains : Routines and functions to perform basic detector calculations 
#
# Author : James Bevins
#
# Last Modified: 26Jan17
#
#######################################################################################################

from math import log10

#-------------------------------------------------------------------------------------------------------------#
def volumeGCF(rSrc, rDet, det2src):
    """
    Calculates the solid angle geometry correction factor for the detector configuration from Knoll p. 119.  
    This is used for a large foil placed close to a detector.  
   
    Parameters
    ==========
    rSrc: float
        Radius of the foil in cm
    rDet : float
        Radius of the detector in cm
    det2src: float
        Distance from the detector to src in cm
        
    Returns
    =======
    gcf : float
        The gcf for the given configuration
    """ 

    alpha=(rSrc/det2src)**2
    beta=(rDet/det2src)**2
    f1=5./16.*(beta/(1+beta)**(7./2.))-35./64.*(beta**2/(1+beta)**(9./2.))
    f2=35./128.*(beta/(1+beta)**(9./2.))-315./256.*(beta**2/(1+beta)**(11./2.))+1155./1028.*(beta**3/(1+beta)**(13./2.))
    gcf=0.5*(1-1./(1+beta)**(1./2.)-3./8.*(alpha*beta/(1+beta)**(5./2.))+alpha**2*f1-alpha**3*f2)

    return gcf

#-------------------------------------------------------------------------------------------------------------#
def CalcRelEff(e,det2src):
    """
    Calculates the relative efficiency based on calibration curves at a given distance  
   
    Parameters
    ==========
    e : float
        Incident gamma ray energy in keV
    det2src: float
        Distance from the detector to src in cm 

    Optional
    ========
        
    Returns
    =======
    eff : float
        The relative efficienty for the given configuration and line
    """ 
    
    if det2src==0:
        eff=10**(0.6368*10-0.1175*10*log10(e)+0.9785*0.1*log10(e)**2-(0.3278*10**4)/(e**2))
    elif det2src==1:
        eff=10**(0.5166*10-0.4065*log10(e)-0.4308*0.1*log10(e)**2-0.1888*10**4/e**2)
    elif det2src==2:
        eff=10**(0.5722*10-0.8004*log10(e)+0.1418*0.1*log10(e)**2-0.2112*10**4/e**2)
    elif det2src==3:
        eff=10**(0.585*10-0.9388*log10(e)+0.3163*0.1*log10(e)**2-0.2266*10**4/e**2)
    elif det2src==4:
        eff=10**(0.5736*10-0.9279*log10(e)+0.2328*0.1*log10(e)**2-0.2373*10**4/e**2)
    elif det2src==5:
        eff=10**(0.6114*10-0.1225*10*log10(e)+0.6998*0.1*log10(e)**2-0.2768*10**4/e**2)
    else:
        print "No calibration exists for that detector to source distance"

    return eff