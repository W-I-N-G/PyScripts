## @file GeneralNuclear/Counting.py
# @package GeneralNuclear
#  
# @defgroup Counting Counting
#
# @brief Tools to be used for calculations involving counting systems.  
#
# This can involve foil activation analysis, detection of unknown sources, etc.  
#
# @author James Bevins
#
# @date 12Feb17

import peakutils

import numpy as np
import pandas as pd

from math import log10, sqrt
from datetime import datetime

#-------------------------------------------------------------------------------------------------------------#
def volume_solid_angle(rSrc, rDet, det2src):    
    """!    
    @ingroup Counting
    Calculates the fractional solid angle for the detector configuration from Knoll pg 119.  
    This is used for a large foil placed close to a detector (but not too close).  

    @param rSrc: <em> integer or float </em>  \n
        Radius of the foil in cm \n
    @param rDet: <em> integer or float </em>  \n
        Radius of the detector in cm \n
    @param det2src: <em> integer or float </em>  \n
        Distance from the detector to src in cm.  This value must be greater than 1 cm. \n

    @return \e float: The gcf for the given configuration \n
    """
    assert det2src>=1.0, "ERROR: The distance between the source and detector must be at least 1.0 cm."
    assert rSrc>=0 and rDet>=0, "ERROR: The radius of the source and/or detector cannot be less than 0.0 cm."
    assert type(rSrc)==int or type(rSrc)==float, "ERROR: The allowed types are integers or floats for the source radius."
    assert type(rDet)==int or type(rDet)==float, "ERROR: The allowed types are integers or floats for the detector radius."
    assert type(det2src)==int or type(det2src)==float, "ERROR: The allowed types are integers or floats for the source to detector distance."
    
    alpha=(rSrc/det2src)**2
    beta=(rDet/det2src)**2
    f1=5./16.*(beta/(1+beta)**(7./2.))-35./64.*(beta**2/(1+beta)**(9./2.))
    f2=35./128.*(beta/(1+beta)**(9./2.))-315./256.*(beta**2/(1+beta)**(11./2.))+1155./1028.*(beta**3/(1+beta)**(13./2.))
    gcf=0.5*(1-1./(1+beta)**(1./2.)-3./8.*(alpha*beta/(1+beta)**(5./2.))+alpha**2*f1-alpha**3*f2)

    return gcf

#-------------------------------------------------------------------------------------------------------------#
def germanium_eff(e,a=0.03279101,b=0.01462466,c=0.15007903,d=-0.0159574):
    """!
    @ingroup Counting
    Calculates the relative efficiency of a germanium detector based on [ref needed]. Defaults are for detector
    #2 in bldg 88 rm 131 at 1 cm.  

    @param e: <em> scalar float/integer or array of floats/integers </em> \n
        Incident gamma ray energy in keV \n
    @param a: \e float \n
        Fit parameter #1  \n
    @param b: \e float \n
        Fit parameter #2  \n
    @param c: \e float \n
        Fit parameter #3  \n
    @param d: \e float \n
        Fit parameter #4  \n
                    
    @return \e float: The relative efficienty for the given configuration and line \n
    """
    
    return a*10-b*10*np.log10(e)+c*0.1*np.log10(e)**2-d*1E4/e**2

#-------------------------------------------------------------------------------------------------------------#
def parse_spe(fname):
    """!
    @ingroup Counting
    Reads in a .Spe spectrum generated from GammaVision or Maestro.  It determines the count time and enegy 
    calibration assuming a linear calibration and places the spectrum into a dataframe.

    @param fname: \e string \n
        The name and path to the .Spe file \n
                    
    @return \e float: The live counting time \n
        \e float: The measurement date and time \n
        \e float: The quadratic term of the energy calibration \n
        \e float: The linear term of the energy calibration \n
        \e float: The intercept of the energy calibration \n
        \e dataframe: The histogram of data for the spectrum \n
    """
    
    try:
        data=pd.read_table(fname,header=0,skipfooter=0)
        data.columns=['counts']

        # Determine energy calibration and live time
        lt=float(data.counts[np.where(data.counts=='$MEAS_TIM:')[0][0]+1].split()[0])
        date=datetime.strptime(data.counts[np.where(data.counts=='$DATE_MEA:')[0][0]+1].split()[0],'%m/%d/%Y')
        time=datetime.strptime(data.counts[np.where(data.counts=='$DATE_MEA:')[0][0]+1].split()[1],'%H:%M:%S')
        c=float(data.counts[np.where(data.counts=='$MCA_CAL:')[0][0]+2].split()[0])
        b=float(data.counts[np.where(data.counts=='$MCA_CAL:')[0][0]+2].split()[1])
        a=float(data.counts[np.where(data.counts=='$MCA_CAL:')[0][0]+2].split()[2])

        # Prune descriptive data
        data=data.drop(data.index[0:np.where(data.counts=='$DATA:')[0][0]+2])
        data=data.drop(data.index[np.where(data.counts=='$ROI:')[0][0]:])

        # Renumber indices
        data.index=range(len(data.counts))
        

        return (lt, datetime.combine(date.date(), time.time()), a, b, c, data)
    
    except IOError: 
        print "WARNING: File does not exist."
        
#-------------------------------------------------------------------------------------------------------------#
def peak_counts(channels, counts, peak, width=25):
    """!
    @ingroup Counting
    Calculate the total number of counts in a peak by assumming a gaussian fit and a linear continuum.  

    @param channels: <em> array/list of integers or floats </em> \n
        The channel index locations \n
    @param counts: <em> array/list of integers or floats </em> \n
        The number of counts in a given channel \n
    @param peak: <em> integer or float </em> \n
        The location of the peak being fitted \n
    @param width: <em> integer or float </em> \n
        The total width of the channel space to be fit.  This should be several times the FWHM resolution \n
                    
    @return \e float: The number of counts in the peak \n
            \e float: The uncertainty of counts in the peak \n
    """
    
    # Fit the peak
    (a,b,c)=peakutils.peak.gaussian_fit(channels[peak-width:peak+width],counts[peak-width:peak+width],center_only=False)
    roiFit=peakutils.peak.gaussian(channels[peak-width:peak+width],a,b,c)
    roiCounts=sum(roiFit) 
    # Determine the continuum to subtract
    roiMap=[0 if item > 0.1 else 1 for item in roiFit]
    baseLine=[]
    for c,r in zip(counts[peak-width:peak+width],roiFit):
        if r<max(roiMap*counts[peak-width:peak+width]):
            baseLine.append(c)
    
    if len(baseLine)>0:
        baseCounts=(len(roiFit)-len(baseLine))*float(sum(baseLine))/len(baseLine)
        return (roiCounts-baseCounts), sqrt(sqrt(roiCounts)**2+sqrt(baseCounts)**2)
    else:
        return (0,0)
    