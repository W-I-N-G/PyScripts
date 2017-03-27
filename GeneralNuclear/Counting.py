"""!
@file GeneralNuclear/Counting.py
@package GeneralNuclear

@defgroup Counting Counting

@brief Tools to be used for calculations involving counting systems.

This can involve foil activation analysis, detection of unknown sources, etc.

@author James Bevins

@date 25Feb17
"""

import os
import sys
import peakutils

import numpy as np
import pandas as pd
import copy as cp

from math import sqrt, ceil
from datetime import datetime
from itertools import permutations
from scipy.integrate import quad
from scipy.special import gamma
from scipy.optimize import curve_fit
from BasicNuclearCalcs import activity, decay, fractional_solid_angle

sys.path.insert(0,os.path.abspath(
 '/home/pyne-user/Dropbox/UCB/Computational_Tools/Scripts/Python/DataAnalysis'))
from Math import gauss, smeared_step, skew_gauss, quadratic
from Stats import red_chisq

#------------------------------------------------------------------------------#
def volume_solid_angle(rSrc, rDet, det2src):
    """!
    @ingroup Counting
    Calculates the fractional solid angle for the detector configuration from
    Knoll pg 119. This is used for a large foil placed close to a detector (but
    not too close).

    @param rSrc: <em> integer or float </em>  \n
        Radius of the foil in cm \n
    @param rDet: <em> integer or float </em>  \n
        Radius of the detector in cm \n
    @param det2src: <em> integer or float </em>  \n
        Distance from the detector to src in cm.  This value must be greater
        than 1 cm. \n

    @return \e float: The gcf for the given configuration \n
    """
    assert det2src >= 1.0, "ERROR: The distance between the source and detector\
                          ({} cm) must be at least 1.0 cm.".format(det2src)
    assert rSrc >= 0 and rDet >= 0, "ERROR: The radius of the source and/or\
                          detector cannot be less than 0.0 cm."

    alpha = (rSrc/det2src)**2
    beta = (rDet/det2src)**2
    f1 = 5./16.*(beta/(1+beta)**(7./2.))-35./64.*(beta**2/(1+beta)**(9./2.))
    f2 = 35./128.*(beta/(1+beta)**(9./2.))-315./256.\
         *(beta**2/(1+beta)**(11./2.))+1155./1028.\
         *(beta**3/(1+beta)**(13./2.))
    gcf = 0.5*(1-1./(1+beta)**(1./2.)-3./8.*(alpha*beta/(1+beta)**(5./2.)) \
         +alpha**2*f1-alpha**3*f2)

    return gcf

#------------------------------------------------------------------------------#
def germanium_eff(e, a=0.03279101, b=0.01462466, c=0.15007903, d=-0.0159574):
    """!
    @ingroup Counting
    Calculates the efficiency of a germanium detector based on
    [ref needed]. Defaults are for detector #2 in bldg 88 rm 131 at 1 cm.

    @param e: <em> scalar float/integer or array of floats/integers </em> \n
        Incident gamma ray energy in keV \n
    @param a: \e float \n
        Fit parameter #1 \n
    @param b: \e float \n
        Fit parameter #2 \n
    @param c: \e float \n
        Fit parameter #3 \n
    @param d: \e float \n
        Fit parameter #4 \n

    @return \e float: The relative efficienty for the given configuration and
    line \n
    """

    return a*10-b*10*np.log10(e)+c*0.1*np.log10(e)**2-d*1E4/e**2

#------------------------------------------------------------------------------#
def germanium_eff_exp(e, a=6.00768900e-01, b=5.84842744e-01, c=3.11757094e-11,
                      d=3.76081347e+00):
    """!
    @ingroup Counting
    Calculates the efficiency of a germanium detector based on the four factor
    or exponential formula from

    http://www.ezag.com/fileadmin/ezag/user-uploads/isotopes/pdf/
    Behavior_of_Several_Germanium_Detector_Full_Energy_Peak.pdf

    This fit is not valid for low (~<100 keV) peaks.
    Defaults are for detector #2 in bldg 88 rm 131 at 1 cm.

    @param e: <em> scalar float/integer or array of floats/integers </em> \n
        Incident gamma ray energy in keV \n
    @param a: \e float \n
        Fit parameter #1 \n
    @param b: \e float \n
        Fit parameter #2 \n
    @param c: \e float \n
        Fit parameter #3 \n
    @param d: \e float \n
        Fit parameter #4 \n

    @return \e float: The relative efficienty for the given configuration and
    line \n
    """
    return 1 / (a*e**b + c*e**d)

#------------------------------------------------------------------------------#
def germanium_eff_poly(e, a=-5.86828677e+01, b=5.19051212e+01, \
                       c=-1.81078895e+01, d=3.12451264e+00, f=-2.67044186e-01, \
                       g=9.05096028e-03):
    """!
    @ingroup Counting
    Calculates the efficiency of a germanium detector based on the six factor
    polynomial formula from

    http://www.ezag.com/fileadmin/ezag/user-uploads/isotopes/pdf/
    Behavior_of_Several_Germanium_Detector_Full_Energy_Peak.pdf

    This fit is not valid for low (~<100 keV) peaks.
    Defaults are for detector #2 in bldg 88 rm 131 at 1 cm.

    @param e: <em> scalar float/integer or array of floats/integers </em> \n
        Incident gamma ray energy in keV \n
    @param a: \e float \n
        Fit parameter #1 \n
    @param b: \e float \n
        Fit parameter #2 \n
    @param c: \e float \n
        Fit parameter #3 \n
    @param d: \e float \n
        Fit parameter #4 \n
    @param f: \e float \n
        Fit parameter #5 \n
    @param g: \e float \n
        Fit parameter #6 \n

    @return \e float: The relative efficienty for the given configuration and
    line \n
    """
    eff = a + b*np.log(e) + c*np.log(e)**2 + d*np.log(e)**3 + f*np.log(e)**4 \
             + g*np.log(e)**5
    return eff

#------------------------------------------------------------------------------#
def parse_spe(fname):
    """!
    @ingroup Counting
    Reads in a .Spe spectrum generated from GammaVision or Maestro.  It
    determines the count time and enegy calibration assuming a linear
    calibration and places the spectrum into a dataframe.

    @param fname: \e string \n
        The name and path to the .Spe file \n

    @return \e float: The real counting time \n
        \e float: The live counting time \n
        \e float: The measurement date and time \n
        \e float: The quadratic term of the energy calibration \n
        \e float: The linear term of the energy calibration \n
        \e float: The intercept of the energy calibration \n
        \e dataframe: The histogram of data for the spectrum \n
    """

    try:
        data = pd.read_table(fname, header=0, skipfooter=0)
        data.columns = ['counts']

        # Determine energy calibration and live time
        lt = float(data['counts'][np.where(data['counts'] == \
                                    '$MEAS_TIM:')[0][0]+1].split()[0])
        rt = float(data['counts'][np.where(data['counts'] == \
                                    '$MEAS_TIM:')[0][0]+1].split()[1])
        date = datetime.strptime(data['counts'][np.where(data['counts'] == \
                                '$DATE_MEA:')[0][0]+1].split()[0], '%m/%d/%Y')
        time = datetime.strptime(data['counts'][np.where(data['counts'] == \
                                '$DATE_MEA:')[0][0]+1].split()[1], '%H:%M:%S')
        c = float(data['counts'][np.where(data['counts'] == \
                                          '$MCA_CAL:')[0][0]+2].split()[0])
        b = float(data['counts'][np.where(data['counts'] == \
                                          '$MCA_CAL:')[0][0]+2].split()[1])
        a = float(data['counts'][np.where(data['counts'] == \
                                          '$MCA_CAL:')[0][0]+2].split()[2])

        # Prune descriptive data
        data = data.drop(data.index[0:np.where(data['counts'] == \
                        '$DATA:')[0][0]+2])
        data = data.drop(data.index[np.where(data['counts'] == '$ROI:')[0][0]:])

        # Renumber indices
        data.index = range(len(data['counts']))

        return (rt, lt, datetime.combine(date.date(), time.time()), a, b, c, data.astype(int))

    except IOError:
        print "WARNING: {} does not exist.".format(fname)

#------------------------------------------------------------------------------#
def find_best_fit(*args, **kwargs):
    """!
    @ingroup Counting
    Finds the best fitting function according to \f$\frac{\chi^2}{\nu}\f$

    @param args \n
        N potential fitting functions \n
    @param kwargs \n
        Keyword arguments for the fitting routine scipy curve_fit \n

    @return \e function: The best functional fit \n
    @return \e list: The parameters for the best functional fit \n
    @return \e list: The covariance for the best functional fit \n
    @return \e float: The \f$\frac{\chi^2}{\nu}\f$ of the best functional
        fit \n
    """

    bestChiSq = np.inf
    # Catch when all cases fail
    try:
        for func in args:
            # Test for correct arguments
            try:
                popt, pcov = curve_fit(func, **kwargs)
                yModel = map(lambda y: func(y, *popt), kwargs['xdata'])
                # Catch if user didn't provide sigma values
                try:
                    redChiSq = red_chisq(kwargs['ydata'], yModel,
                                         kwargs['sigma'],
                                         freeParams=len(popt))
                except KeyError:
                    redChiSq = red_chisq(kwargs['ydata'], yModel,
                                         freeParams=len(popt))
            except TypeError:
                print ("ERROR: Either the incorrect minimum arguments or an "
                       "incorrect argument was specified for the ",
                       "scipy.optimize.curve_fit function.")
            except RuntimeError:
                pass
            
            # Save results if a new best fit found
            if redChiSq < bestChiSq:
                bestFunc = func
                params = popt
                covar = pcov
                bestChiSq = redChiSq
        return bestFunc, params, covar, bestChiSq
    except UnboundLocalError:
        print "WARNING: No fit was found."
        return None, [], [], 0.0

#------------------------------------------------------------------------------#
def simple_peak_counts(channels, counts, peak, width=25):
    """!
    @ingroup Counting
    Calculate the total number of counts in a peak by assumming a gaussian fit
    and a linear continuum.

    @param channels: <em> array/list of integers or floats </em> \n
        The channel index locations \n
    @param counts: <em> array/list of integers or floats </em> \n
        The number of counts in a given channel \n
    @param peak: <em> integer or float </em> \n
        The location of the peak being fitted \n
    @param width: <em> integer or float </em> \n
        The total width of the channel space to be fit.  This should be several
        times the FWHM resolution \n

    @return \e float: The number of counts in the peak \n
            \e float: The uncertainty of counts in the peak \n
    """

    # Fit the peak
    (a, b, c) = peakutils.peak.gaussian_fit(channels[peak-width:peak+width], \
                   counts[peak-width:peak+width], center_only=False)
    roiFit = peakutils.peak.gaussian(channels[peak-width:peak+width], a, b, c)
    roiCounts = sum(roiFit)

    # Determine the continuum to subtract
    roiMap = [0 if item > 0.1 else 1 for item in roiFit]
    baseLine = []
    for c, r in zip(counts[peak-width:peak+width], roiFit):
        if r < max(roiMap*counts[peak-width:peak+width]):
            baseLine.append(c)

    if len(baseLine) > 0:
        baseCounts = (len(roiFit)-len(baseLine))*float(sum(baseLine)) \
                     /len(baseLine)
        return (roiCounts-baseCounts), \
                sqrt(sqrt(roiCounts)**2+sqrt(baseCounts)**2)
    else:
        return (0, 0)

#------------------------------------------------------------------------------#
def ge_bincounts(x, p1, p2, p3, p4, p5, p6, p7, p8, p9):
    """!
    @ingroup Counting
    Calculate the total number of counts at a specified bin location given the
    peak fitting parameters.

    "Analaytic Peak Fitting for Gamma-Ray Spectrum Analysis with Ge Detectors"
    by L.C. Longoria

    This formulation dropped the upper exponential and added a quadratic term
    to the background.

    @param x: <em> scalar float/integer </em> \n
        Channel number \n
    @param p1: \e float \n
        Gaussian amplitude \n
    @param p2: \e float \n
        Gaussian centroid \n
    @param p3: \e float \n
        Gaussian width \n
    @param p4: \e float \n
        Skew Gaussian amplitude  \n
    @param p5: \e float \n
        Skew Gaussian range \n
    @param p6: \e float \n
        Smeared step amplitude \n
    @param p7: \e float \n
        Background quadratic term \n
    @param p8: \e float \n
        Background linear term \n
    @param p9: \e float \n
        Background offset \n

    @return \e float: The number of counts in the specified bin
    line \n
    """
    f1 = gauss(x, p1, p2, p3)
    f2 = smeared_step(x, p2, p3, p6)
    f3 = skew_gauss(x, p2, p3, p4, p5)
    f5 = quadratic(x, p7, p8, p9)
    return f1+f2+f3+f5

#------------------------------------------------------------------------------#
def ge_peakcounts(p1, p3, p4, p5):
    """!
    @ingroup Counting
    Calculate the total number of counts at a specified peak given the peak
    fitting parameters using an analytic function for the integration from:

    "Analaytic Peak Fitting for Gamma-Ray Spectrum Analysis with Ge Detectors"
    by L.C. Longoria

    This formulation dropped the upper exponential term.

    @param x: <em> scalar float/integer </em> \n
        Channel number \n
    @param p1: \e float \n
        Gaussian amplitude \n
    @param p3: \e float \n
        Gaussian width \n
    @param p4: \e float \n
        Skew Gaussian amplitude  \n
    @param p5: \e float \n
        Skew Gaussian range \n

    @return \e float: The number of counts in the specified peak \n
    """
    t1 = p1*p3*(2*np.pi)**0.5
    t2 = p3*p4*(gamma(p5)*gamma(4-p5))/6.
    if t2 > t1:
        return t1
    else:
        return t1+t2

#------------------------------------------------------------------------------#
def ge_peakfit(channels, counts, countStd=[], peakWidth=20):
    """!
    @ingroup Counting
    Calculate the total number of counts in a peak. Fits to an the 9 parameter
    function frunction from:

    "Analaytic Peak Fitting for Gamma-Ray Spectrum Analysis with Ge Detectors"
    by L.C. Longoria

    by first assumming a gaussian fit.

    The fits used could be expanded within this framework to test better
    algorithms.  This framework could also be adapted to provide an adaptable
    fit that covers some of the corner cases.

    Currently only handles one peak at a time.

    @param channels: <em> array of integers or floats </em> \n
        The channel index locations for the region of interest to fit \n
    @param counts: <em> array of integers or floats </em> \n
        The number of counts in a given channel \n
    @param countStd: <em> array of integers or floats </em> \n
        The 1\f$\sigma\f$ uncertainty in the counts for each bin \n
    @param countStd: <em> integer or float </em> \n
        The full width at base of the peak in channels \n

    @return \e float: The number of counts in the peak \n
            \e float: The uncertainty of counts in the peak \n
            \e float: \f$\frac{\chi^2}{\nu}\f$ \n
    """

    # Get initial estimate of fitting parameters by simple gaussian fit
    peak = peakutils.indexes(counts, thres=0.25, min_dist=10)[0]
    (a, b, c) = peakutils.peak.gaussian_fit(channels[peak-peakWidth: 
                                                  peak+peakWidth],
                                            counts[peak-peakWidth: 
                                                  peak+peakWidth], center_only=False)
    # Calculate standard deviations if not supplied
    if len(countStd) != len(counts):
        countStd = np.sqrt(np.asarray(counts))
        countStd[countStd == 0] = 1
    # Fit using initial estimates from gaussian fit
    popt, pcov = curve_fit(lambda x, p3, p4, p5, p6, p7, p8, p9: \
                        ge_bincounts(x, a, b, p3, p4, p5, p6, p7, p8, p9),
                        channels, counts, bounds=(0, 5E5), sigma=countStd,
                        absolute_sigma=True,
                        p0=[abs(c), a/20., abs(c)/10., a/100., .01, 1, 10])
    popt = np.insert(popt, 0, b)
    popt = np.insert(popt, 0, a)
    peakCounts = ge_peakcounts(popt[0], popt[2], popt[3], popt[4])
    peakStd = sqrt(peakCounts)

    # Get the bin by bin model data and perform chi squared test
    modelCounts = []
    for ch in channels:
        modelCounts.append(ge_bincounts(ch, *popt))
    redChiSq = red_chisq(counts, modelCounts, countStd, freeParams=9)

    return peakCounts, peakStd, redChiSq

#------------------------------------------------------------------------------#
def get_peak_windows(ch, maxWindow=100, peakWidth=15, minWindow=20):
    """!
    @ingroup Counting
    Find the fitting window for a list of peaks by locating the nearest peak

    @param ch: <em> array/list of integers or floats </em> \n
        The channel index locations \n
    @param maxWindow: \e integer \n
        The half size of the window surrounding the peak. If no other peak is
        found, the window will extend from peak-window : peak+window \n
    @param peakWidth: \e integer \n
        The nominal maximum full width at base of a peak in the spectrum of
        interest \n
    @param minWindow: \e integer \n
        The minimum half window size to be used.  This should be greater than
        peakWidth \n

    @return <em> dictionary of lists </em>: [lower channel, upper channel]
        of the window \n
    """
    windows = {}
    for i in range(0, len(ch)):
        windows[ch[i]] = [0, 0]
        if i == 0 and i != len(ch)-1:
            windows[ch[i]][0] = ch[i] - maxWindow
            if (ch[i] + maxWindow) > (ch[i+1] - peakWidth):
                windows[ch[i]][1] = ch[i] + \
                                      max(ch[i+1]-peakWidth-ch[i], minWindow)
            else:
                windows[ch[i]][1] = ch[i] + maxWindow
        if i == 0 and i == len(ch)-1:
            windows[ch[i]][0] = ch[i] - maxWindow
            windows[ch[i]][1] = ch[i] + maxWindow
        if i != 0 and i != len(ch)-1:
            if (ch[i] - maxWindow) < (ch[i-1] + peakWidth):
                windows[ch[i]][0] = ch[i] - \
                                      max(ch[i]-peakWidth-ch[i-1], minWindow)
            else:
                windows[ch[i]][0] = ch[i] - maxWindow
            if i != 0 and i != len(ch)-1 \
                      and (ch[i] + maxWindow) > (ch[i+1] - peakWidth):
                windows[ch[i]][1] = ch[i] + \
                                       max(ch[i+1]-peakWidth-ch[i], minWindow)
            else:
                windows[ch[i]][1] = ch[i] + maxWindow
        if i == len(ch)-1:
            windows[ch[i]][1] = ch[i] + maxWindow
            if (ch[i] - maxWindow) < (ch[i-1] + peakWidth):
                windows[ch[i]][0] = ch[i] - \
                                       max(ch[i]-peakWidth-ch[i-1], minWindow)
            else:
                windows[ch[i]][0] = ch[i] - maxWindow

    return windows

#------------------------------------------------------------------------------#
def counts(initActivity, halfLife, countTime, units='Bq', countUnits='s'):
    """!
    @ingroup Counting
    Determine the number of counts over a set counting interval assuming no
    background.

    @param initActivity: <em> integer or float </em> \n
        The initial activity.  The default units are Bq but can be changed by
        setting the units flag. \n
    @param halfLife: <em> integer or float </em> \n
        The half life of the decaying isotope in seconds \n
    @param countTime: <em> integer or float </em> \n
        The counting interval in default units of seconds. This can be changed
        by setting the countUnits flag. \n
    @param units: \e string \n
       This determines the units for the activity. Options are "uCi", "Ci",
       or "Bq".  \n
    @param countUnits: \e string \n
       This determines the units provided for the count time. Options are 
       "s", "h", "d", or "y".  \n

    @return \e float: The number of counts \n
    """
    assert halfLife > 0, "The halfLife must be greater than zero."
    assert countTime > 0, "The countTime must be greater than zero."
    assert activity >= 0, "The initial activity must be greater than zero."

    # Get  count time into seconds
    if countUnits == 'h':
        countTime = countTime*3600
    elif countUnits == 'd':
        countTime = countTime*24*3600
    elif countUnits == 'y':
        countTime = countTime*365*24*3600
    elif countUnits != 's':
        print "WARNING: Invalid countUnits specified. Assuming seconds."

    # Get activity in Bq
    if units == "uCi":
        initActivity = initActivity*1E-6*3.7E10
    elif units == "Ci":
        initActivity = initActivity*3.7E10
    elif units != "Bq":
        print "WARNING: Invalid activity units specified. Assuming Bq."

    def integrand(t):
        """ !
        Define the activity integrand.

        @param t: <em> integer or float </em> \n
            The total decay time in seconds \n

        @return \e float: The decays observed \n
        """
        return decay(halfLife, initActivity, t, units='Bq')

    return quad(integrand, 0, countTime)[0]
    
#------------------------------------------------------------------------------#
def foil_count_time(sigma, halfLife, init, efficiency, background=0.001, \
                    units="atoms", precision=30):
    """!
    @ingroup Counting
    Approximate the optimal foil counting time using an average count rate

    @param sigma: \e float \n
        The desired level of relative statistics in fractional form
        (i.e. 1% = 0.01) \n
    @param halfLife: <em> integer or float </em> \n
        The half life of the decaying isotope in seconds \n
    @param init: <em> integer or float </em> \n
        The initial number of atoms.  If the initial activity is provided
        instead, the activity flag must be set to true.  The activity should be
        provided in Bq. \n
    @param efficiency: \e float \n
        The detector efficiency in fractional form (i.e. 1% = 0.01) \n
    @param background: <em> integer or float </em> \n
        The background count rate at the line of interest
    @param units: \e string \n
       This determines the units provided and whether initial atoms or
       activity is provided. Options are "uCi", "Ci", or "Bq", or "atoms"  \n
    @param precision: \e integer \n
       The precision to which to determine the count time. This stops the
       iterative integration once the count time is less than this value. \n

    @return \e float: The number of counts in the peak \n
            \e float: The uncertainty of counts in the peak \n
    """
    assert sigma <= 1, "The relative statistic level must be specified as a \
                        float less than or equal to 1."
    assert halfLife > 0, "The halfLife must be greater than zero."
    assert init >= 0, "The initial atoms/activity must be greater than zero."
    assert efficiency <= 1, "The efficiency must be specified as a float less \
                             than or equal to 1."
    assert background >= 0, "The background must be greater than or equal to \
                             zero."

    def integrand(t):
        """ !
        @ingroup Counting

        Define the activity integrand accounting for all of the efficiencies

        @param t: <em> integer or float </em> \n
            The total decay time in seconds \n

        @return \e float: The decays observed \n

        """
        if units == "atoms":
            return activity(halfLife, init, t)*efficiency
        else:
            return decay(halfLife, init, t, units)*efficiency

    # If the activity is too high, then the dead time will be high; warn user.
    # This assumes 5% dead time on a germanium as determined with a Co60 and 
    # Eu152 src; it is not perfect.
    if units == "atoms" and activity(halfLife, init, 0) > 12000 or \
       units != "atoms" and decay(halfLife, init, 0, units) > 12000:
        print "WARNING: The Dead time may be > 5% with this set-up."

    # Approximate the optimal foil counting time using an average count rate
    tf = 1
    diff = 1000
    try:
        while diff > precision:
            prevt = tf
            s = quad(integrand, 0, tf)[0]/tf
            tf = ((sqrt(s+background)+sqrt(background))**2/(sigma**2*s**2)) \
                 /(1+1/sqrt((s+background)/background))  #Knoll eqn 3.54/55
            diff = tf-prevt
        tb = tf/sqrt((s+background)/background)
        if tf == np.inf:
            tf = 1E99
        return (tf, tb)
    except (ZeroDivisionError, RuntimeWarning):
        return 1E99, 1E99

#------------------------------------------------------------------------------#
def optimal_count_plan(foilParams, handleTime=60, detR=5, background=0.001,
                       units="Bq", toMinute=False,  funcDict={},
                       funcParamDict={}, func=germanium_eff_exp, **kwargs):
    """!
    @ingroup Counting
    Calculates the best order for counting a set of foils by considering all
    possible permutations. The assumption is made that a single foil as named
    will be counted together.  This assumption can be broken simply my giving
    each reaction channel a unique foil name.

    This function assumes that you have a correctly labeled dataframe.  The
    required columns are:

    ['foil','gammaEnergy','halfLife','initActivity','activityUncert',
    'det2FoilDist','relStat','foilR']

    There can be additional columns.

    @param foilParams: \e dataframe \n
        A dataframe containing all of the required data for the count time
        calculation.  NOTE: the columns must be labeled ['foil', 'gammaEnergy',
        'halfLife', 'initActivity', 'activityUncert', 'foilR', 'relStat'] \n
    @param handleTime: <em> integer or float </em> \n
        The estimated handling time for each foil in sec \n
    @param detR: <em> integer or float </em> \n
        The detector radius. Used to correct solid angle for large foils.
        The default values are sufficient if dealing w/ point source foils. \n
    @param background: <em> integer or float </em> \n
        The background count rate at the line of interest \n
    @param units: \e string \n
       This determines the units provided and whether initial atoms or
       activity is provided. Options are "uCi", "Ci", or "Bq"  \n
    @param toMinute: \e boolean \n
       If this flag is set to True, then the count times are rounded to the
       nearest minute. Results returned are still in seconds. \n
    @param funcDict: <em> disctionary of functions </em> \n
       If multiple positions are used, this is a dictionary of the fitting 
       function used vs position where the key is the position name, and the
       value is the fitting function. funcParamDict must be specified if this 
       argument is used. Priority is given to the dictionaries if both methods
       are given as inputs. \n
    @param funcParamDict: <em> disctionary of functions </em> \n
       If multiple positions are used, this is a dictionary of the fitting 
       function parameters vs position where the key is the position name,
       and the value is the N fitting paramters for the function specified in
       funcDict. funcDict must be specified if this argument is used.\n
    @param func: \e function \n
       The effieciency fitting function to calculate the detector absolute
       efficiency. This argument is only valid if there is only one counting
       position being counsidered. Do not specify funcDict and funcParamDict
       if this argument is used. Do specify the kwargs appropriate for this 
       function. \n
    @param kwargs \n
        Keyword arguments for the fitting function. This argument is only
        valid if there is only one counting position being counsidered.

    @return \e dataframe: a copy of the original dataframe with a count time
                         column added (in seconds)  \n
            \e list: a list of the counting order for the foils considered \n
            \e float: the total count time in seconds
    """
    assert hasattr(func, '__call__'), 'Invalid function handle'

    totalTime = np.inf

    # Get activity in Bq
    if units == "uCi":
        foilParams['initActivity'] = foilParams['initActivity']*1E-6*3.7E10
        foilParams['activityUncert'] = foilParams['activityUncert']*1E-6*3.7E10
    elif units == "Ci":
        foilParams['initActivity'] = foilParams['initActivity']*3.7E10
        foilParams['activityUncert'] = foilParams['activityUncert']*3.7E10

    # Consider each possible permutation of the foils
    for order in list(permutations(set(foilParams.foil.tolist()))):
        # Initialize local variables
        df = cp.deepcopy(foilParams)
        df['countTime'] = 0.0
        df['countOrder'] = 0
        df['countActivity'] = cp.deepcopy(df['initActivity'])
        df['countActUncert'] = cp.deepcopy(df['activityUncert'])
        tmpTotal = 0.0

        # Determine count time for each foil in the ordered list
        for f in order:
            ct = 0
            for rx in df.groupby("foil").get_group(f).index:
                pos = df.at[rx, 'det2FoilDist']
                df.at[rx, 'countOrder'] = max(df['countOrder']) + 1

                if funcDict != {} and funcParamDict != {}:
                    absEff = funcDict[pos](df.at[rx, 'gammaEnergy'],
                                           *funcParamDict[pos]) \
                     *(volume_solid_angle(df.at[rx, 'foilR'], detR, pos)) \
                     / fractional_solid_angle(detR, pos) 
                elif funcDict != {} and funcParamDict != {}:
                    print ("WARNING: Both funcDict and funcParamDict were not ",
                          "specified. A single efficieny fit will be used.")
                elif kwargs:
                    absEff = func(df.at[rx, 'gammaEnergy'], **kwargs) \
                     *(volume_solid_angle(df.at[rx, 'foilR'], detR, pos)) \
                     / fractional_solid_angle(detR, pos)
                else:
                    print ("WARNING: Kwargs were not specified for the fitting ",
                          "function. Function defaults will be used, but may ",
                          "not be appropriate.")
                    absEff = func(df.at[rx, 'gammaEnergy']) \
                     *(volume_solid_angle(df.at[rx, 'foilR'], detR, pos)) \
                     / fractional_solid_angle(detR, pos)
                try:
                    if toMinute:
                        df.at[rx, 'countTime'] = ceil(foil_count_time( \
                                          df.at[rx, 'relStat'], \
                                          df.at[rx, 'halfLife'], \
                                          df.at[rx, 'countActivity']- \
                                          3*df.at[rx, 'countActUncert'], \
                                          absEff, background=background, \
                                          units='Bq')[0]/60.)*60
                    else:
                        df.at[rx, 'countTime'] = foil_count_time( \
                                          df.at[rx, 'relStat'], \
                                          df.at[rx, 'halfLife'], \
                                          df.at[rx, 'countActivity']- \
                                          3*df.at[rx, 'countActUncert'], \
                                          absEff, background=background, \
                                          units='Bq')[0]
                except AssertionError:
                    df.at[rx, 'countTime'] = 1E99
                    break

                if df.at[rx, 'countTime'] > ct:
                    ct = df.at[rx, 'countTime']
                df.at[rx, 'countTime'] = max(df['countTime']) + 1

            # Update total counting time for this order
            tmpTotal += ct

            # Update counting times to longest for a given set of reactions
            # within a foil
            for rx in df.groupby("foil").get_group(f).index:
                df.at[rx, 'countTime'] = ct

            # Decay remaining foils by the count time of the current foil
            for rx in df.index:
                if df.at[rx, 'countTime'] == 0.0:
                    df.at[rx, 'countActivity'] = decay(df.at[rx, 'halfLife'], \
                                               df.at[rx, 'countActivity'], \
                                               ct+handleTime, units='Bq')
                    df.at[rx, 'countActUncert'] = decay(df.at[rx, 'halfLife'], \
                                                 df.at[rx, 'countActUncert'], \
                                                 ct+handleTime, units='Bq')

        # Determine if a better solution has been found
        if tmpTotal < totalTime:
            bestOrder = cp.copy(order)
            totalTime = tmpTotal
            bestDF = cp.deepcopy(df)
    
    
    return bestDF.sort_values(by='countOrder'), bestOrder, totalTime

#------------------------------------------------------------------------------#
def channel_statistics(df, countTime, detR=5, units='Bq', 
                       countUnits='s', func=None, **kwargs):
    """!
    @ingroup Counting
    Calculates the counting statistics achieved for all of the channels of a
    given foil assuming a set counting time.  Uses simple integration of the
    activity to get the counts over the period.

    This does take into account the detector efficiency if a efficiency
    function and the matching arguments are provided.  If they are not
    provided, the counts returned will be the total counts into 4PI.

    This function assumes that you have a correctly labeled dataframe.  The
    required columns for use with a detector efficiency function are:

    ['foil','gammaEnergy','halfLife','initActivity','activityUncert',
    'det2FoilDist', 'foilR']

    Otherwise, only the following columns are required:

    ['foil','gammaEnergy','halfLife','initActivity','activityUncert']

    There can be additional columns.

    @param df: \e dataframe \n
        A dataframe containing all of the required data for the count time
        calculation.  NOTE: the columns must be labeled ['foil', 'gammaEnergy',
        'halfLife', 'initActivity', 'activityUncert', 'foilR']. \n
    @param countTime: <em> integer or float </em> \n
        The fixed foil count time. \n
    @param detR: <em> integer or float </em> \n
        The detector radius. Used to correct solid angle for large foils.
        The default values are sufficient if dealing w/ point source foils. \n
    @param units: \e string \n
       This determines the units provided and whether initial atoms or
       activity is provided. Options are "uCi", "Ci", or "Bq"  \n
    @param countUnits: \e string \n
       This determines the units provided for the count time. Options are 
       "s", "h", "d", or "y".  \n
    @param func: \e function \n
       The effieciency fitting function to calculate the detector absolute
       efficiency. This argument is only valid if there is only one counting
       position being counsidered. Do not specify funcDict and funcParamDict
       if this argument is used. Do specify the kwargs appropriate for this 
       function. \n
    @param kwargs \n
        Keyword arguments for the fitting function. This argument is only
        valid if there is only one counting position being counsidered.

    @return \e dataframe: a copy of the original dataframe with a count, 
        count incertainty, and statistics column added (in fractional form) \n
    """

    if func != None:
        assert hasattr(func, '__call__'), 'Invalid function handle'

    # Get  count time into seconds
    if countUnits == 'h':
        countTime = countTime*3600
    elif countUnits == 'd':
        countTime = countTime*24*3600
    elif countUnits == 'y':
        countTime = countTime*365*24*3600
    elif countUnits != 's':
        print "WARNING: Invalid countUnits specified. Assuming seconds."

    # Get activity in Bq
    if units == "uCi":
        df['initActivity'] = df['initActivity']*1E-6*3.7E10
        df['activityUncert'] = df['activityUncert']*1E-6*3.7E10
    elif units == "Ci":
        df['initActivity'] = df['initActivity']*3.7E10
        df['activityUncert'] = df['activityUncert']*3.7E10
    elif units != "Bq":
        print "WARNING: Invalid activity units specified. Assuming Bq."
    
    # Calculate the total counts
    for ind in df.index:
        if func != None:
            df.at[ind, 'counts'] = counts(df.at[ind, 'initActivity'],
                                       df.at[ind, 'halfLife'],
                                       countTime)\
                                *func(df.at[ind, 'gammaEnergy'],
                                      **kwargs)\
                        *(volume_solid_angle(df.at[ind, 'foilR'], detR,
                                         df.at[ind, 'det2FoilDist']))\
                        /fractional_solid_angle(detR, 
                                           df.at[ind, 'det2FoilDist'])
        else:
            df.at[ind, 'counts'] = counts(df.at[ind, 'initActivity'],
                                       df.at[ind, 'halfLife'],
                                       countTime)
        df.at[ind, 'countsUncert'] = df.at[ind, 'activityUncert']\
                                  /df.at[ind, 'initActivity'] \
                                  *df.at[ind, 'counts']
        df.at[ind, 'countingStat'] = np.sqrt(df.at[ind, 'counts'])\
                              /df.at[ind, 'counts']
    return df


