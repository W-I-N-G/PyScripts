"""!
@file GeneralNuclear/Detectors.py
@package GeneralNuclear

@defgroup Detectors Detectors

@brief Tools to be used for calculations involving general detector systems.

Tools and methods that are generally useful in characterizing detector systems.

@author James Bevins

@date 7Mar17
"""

from math import exp

#------------------------------------------------------------------------------#
def nonparalyzableDeadTime(obsCountRate, deadTime):
    """!
    @ingroup Detectors
    Calculates the true count rate given a measured count rate and the known
    system dead time according to the nonparalyzable model from Knoll (p.120):

    \f$ n = \frac{m}{1-m \tau} \f$

    @param obsCountRate: <em> integer or float </em> \n
        The recorded count rate for the system in units of [\f$s^{-1}\f$] \n
    @param deadTime: <em> integer or float </em>  \n
        System dead time in untis of [s]  \n

    @return float: The actual interaction rate \n
    """

    return obsCountRate/(1-obsCountRate*deadTime)

#------------------------------------------------------------------------------#
def paralyzableDeadTime(obsCountRate, deadTime):
    """!
    @ingroup Detectors
    Calculates the true count rate given a measured count rate and the known
    system dead time according to the paralyzable model from Knoll (p.121):

    \f$ m = n \exp{-n \tau} \f$

    This is solved iteratively since n cannot be solved for exlicitly.

    @param obsCountRate: <em> integer or float </em> \n
        The recorded count rate for the system in units of [\f$s^{-1}\f$] \n
    @param deadTime: <em> integer or float </em>  \n
        System dead time in untis of [s]  \n

    @return float: The actual interaction rate \n
    """

    # Create initial guess
    trueRate = nonparalyzableDeadTime(obsCountRate, deadTime)

    # Solve iteratively
    while abs(trueRate*exp(-trueRate*deadTime) - obsCountRate) > 1:
        trueRate += 1

    return trueRate
