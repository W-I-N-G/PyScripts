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
def nonparalyzable_dead_time(obsCountRate, tauDetector):
    """!
    @ingroup Detectors
    Calculates the true count rate given a measured count rate and the known
    system dead time according to the nonparalyzable model from Knoll (p.120):

    \f$ n = \frac{m}{1-m \tau} \f$

    @param obsCountRate: <em> integer or float </em> \n
        The recorded count rate for the system in units of [\f$s^{-1}\f$] \n
    @param tauDetector: <em> integer or float </em>  \n
        System dead time in untis of [s]  \n

    @return float: The actual interaction rate \n
            float: The fractional dead time \n
    """

    trueRate = obsCountRate/(1.0-obsCountRate*tauDetector)
    deadTime = (trueRate-obsCountRate)/float(trueRate)
    return trueRate, deadTime

#------------------------------------------------------------------------------#
def paralyzable_dead_time(obsCountRate, tauDetector):
    """!
    @ingroup Detectors
    Calculates the true count rate given a measured count rate and the known
    system dead time according to the paralyzable model from Knoll (p.121):

    \f$ m = n \exp{-n \tau} \f$

    This is solved iteratively since n cannot be solved for exlicitly.

    @param obsCountRate: <em> integer or float </em> \n
        The recorded count rate for the system in units of [\f$s^{-1}\f$] \n
    @param tauDetector: <em> integer or float </em>  \n
        System dead time in untis of [s]  \n

    @return float: The actual interaction rate \n
            float: The fractional dead time \n
    """

    # Create initial guess
    trueRate = nonparalyzableDeadTime(obsCountRate, tauDetector)

    # Solve iteratively
    while abs(trueRate*exp(-trueRate*tauDetector) - obsCountRate) > 1:
        trueRate += 1

    deadTime = (trueRate-obsCountRate)/float(trueRate)
    return trueRate, deadTime

#------------------------------------------------------------------------------#
def nonparalyzable_beam_dead_time(obsCountRate, tauDetector, tauBeam):
    """!
    @ingroup Detectors
    Calculates the true count rate and dead time fraction using a nonparalyzable
    model given a measured count rate and the known system dead time according
    to the bunched synchrotron radiation dead time model from: "Success and
    failure of dead- time models as applied to hybrid pixel detectors in
    high-flux applications."

    \f$ N_{out} = \frac{[1- \exp(-N_{in}* \tau_b]/ \tau_b}
                       {1+[1- \exp(-N_{in}* \tau_b]*n} \f$

    where

    \f$ n=Int( \frac{\tau_s}{\tau_b}) \f$

    This is solved iteratively since \f$N_{in}\f$ cannot be solved for
    exlicitly.

    @param obsCountRate: <em> integer or float </em> \n
        \f$N_{out}\f$: The recorded count rate for the system in units of
        [\f$s^{-1}\f$] \n
    @param tauDetector: <em> integer or float </em>  \n
        \f$\tau_s\f$: System dead time in untis of [s]  \n
    @param tauBeam: <em> integer or float </em>  \n
        \f$\tau_b\f$: The time between beam bunches in untis of [s]  \n

    @return float: \f$N_{in}\f$: The actual interaction rate \n
            float: The fractional dead time \n
    """

    # Create initial guess
    trueRate = obsCountRate
    n = int(tauDetector/float(tauBeam))
    tmpObsRate = (1-exp(-trueRate*tauBeam))/float(tauBeam)\
                 /(1+(1-exp(-trueRate*tauBeam))*n)

    # Solve iteratively
    while abs(tmpObsRate - obsCountRate) > 1:
        trueRate += 1
        tmpObsRate = (1-exp(-trueRate*tauBeam))/float(tauBeam)\
                     /(1+(1-exp(-trueRate*tauBeam))*n)

    deadTime = (trueRate-obsCountRate)/float(trueRate)
    return trueRate, deadTime

#------------------------------------------------------------------------------#
def paralyzable_beam_dead_time(obsCountRate, tauDetector, tauBeam):
    """!
    @ingroup Detectors
    Calculates the true count rate and dead time fraction using a paralyzable
    model given a measured count rate and the known system dead time according
    to the bunched synchrotron radiation dead time model from: "Success and
    failure of dead- time models as applied to hybrid pixel detectors in
    high-flux applications."

    \f$ N_{out} = N_{in}* \exp(-N_{in}* \tau_b*(2n+1)) \f$

    where

    \f$ n=Int(\frac{\tau_s}{\tau_b})\f$

    This is solved iteratively since \f$N_{in}\f$ cannot be solved for
    exlicitly.

    @param obsCountRate: <em> integer or float </em> \n
        \f$N_{out}\f$: The recorded count rate for the system in units of
        [\f$s^{-1}\f$] \n
    @param tauDetector: <em> integer or float </em>  \n
        \f$\tau_s\f$: System dead time in untis of [s]  \n
    @param tauBeam: <em> integer or float </em>  \n
        \f$\tau_b\f$: The time between beam bunches in untis of [s]  \n

    @return float: \f$N_{in}\f$: The actual interaction rate \n
            float: The fractional dead time \n
    """

    # Create initial guess
    trueRate = obsCountRate
    n = int(tauDetector/float(tauBeam))

    # Solve iteratively
    while abs(trueRate*exp(-trueRate*tauBeam*(n+1)) - obsCountRate) > 1:
        trueRate += 1

    deadTime = (trueRate-obsCountRate)/float(trueRate)
    return trueRate, deadTime
