## @file GeneralNuclear/BasicNuclearCalcs.py
# @package GeneralNuclear
#
# @defgroup BasicNuclearCalcs BasicNuclearCalcs
#
# @brief Routines and functions to perform basic nuclear calculations 
#
# @author James Bevins
#
# @date 13Feb17

from math import log, pi, sqrt, exp

#-------------------------------------------------------------------------------------------------------------#
def production_decay(halfLife, n, t, rate, src, vol=1, tt=0.0):
    """!
    @ingroup BasicNuclearCalcs
    Calculates the initial population of atoms post-irradiation.  Accounts for production and decay during the 
    irradiation period. It optionally can calculate the number of atoms after a cool down period. Care should be taken with
    the units used for the input parameters as the function can handle both experimental and simulated calculations that use 
    different units and normalizations. \n\n
    
    The two typical sets of units for the production terms are: \n
    
    Eqn:        \f$n(t)=\phi * \Sigma * V * t \f$ \n
    Local Vars: n0= src * rate * vol * t \n
    Units:      [atoms]=[\f$\frac{n}{cm^2 s}\f$][\f$\frac{1}{cm}\f$][\f$cm^3\f$][s] \n\n
    
    Eqn:        \f$n(t)=I * F4Tally (\phi * \Sigma) * V * t \f$ \n
    Local Vars: n0= src * rate * vol * t \n
    Units:      [atoms]=[\f$\frac{n}{cm^2 s}\f$][\f$\frac{1}{cm}\f$][\f$cm^3\f$][s] \n

    @param halfLife: <em> integer or float </em>   \n
        The half life of the decaying isotope in seconds  \n
    @param n: <em> integer or float </em>  \n
        Initial number of parent atoms from previous irradiations  \n
    @param t: <em> integer or float </em> \n
        Irradiation time in seconds \n
    @param rate: <em> integer or float </em>  \n
        The reaction rate. The reaction rate is typically in rx per src per cm$^3$ if volume is specified 
        (for simulated output for example) or cm$^-1$ ($Sigma$) if experimental data is used.  n\n
    @param src: <em> integer or float </em>  \n
        The source strength in n per sec \n
    @param vol: <em> integer or float </em>  \n
        The volume of the foil in cm3 \n  
    @param tt: <em> integer or float </em> \n
        Post-irradiation transfer time in sec \n
        
    @return \e float: The number of atoms after irradiation time t \n
    """ 
    assert halfLife>0, "The half life specified must be greater than zero."
    assert n>=0, "The initial number of atoms specified must be greater than or equal to zero."
    assert t>=0, "The sample decay time specified must be greater than or equal to zero."
    assert rate>=0, "The production rate specified must be greater than or equal to zero."
    assert src>=0, "The src specified must be greater than or equal to zero."
    assert vol>0, "The volume specified must be greater than zero."
    assert tt>=0, "The sample decay time specified must be greater than or equal to zero."
    
    lam=get_decay_const(halfLife)
    # Production and decay during irradiation
    n0=rate*vol*src/lam*(1-exp(-lam*t))+n*exp(-lam*t)
    
    # Decay post irradiation
    n0=n0*exp(-lam*tt)
    
    return n0

#-------------------------------------------------------------------------------------------------------------#
def decay(halfLife, n, t, units='uCi'):
    """!
    @ingroup BasicNuclearCalcs
    Calculates the activity or numnber of atoms following a decay period. Cannot account for ingrowth terms.

    @param halfLife: <em> integer or float </em>   \n
        The half life of the decaying isotope in seconds  \n
    @param n: <em> integer or float </em>  \n
        Initial number of parent atoms  \n
    @param t: <em> integer or float </em>  \n
        The sample decay time \n
    @param units: \e string \n
        If the initial activity is specified, this determines the units provided.  Options are "uCi", "Ci", or "Bq", or "atoms"  \n
        
    @return float: The number of atoms (or activity in Bq) after decay time t \n
    """ 
    assert halfLife>0, "The half life specified must be greater than zero."
    assert n>=0, "The initial number of atoms specified must be greater than or equal to zero."
    assert t>=0, "The sample decay time specified must be greater than or equal to zero."
    
    if units=="uCi":
        n=n*1E-6*3.7E10
    elif units=="Ci":
        n=n*3.7E10
    elif units != "Bq" and units != "atoms":
        print "WARNING: Unknown activity units specified.  Valid specfications are: 'uCi', 'Ci', or 'Bq', or 'atoms'. \
               Decay calculations are not to be trusted."
        
    return n*exp(-get_decay_const(halfLife)*t)

#-------------------------------------------------------------------------------------------------------------#
def get_decay_const(halfLife):
    """!
    @ingroup BasicNuclearCalcs
    Calculates the decay constant given a half life.
   
    @param halfLife: <em> integer or float </em> \n
        The half life of the decaying isotope in seconds \n
        
    @return \e float: The decay constant  \n
    """ 
    assert halfLife>0, "The half life specified must be greater than zero."
    
    return log(2)/halfLife

#-------------------------------------------------------------------------------------------------------------#
def get_halflife(decayConst):
    """!
    @ingroup BasicNuclearCalcs
    Calculates the half life given a decay constant.
    
    @param decayConst: <em> integer or float </em> \n
        The decay constant  \n
        
    @return \e float: The half life of the decaying isotope in seconds \n
    """ 
    assert decayConst>0, "The decay constant specified must be greater than zero."
    
    return log(2)/decayConst

#-------------------------------------------------------------------------------------------------------------#
def activity(halfLife, n, t=0):
    """!
    @ingroup BasicNuclearCalcs
    Calculates the activity of a given isotope with n atoms at time t after production.  An initial activity can
    be specified instead of the initial number of atoms. 
   
    @param halfLife: <em> integer or float </em> \n
        The half life of the decaying isotope in seconds \n
    @param n: <em> integer or float </em> \n
        Initial number of parent atoms \n
    @param t: <em> integer or float </em> \n
        Time post irradiation \n
        
    @return float: The activity at time t in decays/s \n
    """ 
    assert halfLife>0, "The half life specified must be greater than zero."
    assert n>=0, "The initial number of atoms specified must be greater than or equal to zero."
    assert t>=0, "The sample decay time specified must be greater than or equal to zero."
    
    lam=get_decay_const(halfLife)
    
    return lam*n*exp(-lam*t)

#-------------------------------------------------------------------------------------------------------------#
def solid_angle(a, d):
    """!
    @ingroup BasicNuclearCalcs
    Calculates the solid angle assuming a point source (Knoll 4.21).    
   
    @param a: <em> integer or float </em> \n
        Radius of the detector in cm \n
    @param d: <em> integer or float </em> \n
        Distance from the detector to src in cm  \n
        
    @return float: The fractional solid angle for the given configuration \n
    """ 
    assert a>=0, "The detector radius specified must be greater than or equal to zero."
    assert d>=0, "The distance to src specified must be greater than or equal to zero."
    
    return 2.0*pi*(1.0-d/sqrt(d**2+a**2))

#-------------------------------------------------------------------------------------------------------------#
def fractional_solid_angle(a, d):
    """!
    @ingroup BasicNuclearCalcs
    Calculates the fractional solid angle assuming a point source.    
   
    @param a: <em> integer or float </em> \n
        Radius of the detector in cm \n
    @param d: <em> integer or float </em> \n
        Distance from the detector to src in cm  \n
        
    @return float: The fractional solid angle for the given configuration \n
    """ 
    
    return solid_angle(a,d)/4.0/pi