"""!
@file DataAnalysis/DataManipulation.py
@package DataAnalysis

@defgroup DataManipulation DataManipulation

@brief Routines to manipulate lists or arrays of data.

@author James Bevins

@date 17July17
"""

from math import sqrt, log

#------------------------------------------------------------------------------#
def bin_integration(edges, data, edgeLoc="low", lethargy=False):
    """!
    @ingroup DataManipulation

    Integrates binned data.  Can be used, for example, to convert from
    differential flux to flux.  Valid for binned data with edge or midpoint
    values if properly specified in the inputs.

    @param edges: <em> list or array of floats </em> \n
        The lower or upper bin energies.  This list should have a size that is
        one greater than the size of the data. \n
    @param data: <em> list or array of floats </em> \n
        The data corresponding to the bin structure.  This should be of
        len(edges)-1 \n
    @param edgeLoc: \e string \n
        Indicator for the location of the energy boundary edges.  Options
        are "low", "mid", or "up"
    @param lethargy: \e boolean \n
        Specifies whether the integration is done with logarithmic bins. \n

    @return <em> list of floats </em>: The integrated data
    """

    f = []

    edges = check_data(edges, data, edgeLoc)

    # Integrate according to the type of binning provided
    if lethargy:
        for i in range(0, len(data)):
            f.append(log(edges[i+1]/edges[i])*data[i])
    else:
        for i in range(0, len(data)):
            f.append((edges[i+1]-edges[i])*data[i])

    return f

#------------------------------------------------------------------------------#
def bin_differentiation(edges, data, edgeLoc="low", lethargy=False):
    """!
    @ingroup DataManipulation

    Differentiates binned data.  Can be used, for example, to convert from
    flux to differntial flux.  Valid for binned data with edge or midpoint
    values if properly specified in the inputs.

    @param edges: <em> list or array of floats </em> \n
        The lower or upper bin energies.  This list should have a size that is
        one greater than the size of the data. \n
    @param data: <em> list or array of floats </em> \n
        The data corresponding to the bin structure.  This should be of
        len(edges)-1. \n
    @param edgeLoc: \e string \n
        Indicator for the location of the energy boundary edges.  Options
        are "low", "mid", or "up". \n
    @param lethargy: \e boolean \n
        Specifies whether to take the differential as a function of
        ln(dBin). \n

    @return <em> list of floats </em>: The differentiated data
    """

    f = []

    edges = check_data(edges, data, edgeLoc)

    # Differntiate bin by bin
    if edgeLoc == 'low':
        for i in range(0, len(data)):
            if lethargy:
                f.append(data[i]/log(edges[i+1]/float(edges[i])))
            else:
                f.append(data[i]/(edges[i+1]-edges[i]))
    if edgeLoc == 'up':
        if lethargy:
            if edges[0] <=0:
                edges[0] = 1E-10
        for i in range(1, len(data)+1):
            if lethargy:
                f.append(data[i-1]/log(edges[i]/float(edges[i-1])))
            else:
                f.append(data[i-1]/(edges[i]-edges[i-1]))
        
    return f

#------------------------------------------------------------------------------#
def normAUBC(data):
    """!
    @ingroup DataManipulation

    Normalized the area under a binned curve or distribution (AUBC) to equal 1.

    @param data: <em> list or array of floats </em>
        The input data

    @return <em> list or array of floats </em>: The normalized data
    """

    try:
        f = data/sum(data)
    except TypeError:
        sumData = 0
        f = []
        for i in range(0, len(data)):
            sumData += data[i]
        for i in range(0, len(data)):
            f.append(data[i]/sumData)

    return f

#------------------------------------------------------------------------------#
def bin_intensity_to_flux(edges, data, edgeUnits='MeV', edgeLoc="low",
                          detectionVolume=1):
    """!
    @ingroup DataManipulation

    Converts from a neutron intensity to flux according to the formula:
    
    \f$ \phi = nv \f$
    
    \f$ \phi \f$ is the flux, \f$ n \f$ is the neutron density, and \f$ v \f$
    is the neutron velocity.  If the data is not provided in terms of neutron
    density, the detectionVolume optional parameter can be used to normalize
    the data to the correct units. The neutron velocity is calculated from the
    mid-point of each bin.

    @param edges: <em> list or array of floats </em> \n
        The lower or upper bin energies.  This list should have a size that is
        one greater than the size of the data. \n
    @param data: <em> list or array of floats </em> \n
        The data corresponding to the bin structure.  This should be of
        len(edges)-1. \n
    @param edgeUnits: \e string \n
        Specifies the energy units.  Valid options are 'MeV', 'keV', or
        'eV.' \n
    @param edgeLoc: \e string \n
        Indicator for the location of the energy boundary edges.  Options
        are "low", "mid", or "up." \n
    @param detectionVolume: <em> integer or float </em> \n
        Normalization to convert the data to neutron density. \n

    @return <em> list of floats </em>: The flux data.
    """

    assert edgeUnits == "MeV" or edgeUnits == "keV" or edgeUnits == "eV", \
      "Valid specifications for the edge units are 'MeV', 'keV', or 'eV'."

    # Initialize variables
    f = []
    neutronMass = 1.674929E-27
    if edgeUnits == 'MeV':
        energy = 1.602E-13
    elif edgeUnits == 'keV':
        energy = 1.602E-16
    elif edgeUnits == 'eV':
        energy = 1.602E-19

    edges = check_data(edges, data, edgeLoc)

    # Convert to flux
    for i in range(0, len(data)):
        v = sqrt(2*((edges[i+1]+edges[i])/2*energy)/neutronMass)
        f.append(v*data[i]/detectionVolume)

    return f

#------------------------------------------------------------------------------#
def check_data(edges, data, edgeLoc="low"):
    """!
    @ingroup DataManipulation

    Checks the data, converts to lists if necessary, and adds the appropriate
    edges for the calculation.

    @param edges: <em> list or array of floats </em> \n
        The lower or upper bin energies.  This list should have a size that is
        one greater than the size of the data. \n
    @param data: <em> list or array of floats </em> \n
        The data corresponding to the bin structure.  This should be of
        len(edges)-1. \n
    @param edgeLoc: \e string \n
        Indicator for the location of the energy boundary edges.  Options
        are "low", "mid", or "up." \n

    @return <em> list of floats </em>: The modified edges.
    """

    assert edgeLoc == "low" or edgeLoc == "mid" or edgeLoc == "up", \
      "Valid specifications for the edge location are 'low', 'mid', or 'up'."

    f = []

    # Check for expected data consistency
    if edgeLoc == "low":
        if len(edges) == len(data):
            try:
                edges.append(edges[-1]+(edges[-1]-edges[-2]))
            except (TypeError, AttributeError) as e:
                edges = edges.tolist()
                edges.append(edges[-1]+(edges[-1]-edges[-2]))
    if edgeLoc == "up":
        if len(edges) == len(data):
            try:
                edges.insert(0, edges[0]-(edges[1]-edges[0]))
            except (TypeError, AttributeError) as e:
                edges = edges.tolist()
                edges.insert(0, edges[0]-(edges[1]-edges[0]))
    if edgeLoc == "mid":
        tmp = []
        width = edges[1]-edges[0]
        tmp.append(edges[0]-0.5*width)
        tmp.append(edges[0]+0.5*width)
        for i in range(1, len(edges)):
            width = edges[i]-edges[i-1]
            tmp.append(edges[i]+0.5*width)
        edges = tmp
    return edges
