"""!
@file DatAnalysis/DataManipulation.py
@package DatAnalysis

@defgroup DataManipulation DataManipulation

@brief Routines to manipulate lists or arrays of data.

@author James Bevins

@date 7Mar17
"""

#------------------------------------------------------------------------------#
def binIntegration(edges, data, edgeLoc="low"):
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

    @return <em> list of floats </em>: The integrated data
    """

    assert edgeLoc == "low" or edgeLoc == "mid" or edgeLoc == "up", \
      "Valid specifications for the edge location are 'low', 'mid', or 'up'."

    f = []

    # Check for expected data consistency
    if edgeLoc == "low":
        if len(edges) == len(data):
            try:
                edges.append(edges[-1]+(edges[-1]-edges[-2]))
            except TypeError:
                edges = edges.tolist()
                edges.append(edges[-1]+(edges[-1]-edges[-2]))
    if edgeLoc == "up":
        if len(edges) == len(data):
            try:
                edges.insert(0, 0)
            except TypeError:
                edges = edges.tolist()
                edges.insert(0, 0)
    if edgeLoc == "mid":
        tmp = []
        width = edges[1]-edges[0]
        tmp.append(edges[0]-0.5*width)
        for i in range(0, len(edges)):
            tmp.append(edges[i]+0.5*width)
        edges = tmp

    # Integrate according to the type of binning provided
    for i in range(0, len(data)):
        f.append((edges[i+1]-edges[i])*data[i])

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
