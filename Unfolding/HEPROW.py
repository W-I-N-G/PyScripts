"""!
@file Unfolding/HEPROW.py
@package Unfolding

@defgroup HEPROW HEPROW

@brief Routines support data input and output for the Bevins PyScripts package.

@author James Bevins

@date 1May17
"""

import numpy as np
import pandas as pd

#------------------------------------------------------------------------------#
def readGru(path, **kwargs):
    """!
    @ingroup HEPROW
    Reads in a HEPROW .gru output file and returns a pandas data frame
    containing the low bin edges, the absolute flux data, and uncertainty.
    This does not read the correlation coefficient matrix.

    @param path: \e string \n
        Absolute path to the file \n
    @param kwargs: \n
        Keyword arguments for pandas.read_table() \n

    @return <em> pandas data frame </em>: A data frame containing the lower
        bin edges, the absolute flux, and its uncertainty \n
    """

    df = pd.read_table(path, **kwargs)

    # Find the row for he .gru separator for the flux and correlation matrix
    dataStop = "*********format(16i5)*********"
    loc = df[df.ix[:, 0] == dataStop].index.tolist()[0]
    df = df.drop(df.index[loc:])

    return df.apply(pd.to_numeric)

#------------------------------------------------------------------------------#
def readFlu(path, **kwargs):
    """!
    @ingroup HEPROW
    Reads in a HEPROW .fru output file and returns a pandas data frame
    containing the low bin edges, the absolute flux data, and uncertainty.
    This does not read the correlation coefficient matrix.

    @param path: \e string \n
        Absolute path to the file \n
    @param kwargs: \n
        Keyword arguments for pandas.read_table() \n

    @return <em> pandas data frame </em>: A data frame containing the lower
        bin edges, the absolute flux, and its uncertainty \n
    """

    df = pd.read_table(path, **kwargs)

    return df.apply(pd.to_numeric)

#------------------------------------------------------------------------------#
def readMTX(path):
    """!
    @ingroup HEPROW
    Reads in a HEPROW .MTX covariance output file and returns a the sqrt of the
    diagonal times the inverse survival function of the \f$\chi^2\f$
    distribution.

    @param path: \e string \n
        Absolute path to the file \n

    @return <em> array of floats </em>: An array containing the standard
    deviations \n
    """

    data = []

    # Open file
    try:
        f = open(path, 'r')
        # Read header lines
        bins = int(f.next().split()[0])

        # Read the file line by line and store the values
        tmp = []
        for line in f:
            if len(tmp) < bins:
                tmp.extend(line.split())
                if len(tmp) == bins:
                    data.append(tmp)
                    tmp = []
            else:
                print 'ERROR: MTX read error. Covariance length !=  # bins.'

        # Close the file
        f.close()
    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror)

    # Determine the standard deviation
    data = np.sqrt(np.diag(np.asarray(data).astype(float)))

    return data

#------------------------------------------------------------------------------#
def readRSP(path, minE=None, maxE=None, minPH=None, maxPH=None):
    """!
    @ingroup HEPROW
    Reads in a HEPROW .rsp response matrix file and returns an array
    representing the response matrix.

    @param path: \e string \n
        Absolute path to the file \n
    @param minE: <em> integer or float </em> \n
        Optional specifier for the minimum energy to read. \n
    @param maxE: <em> integer or float </em> \n
        Optional specifier for the maximum energy to read. \n
    @param minPH: <em> integer or float </em> \n
        Optional specifier for the minimum pulse height to read. \n
    @param maxPH: <em> integer or float </em> \n
        Optional specifier for the maximum pulse height to read. \n

    @return <em> array of floats </em>: An array containing the standard
    deviations \n
    """

    # Initialize Variables
    data = []
    allE = []
    if minE == None:
        minE = 0
    if maxE == None:
        maxE = 1000

    # Open file
    try:
        f = open(path, 'r')
        # Skip the first header
        line = f.next()
        phScale = float(line.rstrip().split()[0])
        line = f.next()
        curE = float(line.rstrip().split()[0])
        phLowBound = float(line.rstrip().split()[2])
        if minPH == None:
            minPH = phLowBound
        phUpBound = float(line.rstrip().split()[3])
        if maxPH == None:
            maxPH = phUpBound
        count = 1

        # Read the file line by line and store the values
        tmp = []
        for line in f:
            splitList = line.rstrip().split()
            if len(splitList) == 4 and float(splitList[3]) == phUpBound:
                curE = float(line.rstrip().split()[0])
                print curE
                if tmp != []:
                    data.append(tmp)
                    tmp = []
                count = 1
            elif curE >= minE and curE <= maxE:
                if curE not in allE:
                    allE.append(curE)
                for item in splitList:
                    if count*phScale >= minPH and count*phScale <= maxPH:
                        tmp.append(float(item))
                    count += 1

        # Append last data set
        data.append(tmp)

        # Close the file
        f.close()
    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror)

    for row in data:
        while len(row) < len(data[-1]):
            row.append(0.0)

    return np.transpose(np.asarray(data))
