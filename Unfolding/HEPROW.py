"""!
@file Unfolding/HEPROW.py
@package Unfolding

@defgroup HEPROW HEPROW

@brief Routines support data input and output for the Bevins PyScripts package.

@author James Bevins

@date 7Mar17
"""

import numpy as np
import pandas as pd

from scipy.stats import chi2

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
    data = np.sqrt(np.diag(np.asarray(data).astype(float))) \
            * chi2.isf(1 - .68, bins)

    return data
