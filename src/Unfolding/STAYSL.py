"""!
@file Unfolding/STAYSL.py
@package Unfolding

@defgroup STAYSL STAYSL

@brief Routines support data input and output for STAYSL.

@author James Bevins

@date 27Dec17
"""

import numpy as np

from datetime import datetime

#------------------------------------------------------------------------------#
def bcmToBCF(bcmPath, outPath='flux_history.dat', timeOut='cumulative', measOut='differential', measType='flux', title='PyScripts Generated BCF file'):
    """
    @ingroup STAYSL
    Converts from a bcm file with the following format 
    
    date time current
    4/26/2017 2:45:16 PM 0.000000014
    
    to a BCF input file format.
    
    Parameters
    ==========
    @param bcmPath: \e string \n
        Absolute path to the input BCM file. \n
    @param outPath: \e string \n
        Absolute path to the created BCF input file. \n
    @param timeOut: \e string \n
        Type of time data output for the BCF input file.  Either 'differential'
        or 'cumulative'. \n
    @param measOut: \e string \n
        Type of measurement output for the BCF input file.  Either 
        'differential' or 'cumulative'. \n
    @param measType: \e string \n
        Type of measurment output in the BCM file.  Either 'flux' or 
        'fluence'. \n
    @param title: \e string \n
        A STAYSL BCF input file title. \n

    """

    # Initialize variables and save header lines
    out = ''
    totTime = 0
    totMeas = 0

    if title.endswith('\n'):
        ntitl = title.count('\n')
    else:
        ntitl = title.count('\n')+1
        title += '\n'
        
    if timeOut.lower() == 'cumulative':
        if measOut.lower() == 'cumulative':
            ntype = 3
        elif measOut.lower() == 'differential':
            ntype = 2
        else:
            print 'ERROR: invalid entry for measOut parameter.'
            return
    elif timeOut.lower() == 'differential':
        if measOut.lower() == 'cumulative':
            ntype = 1
        elif measOut.lower() == 'differential':
            ntype = 0
        else:
            print 'ERROR: invalid entry for measOut parameter.'
            return
    else:
        print 'ERROR: invalid entry for timeOut parameter.'
        return

    if measType.lower() == 'flux':
        mfee = 0
    elif measType.lower() == 'fluence':
        mfee = 1
    else:
        print 'ERROR: Invalid entry for measType paramtere.'
        return
    
    out += '{} {} {} {} {} \n'.format(ntitl, ntype, 0, mfee, 'S')
    out += title
    
    # Open BCM file
    try:
        nrec = 0
        tmpOut = ''
        f = open(bcmPath, 'r')

        # Store first time slice
        line = f.next().rstrip().split('\t')
        nrec += 1
        prevTime = datetime.strptime(line[0]+line[1], '%m/%d/%Y%I:%M:%S %p')

        for line in f:
            line = line.rstrip().split('\t')
            nrec += 1
            curTime = datetime.strptime(line[0]+line[1],
                                        '%m/%d/%Y%I:%M:%S %p')
            curMeas = float(line[2])

            # Store in BCF format
            deltaT = (curTime-prevTime).total_seconds()
            totTime += deltaT
            totMeas += curMeas
            if ntype == 0:
                tmpOut += '     {:.2f}          {:.9f}\n'.format(deltaT, curMeas)
            if ntype == 1:
                tmpOut += '     {:.2f}          {:.9f}\n'.format(deltaT, totMeas)
            if ntype == 2:
                tmpOut += '     {:.2f}          {:.9f}\n'.format(totTime, curMeas)
            if ntype == 3:
                tmpOut += '     {:.2f}          {:.9f}\n'.format(totTime,
                                                                 totMeas)

            prevTime = curTime

        # Close the file
        f.close()
    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror)

    # Finalize the output file and save
    out += '{}\n'.format(nrec)
    out += tmpOut
    out += '     {:.2f}          {:.4f}\n'.format(0.0, -100)
    try:
        nrec = 0
        tmpOut = ''
        f = open(outPath, 'w')
        f.write(out)

        # Close the file
        f.close()
    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror)
    
    # Print output
    print 'The total measurment time was {} seconds with an integrated ' \
          'measurement of {}.'.format(totTime, totMeas)
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
    eBins = []
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
        numPHBins = float(line.rstrip().split()[1])
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
                if tmp != []:
                    data.append(tmp)
                    tmp = []
                count = 1
            elif curE >= minE and curE <= maxE:
                if curE not in eBins:
                    eBins.append(curE)
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

    return np.transpose(np.asarray(data)), eBins, \
           np.linspace(phScale, phUpBound, numPHBins)
