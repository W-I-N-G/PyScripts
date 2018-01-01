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
def stayslFlux(df, fluxName='tally', uncertName='uncertainty',
               maxBinAdjust=0, adjFlux=0.0, adjUncert=0.0):
    """
    @ingroup STAYSL
    Converts a tally dataframe to the STAYSL input format to generate the flux
    and flux uncertainy inputs. The results are printed to the screen.

    There are optional parameters to allow for bin adjustments below a specified
    bin if the simulation records a tally of 0.  This is becuase STAYSL will
    not adjust the bin if a zero flux is recorded in the bin.  This can be
    necessary for difficult simulations or where the low-energy portion of the
    spectrum is not well known. This should be used with EXTREME CAUTION.

    Parameters
    ==========
    @param df: <em> pandas dataframe </em> \n
        Absolute path to the input BCM file. \n
    @param fluxName: \e string \n
        Column header for flux in the supplied dataframe. \n
    @param uncertName: \e string \n
        Column header for flux uncertainty in the supplied dataframe. \n
    @param maxBinAdjust: \e integer \n
        The highest bin number to adjust. \n
    @param adjFlux: \e float \n
        The flux to use for bins with zero flux if below maxBinAdjust. \n
    @param adjUncert: \e float \n
        The uncertainty to use for bins with zero flux if below
        maxBinAdjust. \n
    """

    # Output the flux
    out = ' '
    bin = 0
    print "The flux:"
    for f in df[fluxName]:
        if f == 0 and bin < maxBinAdjust:
            f = adjFlux
        out += '{:.4e} '.format(f)
        if len(out)%78==0:
            print out
            out = ' '
        bin += 1
    print out

    # Output the uncertainty
    out = ' '
    bin = 0
    print "The Uncertainty:"
    for f in df[uncertName]:
        if f == 0 and bin < maxBinAdjust:
            f = adjUncert
        out += '{:.4e} '.format(f)
        if len(out)%78==0:
            print out
            out = ' '
        bin += 1
    print out

#------------------------------------------------------------------------------#
class IterativeSTAYSL(object):
    """!
    @ingroup STAYSL
    This class creates an object used to perform an iterative STAYSL solution.
    The user can specify the convergence criteria, how to handle the solution
    uncertainty, and the convergence tolerance.
    """

    ##
    def __init__(self, path, chiConv=0.1, stdConv=0.1, updateStd=False,
                 **kwargs):
        """!
        Constructor to build the IterativeSTYASL class.

        @param self: <em> object pointer </em> \n
            The object pointer. \n
        @param path: \e string \n
            The path to the base directory containing all of the STAYSL_PNNL
            run files. \n
        @param chiConv: \e integer \n
            The convergence criteria for the $\chi^2$ value. \n
        @param stdConv: \e integer \n
            The convergence criteria for the flux uncertainty. These are
            calculated according to the 2-norm by default, but this can
            be changed using the kwargs. \n
        @param updateStd: \e boolean \n
            Optional specified to update the flux uncertainty with each
            iteration.  If False, the flux uncertainty will not update until
            after the $\chi^2$ convergence criteria has been met.  The
            iteration will then proceed until the flux is converged.  \n
        @param kwargs: <em> optional inputs </em> \n
            An optional list of additional inputs to specify optional inputs
            used by np.linalg.norm(). \n
        """

        ## @var path: \e string
        # The path to the base directory containing all of the STAYSL_PNNL
        # run files.
        self.path = path
        ## @var chiConv: \e integer
        # The convergence criteria for the $\chi^2$ value.
        self.chiConv = chiConv
        ## @var stdConv: \e integer 
        # The convergence criteria for the flux uncertainty.
        self.stdConv = stdConv
        ## updateStd: \e boolean 
        # Optional specified to update the flux uncertainty with each
        # iteration.  
        self.updateStd = updateStd
        
        self.norm = kwargs.pop('norm', 2)

    def __repr__(self):
        """!
        IterativeSTYASL print function.

        @param self: <em> IterativeSTYASL pointer </em> \n
            The IterativeSTYASL pointer. \n
        """
        return "IterativeSTYASL({}, {}, {}, {})".format(self.path,
                                                        self.chiConv,
                                                        self.stdConv,
                                                        self.updateStd)

    def __str__(self):
        """!
        Human readable IterativeSTYASL print function.

        @param self: <em> IterativeSTYASL pointer </em> \n
            The IterativeSTYASL pointer. \n
        """

        header = ["\IterativeSTYASL:"]
        header += ["STAYSL Path: {}".format(self.path)]
        header += ["$\chi^2$ Convergence: {}".format(self.chiConv)]
        header += ["Flux Std  Convergence: {}".format(self.stdConv)]
        header += ["Update Flux Std Each Iteration: {}".format(self.updateStd)]
        header = "\n".join(header)+"\n"
        return header
