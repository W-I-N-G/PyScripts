"""!
@file GeneralNuclear/ActivationAnalysis.py
@package GeneralNuclear

@defgroup ActivationAnalysis ActivationAnalysis

@brief Tools to be used to perform activation analysis.

This module contains data structures and methods to performa activation analysis
and all of the corrections required.  It leverages the basic functions in the
counting module extensively.

@author James Bevins

@date 23Mar17
"""

import os
import sys
import peakutils

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from math import sqrt
from scipy.optimize import curve_fit
from matplotlib.ticker import MultipleLocator

from Counting import parse_spe, get_peak_windows, ge_peakfit
sys.path.insert(0,os.path.abspath(
 '/home/pyne-user/Dropbox/UCB/Computational_Tools/Scripts/Python/src/DataAnalysis'))
from Math import scaled_exponential

#------------------------------------------------------------------------------#
class ActivationData(object):
    """!
    @ingroup ActivationAnalysis
    This class stores the raw and metadata associated with a counting dataset.
    It does not attempt to store every dataset associated with an overall
    experiment, but is designed to hold and operate on logical grouping of the
    data.  The relevant class functions will depend on the type of data and
    analysisbeing performed.
    """

    ##
    def __init__(self, name=''):
        """!
        Constructor to build the ActivationData class. The class is
        insubstantiated with empty attributes.

        @param self: <em> object pointer </em> \n
            The object pointer. \n
        @param name: \e string \n
            The identifying name of the dataset for plotting and other
            purposes. \n
        """

        ## @var name: \e string
        # An identifying name to use with the dataset
        self.name = name
        
        ## @var rawData: \e dataframe
        # The raw count spectrum data in the form of channels vs counts
        self.rawData = pd.DataFrame()

        ## @var metaData: \e dataframe
        # The general information corresponding to a specific line and
        # measurement. Also contains all of the processed data after
        # after corrections are made.
        self.metaData = pd.DataFrame()

        ## @var fitFunc: \e function
        # The function used to fit the data
        self.fitFunc = None

        ## @var fitParams: <em> array of floats or integers </em>
        # The parameters associated with the specified fit function
        self.fitParams = np.array([0])

        ## @var fitCovar: <em> array of floats or integers </em>
        # The covariance matrix associated with the specified fit function
        # parameters.
        self.fitCovar = np.array([0])

        ## @var eCalibParams: <em> list of floats or integers </em>
        # The parameters associated with the energy calibration
        self.eCalibParams = []

    def __repr__(self):
        """!
        Object print function.

        @param self: <em> ActivationData object pointer </em> \n
            The object pointer. \n
        """

        return "Activation Data Object({}, {}, {}, {}, {}, {})".format(
                        self.rawData, self.metaData, self.fitFunc,
                        self.fitParams, self.fitCovar, self.eCalibParams)

    def __str__(self):
        """!
        Human readable object print function.

        @param self: <em> ActivationData object pointer </em> \n
            The object pointer. \n
        """
        header = ["\nActivation Data Set"]
        header += ["The Raw Data:\n {}".format(self.rawData)]
        header += ["The Metadata:\n {}".format(self.metaData)]
        header += ["The function used to fit the data was {}"\
                   .format(self.fitFunc.__name__)]
        header += ["with the parameters: {}".format(self.fitParams)]
        header += ["and the cavariance: {}".format(self.fitCovar)]
        header += ["The energy calibration parameters: {}"\
                   .format(self.eCalibParams)]
        header = "\n".join(header)+"\n"
        return header

    def read_pileup_data(self, path, cutCh, *args, **kwargs):
        """!
        Reads in pileup data from a file in a given directory. The method
        assumes all of the files are in .spe format and if more than one
        line was desired, an  ActivationData object was provided for each.
        If no additional ActivationData objects are provided, the first
        major line above the cut channel will be used.

        @param self: <em> ActivationData object pointer </em> \n
            The object pointer. \n
        @param path: \e string \n
            The path to the pileup data files. \n
        @param cutCh: \e integer \n
            The channel specifying where above to look for peaks.  This
            should be higher than the channels associated with the high
            activity peak used to cause the dead time. \n
        @param args \n
            An Activation Data object for each line of interest. \n
        @param kwargs \n
            Keyword arguments for the get_peaks function \n
        """

        # Initialize dataframe
        i = 0
        self.metaData = pd.DataFrame(columns=('Cs137Position', 'liveTime',
                               'realTime', 'deadTime', 'counts', 'countUncert',
                               'countRate', 'countRateUncert'))

        # Loop over all files
        for filename in os.listdir(path):
            if filename.endswith(".Spe"):
                name = os.path.splitext(filename)[0].split('_')[3]

                # Read data from this file and store the raw results
                print "Processing:", filename
                (realTime, liveTime, measDate, a, b, c, rawdf) =\
                                                    parse_spe(path+filename)
                if i == 0:
                    self.rawData = rawdf
                    self.rawData.columns = [name]
                else:
                    self.rawData.loc[:, name] = rawdf['counts']

                # Store energy calibrations
                if self.eCalibParams == []:
                    self.eCalibParams = [a, b, c]
                elif self.eCalibParams != [a, b, c]:
                    print "WARNING: Energy calibration parameters are ",\
                          "shifting in the pileup correction data."

                (channels, counts, energy, peaks, windows) = self.get_peaks(
                                                           col=name,
                                                           cutCh=cutCh,
                                                           **kwargs)

                def set_data(df):
                    """!
                    Internal private function for read_pileup_data to set
                    each ofthe cells in the data frame.

                    @param df \n
                        An Activation Data object to store the information
                        in. \n
                    """
                    df.at[i, 'Cs137Position'] = os.path.splitext(filename)\
                                                          [0].split('_')[3]
                    df.at[i, 'line'] = energy[pk]
                    df.at[i, 'realTime'], df.at[i, 'liveTime'] = realTime,\
                                                                 liveTime
                    df.at[i, 'deadTime'] = (realTime-liveTime)/realTime*100
                    pkChannels = channels[windows[pk][0]:windows[pk][1]]
                    pkCounts = counts[windows[pk][0]:windows[pk][1]]
                    (df.at[i, 'counts'], df.at[i, 'countUncert'], redChiSq)\
                                          = ge_peakfit(pkChannels, pkCounts)
                    df.at[i, 'countRate'] = df.at[i, 'counts']\
                                           /df.at[i, 'liveTime']
                    df.at[i, 'countRateUncert'] = sqrt((df.at[i, 'countUncert']\
                                                /df.at[i, 'counts'])**2 \
                                                +(0.5/df.at[i, 'liveTime'])**2)\
                                                *df.at[i, 'countRate']

                # Store each peak in its own dataframe
                for j in range(0, len(peaks)):
                    pk = peaks[j]
                    if j == 0:
                        set_data(self.metaData)
                    elif j <= len(args):
                        set_data(args[j-1].metaData)
                i += 1

        # Store data as floats
        for col in ['liveTime', 'realTime', 'deadTime', 'counts', 'countUncert',
                    'countRate', 'countRateUncert']:
            self.col_to_float('metaData', col)
        self.metaData=self.metaData.sort_values(by='deadTime')

    def get_peaks(self, col, cutCh=0, threshold=0.5, minDist=10):
        """!
        Prepares the rawData for peak analysis.

        @param self: <em> ActivationData object pointer </em> \n
            The object pointer. \n
        @param col: <em> integer or string </em> \n
            The index or name of the raw data set. \n
        @param cutCh: \e integer \n
            The channel specifying where above to look for peaks. \n
        @param threshold: \e integer \n
            The minimum peak threshold.  Calculated from the most prominent
            peak in the spectrum. \n
        @param minDist: \e integer \n
            The minimum number of channels separating peaks. \n

        @return <em> array of floats </em>: The channels for the raw data
             set. \n
        @return <em> array of floats </em>: The counts for the raw data
             set. \n
        @return <em> array of floats </em>: The channels converted to energy
             for the raw data set. \n
        @return <em> list of ints </em>: The peak channel numbers. \n
        @return <em> dictionary of ints </em>: The windows corresponfing to
            each peak. \n
        """

        channels = np.asarray(self.rawData[col].index[cutCh:])
        counts = np.asarray(self.rawData[col][cutCh:])
        energy = self.eCalibParams[2]+self.eCalibParams[1]*channels\
                 +self.eCalibParams[0]*channels**2
        peaks = peakutils.indexes(counts, thres=threshold, min_dist=minDist)
        windows = get_peak_windows(peaks)

        return channels, counts, energy, peaks, windows

    def fit_data(self, func, xName, yName, sigmaName='', dataSet='metaData',
                 **kwargs):
        """!
        Fits a set of data to a specified function and sets the object's 
        attributes associated with the fit.
        
        @param self: <em> ActivationData object pointer </em> \n
            The object pointer. \n
        @param func: \e function \n
           The fitting function to fit the data \n
        @param xName: \e string \n
           The dataframe column name for the x data. \n
        @param yName: \e string \n
           The dataframe column name for the y data. \n
        @param sigmaName: \e string \n
           The dataframe column name for the standard deviation data. \n
        @param dataSet: <em> ActivationData object attribute </em> \n
           The attribute name of the dataset to be fit. Default is to use
           metaData. \n
        @param kwargs \n
            Keyword arguments for the curve fit fitting function.
        """

        assert hasattr(func, '__call__'), 'Invalid function handle'

        if sigmaName == '':
            popt, pcov = curve_fit(func, getattr(self, dataSet)[xName],
                               getattr(self, dataSet)[yName],
                               **kwargs)
        else:
            popt, pcov = curve_fit(func, getattr(self, dataSet)[xName],
                               getattr(self, dataSet)[yName],
                               sigma=getattr(self, dataSet)[sigmaName], 
                               **kwargs)

        # Set the object attributes
        self.fitFunc = func
        self.fitParams = popt
        self.fitCovar = pcov

    def col_to_float(self, dataSet, colName):
        """!
        Converts a column to floats. 
        
        @param self: <em> ActivationData object pointer </em> \n
            The object pointer. \n
        @param dataSet: <em> ActivationData object attribute </em> \n
           The attribute name of th
        @param colName: \e string \n
           The dataframe column name. \ne dataset to be fit. Default is to use
           metaData. \n
        """

        getattr(self, dataSet)[colName] = np.asarray(getattr(self, 
                                            dataSet)[colName], dtype='float64')

    def plot_fitted_data(self, func, xName, yName, *args, **kwargs):
        """!
        Plots a fitted data set and N other data sets.  
        
        @param self: <em> ActivationData object pointer </em> \n
            The object pointer. \n
        @param func: \e function \n
           The fitting function to fit the data \n
        @param xName: \e string \n
           The dataframe column name for the x data. \n
        @param yName: \e string \n
           The dataframe column name for the y data. \n
        @param args: \e datasets \n
            An optional list of additional datasets to plot. The input is the
            name of the AcivationData object. The same x,y, and sigma names
            are used for this dataset.  \n
        @param kwargs: <em> optional plotting inputs </em> \n
            An optional list of additional plot options.  This is wrapped in
            kwargs because 2.7 doesn't support args and keyword specified
            arguements.  The options are listed as kwargs parameters below \n

        @param sigmaName: \e string \n
           The dataframe column name for the standard deviation data. \n
        @param dataSet: <em> ActivationData object attribute </em> \n
           The attribute name of the dataset to be fit. Default is to use
           metaData. \n
        @param logX: <em> kwargs boolean </em> \n
            Flag to use a log scale on the x axis \n
        @param logY: <em> kwargs boolean </em> \n
            Flag to use a log scale on the y axis \n
        @param title: <em> kwargs string </em> \n
            An optional specification for the plot title. \n
        @param xLabel: <em> kwargs string </em> \n
            An optional specification for the x axis label. \n
        @param yLabel: <em> kwargs string </em> \n
            An optional specification for the y axis label. \n
        """

        # Set defaults if not specified since 2.7 sucks
        if 'logX' not in kwargs.keys():
            kwargs['logX'] = False
        if 'logY' not in kwargs.keys():
            kwargs['logY'] = False
        if 'title' not in kwargs.keys():
            kwargs['title'] = ''
        if 'xLabel' not in kwargs.keys():
            kwargs['xLabel'] = ''
        if 'yLabel' not in kwargs.keys():
            kwargs['yLabel'] = ''
        if 'sigmaName' not in kwargs.keys():
            kwargs['sigmaName'] = ''
        if 'dataSet' not in kwargs.keys():
            kwargs['dataSet']='metaData'
            
        # Set up figure
        fig = plt.figure(figsize=(6, 4))
        ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.85])

        # Set axes
        ax1.axis([0.9*min(getattr(self, kwargs['dataSet'])[xName]),
                  1.1*max(getattr(self, kwargs['dataSet'])[xName]),
                  0.9*min(getattr(self, kwargs['dataSet'])[yName]),
                  1.1*max(getattr(self, kwargs['dataSet'])[yName])])
        if kwargs['logX']:
            ax1.set_xscale('log')
        if kwargs['logY']:
            ax1.set_yscale('log')

        # Set axes labels and plot title.
        ax1.set_title('{}'.format(kwargs['title']), fontsize=30, weight="bold")
        ax1.set_xlabel('{}'.format(kwargs['xLabel']), fontsize=22,
                       weight='bold')
        ax1.set_ylabel('{}'.format(kwargs['yLabel']), fontsize=22,
                       weight="bold")
        ax1.tick_params(axis='both', which='major', labelsize=18, width=3)
        ax1.tick_params(axis='both', which='minor', width=3)
        minorLocator = MultipleLocator(1)
        ax1.xaxis.set_minor_locator(minorLocator)

        # Add self to plot
        def add_dataset(df, label, x, y, sig=''):
            """!
            Adds a dataset to the plot.

            @param df: \e dataframe \n
               A dataframe containing the data to plot. \n
            @param label: \e string \n
               A label for the dataframe being plotted. \n
            @param x: \e string \n
               Name of the x data column. \n
            @param y: \e string \n
               Name of the y data column. \n
            @param sig: \e string \n
               Name of the uncertainty column. \n
            """
            if sig != '':
                plt.errorbar(df[x], df[y], yerr=df[sig], label=label, fmt='o')
            else:
                plt.errorbar(df[x], df[y], label=label, fmt='o')
            xFit = np.arange(0.01,1.1*max(df[x]), 0.01)
            yFit = map(lambda y: func(y, *self.fitParams), xFit)
            plt.plot(xFit, yFit)

        # Plot all specified datasets
        add_dataset(getattr(self, kwargs['dataSet']), self.name, xName, yName,
                    kwargs['sigmaName'])
        for arg in args:
            add_dataset(getattr(arg, kwargs['dataSet']), arg.name, xName, yName,
                        kwargs['sigmaName'])

        # Add and locate legend
        plt.legend(borderaxespad=0.75, loc=1, fontsize=16, handlelength=5,
                   borderpad=0.5, labelspacing=0.75, fancybox=True,
                   framealpha=0.5, numpoints=1)

        plt.show()
#------------------------------------------------------------------------------#
