"""!
@file GeneralNuclear/ActivationAnalysis.py
@package GeneralNuclear

@defgroup ActivationAnalysis ActivationAnalysis

@brief Tools to be used to perform activation analysis.

This module contains data structures and methods to performa activation analysis
and all of the corrections required.  It leverages the basic functions in the
counting module extensively.

@author James Bevins

@date 22Mar17
"""

import os
import peakutils

import pandas as pd
import numpy as np

from math import sqrt

from Counting import parse_spe, get_peak_windows, ge_peakfit

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
    def __init__(self):
        """!
        Constructor to build the ActivationData class. The class is
        insubstantiated with empty attributes.

        @param self: <em> object pointer </em> \n
            The object pointer. \n
        """

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

        ## @var fitParams: <em> list of floats or integers </em>
        # The parameters associated with the specified fit function
        self.fitParams = []

        ## @var eCalibParams: <em> list of floats or integers </em>
        # The parameters associated with the energy calibration
        self.eCalibParams = []

    def __repr__(self):
        """!
        Object print function.

        @param self: <em> ActivationData object pointer </em> \n
            The object pointer. \n
        """

        return "Activation Data Object({}, {}, {}, {})".format(self.rawData,
                                                         self.metaData,
                                                         self.fitFunc,
                                                         self.fitParams)

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
                   .format(self.fitFunc)]
        header += ["with the parameter: {}".format(self.fitParams)]
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
        self.metaData = pd.DataFrame(columns=('Am241Position', 'liveTime',
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
                    df.at[i, 'Am241Position'] = os.path.splitext(filename)\
                                                          [0].split('_')[3]
                    df.at[i, 'line'] = energy[pk]
                    df.at[i, 'realTime'], df.at[i, 'liveTime'] = realTime,\
                                                                 liveTime
                    df.at[i, 'deadTime'] = (realTime-liveTime)/liveTime*100
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

#------------------------------------------------------------------------------#
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
