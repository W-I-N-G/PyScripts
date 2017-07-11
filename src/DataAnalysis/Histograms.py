"""!
@file DataAnalysis/Histograms.py
@package DataAnalysis

@defgroup Histograms Histograms

@brief A histogram class and related support functions to handle histograms of
       data

@author James Bevins

@date 14Jun17
"""

import os

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import MultipleLocator
from matplotlib.colors import LogNorm

#------------------------------------------------------------------------------#
class Histogram(object):
    """!
    @ingroup Histograms
    The class creates a histogram object representing the data.  This defines
    the bin edges and midpoints necessary to manipulate and plot the histogram.
    """

    ##
    def __init__(self, name='', edges=None, data=None, midPtLoc=None,
                 midPtData=None, uncertainty=None):
        """!
        Constructor to build the Histograms class.

        @param self: <em> object pointer </em> \n
            The object pointer. \n
        @param name: \e string \n
            An identifier for the histogram. \n
        @param edges: <em> list or array of integers or floats </em> \n
            The upper and lower edges for each bin.  If there are N bins,
            there should be 2N+1 edges. \n
        @param data: <em> list or array of integers or floats </em> \n
            The y values at the upper and lower edges for each bin. The number
            of data points should equal the number of edges. \n
        @param midPtLoc: <em> list or array of integers or floats </em> \n
            The mid point locations for each bin.  If there are N bins, there
            should be N midpoint x locations. \n
        @param midPtData: <em> list or array of integers or floats </em> \n
            The mid point values for each bin.  If there are N bins, there
            should be N midpoint y values. \n
        @param uncertainty: <em> list or array of integers or floats </em> \n
            The mid point uncertainties for each bin.  If there are N bins,
            there should be N midpoint uncertainties. \n
        """

        ## @var label: \e string
        # A string identifier for the class instance.
        self.label = name
        ## @var xEdges: <em> list or array of integers or floats </em>
        # The upper and lower edges for each bin.  If there are N bins, there
        # will be 2N+1 edges.
        self.xEdges = edges
        ## @var data: <em> list or array of integers or floats </em>
        # The independent variable values at the upper and lower edges for
        # each bin.  The number of values should equal the number of edges.
        self.data = data
        ## @var midPtX: <em> list or array of integers or floats </em>
        # The mid point locations for each bin.  If there are N bins, there
        # should be N midpoint x locations.
        self.midPtX = midPtLoc
        ## @var midPtData: <em> list or array of integers or floats </em>
        # The mid point values for each bin.  If there are N bins, there
        # should be N midpoint values.
        self.midPtData = midPtData
        ## @var sigma: <em> list or array of integers or floats </em>
        # The mid point uncertainties for each bin.  If there are N bins,
        # there should be N midpoint uncertainties.
        self.sigma = uncertainty

    def __repr__(self):
        """!
        Histogram print function.

        @param self: <em> histogram pointer </em> \n
            The histogram pointer. \n
        """
        return "Histogram({}, {}, {}, {}, {})".format(self.xEdges,
                                                             self.data,
                                                             self.midPtX,
                                                             self.midPtData,
                                                             self.sigma)

    def __str__(self):
        """!
        Human readable histogram print function.

        @param self: <em> histogram pointer </em> \n
            The histogram pointer. \n
        """

        header = ["\nHistogram:"]
        header += ["X        Y"]
        header = "\n".join(header)+"\n"
        tmp = ""
        for i in range(0, len(self.xEdges)):
            tmp += "{0:<7}{1}\n".format(self.xEdges[i], self.data[i])
        header = header + tmp
        header += "\nMid Point Data\n"
        header += "X        Y       Sigma\n"
        tmp = ""
        for i in range(0, len(self.midPtX)):
            tmp += "{0:<7} {1} {2}\n".format(self.midPtX[i], self.midPtData[i],
                                             self.sigma[i])
        header = header + tmp
        return header

    def build_histo(self, edges, data, uncert=None, edgeLoc="low", name=''):
        """!
        Builds a histogram object from input tabulated data.

        @param self: <em> histogram pointer </em> \n
            The histogram pointer. \n
        @param edges: <em> list or array of floats </em> \n
            The lower or upper bin energies.  This list should have a size that
            is one greater than the size of the data. \n
        @param data: <em> list or array of floats </em> \n
            The data corresponding to the bin structure. \n
        @param uncert: <em> list or array of floats </em> \n
            The uncertainty corresponding to the data. \n
        @param edgeLoc: \e string \n
            Indicator for the location of the x boundary edges.  Options
            are "low", "mid", or "up" \n
        @param name: \e string \n
            The name or label associated with the histogram data. \n
        """

        assert edgeLoc == "low" or edgeLoc == "mid" or edgeLoc == "up", \
         "Valid specifications for the edge location are 'low', 'mid', or 'up'."

        # Initialize variables
        self.label = name
        self.sigma = uncert
        self.xEdges = []
        self.data = []
        self.midPtX = []
        self.midPtData = []

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
            for i in range(1, len(edges)):
                width = edges[i]-edges[i-1]
                tmp.append(edges[i]-0.5*width)
            tmp.append(edges[-1]+0.5*width)
            edges = tmp

        # Build histogram data
        for i in range(1, len(data)+1):
            self.xEdges.append(edges[i-1])
            self.data.append(data[i-1])
            self.xEdges.append(edges[i])
            self.data.append(data[i-1])
            self.midPtX.append((edges[i-1]+edges[i])/2)
            self.midPtData.append(data[i-1])

        self.xEdges.append(edges[-1])
        self.data.append(0.0)

    def plot(self, *args, **kwargs):
        """!
        Plots a histogram object.

        @param self: <em> histogram pointer </em> \n
            The histogram pointer. \n
        @param args: \e histograms \n
            An optional list of additional histograms to plot. \n
        @param kwargs: <em> optional plotting inputs </em> \n
            An optional list of additional plot MatPlotLib plotoptions.  
            This is wrapped in kwargs because 2.7 doesn't support args
            and keyword specified arguments. The options supported are
            listed as kwargs parameters below. \n
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
        @param savePath: <em> kwargs string </em> \n
            An optional specification for the save location. \n
        @param xMin: <em> kwargs integer or float </em> \n
            An optional specification for the minimum X axis value. \n
        @param xMax: <em> kwargs integer or float </em> \n
            An optional specification for the maximum X axis value. \n
        @param yMin: <em> kwargs integer or float </em> \n
            An optional specification for the minimum Y axis value. \n
        @param yMax: <em> kwargs integer or float </em> \n
            An optional specification for the maximum Y axis value. \n
        @param figSize: <em> kwargs tuple </em> \n
            The (x,y) scale of the plot. \n
        @param grid: <em> kwargs boolean </em> \n
            Specifies whether to plot the tick grid. \n
        @param legendLoc: <em> kwargs int </em> \n
            Specifies the legend location. \n
        @param legend: <em> kwargs int </em> \n
            Specifies wehter to include the legend. \n
        """

        # Set defaults if not specified since 2.7 sucks
        logX = kwargs.pop('logX', False)
        logY = kwargs.pop('logY', False)
        title = kwargs.pop('title', '')
        xLabel = kwargs.pop('xLabel', '')
        yLabel = kwargs.pop('yLabel', '')
        savePath = kwargs.pop('savePath', '')
        xMin = kwargs.pop('xMin', 0)
        xMax = kwargs.pop('xMax', max(self.xEdges)+1)
        yMin = kwargs.pop('yMin', 0.5*min(y for y in self.data if y > 0))
        yMax = kwargs.pop('yMax', 1.5*max(self.data))
        figSize = kwargs.pop('figSize', (12, 8))
        grid = kwargs.pop('grid', True)
        legendLoc = kwargs.pop('legendLoc', 1)
        legend = kwargs.pop('legend', True)

        # Allow use of Tex sybols
        plt.rc('text', usetex=True)

        # Set up figure
        fig = plt.figure(figsize=figSize)
        ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.85])

        # Preset data set format scheme
        linewidth = [2, 4]
        linestyle = ['-', ':', '-.', '--', '-', '-']
        dashes = [[10, 0.1], [2, 2, 2, 2], [10, 5, 2, 5], 
                  [10, 5, 10, 5], [10, 2, 2, 2, 2, 2],
                  [10, 2, 10, 2, 2, 2, 2, 2]]
        ax1.set_prop_cycle(color=['k', 'k', 'k', 'k', 'k', 'k'])

        # Set axes
        ax1.axis([xMin, xMax, yMin, yMax])
        if logX:
            ax1.set_xscale('log')
        if logY:
            ax1.set_yscale('log')

        # Set axes labels and plot title.
        ax1.set_title('{}'.format(title), fontsize=30, weight="bold")
        ax1.set_xlabel('{}'.format(xLabel), fontsize=20, weight='bold',
                       labelpad=10)
        ax1.set_ylabel('{}'.format(yLabel), fontsize=20, weight="bold",
                       labelpad=10)
        
        # Set minor and major gridlines 
        #minorLocator = MultipleLocator(1)
        #ax1.xaxis.set_minor_locator(minorLocator)
        ax1.xaxis.set_tick_params(which='major', width=2, labelsize=20, length=5)
        ax1.yaxis.set_tick_params(which='major', width=2, labelsize=20, length=5)
        ax1.xaxis.set_tick_params(which='minor', width=1.5, length=4)
        ax1.yaxis.set_tick_params(which='minor', width=1.5, length=4)
        ax1.xaxis.grid(b=True, which='major', color='0.2', linestyle='-', alpha=0.5)
        ax1.yaxis.grid(b=True, which='both', color='0.2', linestyle='-', alpha=0.5)

        # Add self to plot
        ax1.plot(self.xEdges, self.data, linewidth=linewidth[0],
                 linestyle=linestyle[0], dashes=dashes[0], marker=None,
                 label=self.label)
        if self.sigma != None:
            ax1.errorbar(self.midPtX, self.midPtData, yerr=self.sigma,
                         marker=None, linestyle='None', capsize=4, capthick=1.5)

        # Plot additional histograms, if specified
        num = 1
        for arg in args:
            ax1.plot(arg.xEdges, arg.data, linewidth=linewidth[num//6],
                     linestyle=linestyle[num%6], dashes=dashes[num%6],
                     marker=None, label=arg.label)
            if arg.sigma != []:
                ax1.errorbar(arg.midPtX, arg.midPtData, yerr=arg.sigma,
                             marker=None, linestyle='None', capsize=4,
                             capthick=1.5)
            num += 1

        # Add and locate legend
        if self.label != '' and legend == True:
            plt.legend(borderaxespad=0.75, loc=legendLoc, fontsize=18,
                       handlelength=5, borderpad=0.5, labelspacing=0.75,
                       fancybox=True, framealpha=0.75, numpoints=1)

        plt.show()

        if savePath != '':
            fig.savefig(savePath, bbox_inches='tight')

    def write(self, path, includeUncert=False, edge=False):
        """!
        Writes a histogram object to a txt file.

        @param self: <em> histogram pointer </em> \n
            The histogram pointer. \n
        @param path: <em> kwargs string </em> \n
            Specification for the save location. Must include both path and
            name. \n
        @param includeUncert: \e string \n
            An option to include the uncertainty in the output. \n
        @param edge: \e string \n
            An option to output the edges or midpoint values. \n
        """

        if os.path.exists(path):
            if os.path.isfile("{}.txt".format(path)):
                os.remove("{}.txt".format(path))

        # Create and open input file
        try:
            with open("{}.txt".format(path), "w") as inpFile:
                if edge == False:
                    for i in range(len(self.midPtX)):
                        if includeUncert == True:
                            inpFile.write("{} {} {}\n".format(self.midPtX[i],
                                                              self.midPtData[i],
                                                              self.sigma[i]))
                        else:
                            inpFile.write("{} {}\n".format(self.midPtX[i],
                                                           self.midPtData[i]))
                else:
                    for i in range(len(self.xEdges)):
                        if includeUncert == True:
                            inpFile.write("{} {} {}\n".format(self.xEdges[i],
                                                             self.data[i],
                                                             self.sigma[i/2]))
                        else:
                            inpFile.write("{} {}\n".format(self.xEdges[i],
                                                           self.data[i]))

            # Close the file
            inpFile.close()

        except IOError as e:
            print "I/O error({0}): {1}".format(e.errno, e.strerror)
            print "File not found was: {0}".format(path)

        # Test that the file closed
        assert inpFile.closed == True, "File did not close properly."

#------------------------------------------------------------------------------#
class Histogram2D(Histogram):
    """!
    @ingroup Histograms
    The class creates a 2Dhistogram object representing the data.  This defines
    the bin edges and midpoints necessary to manipulate and plot the histogram.
    """

    ##
    def __init__(self, name='', xEdges=None, yEdges=None, data=None,
                 xMidPtLoc=None, yMidPtLoc=None, midPtData=None,
                 uncertainty=None):
        """!
        Constructor to build the Histograms2d class.

        @param self: <em> object pointer </em> \n
            The object pointer. \n
        @param name: \e string \n
            An identifier for the histogram. \n
        @param xEdges: <em> list or array of integers or floats </em> \n
            The upper and lower x edges for each bin.  If there are N bins,
            there should be 2N+1 edges. \n
        @param yEdges: <em> list or array of integers or floats </em> \n
            The upper and lower y edges for each bin.  If there are N bins,
            there should be 2N+1 edges. \n
        @param data: <em> list or array of integers or floats </em> \n
            The z values at the upper and lower edges for each bin. The number
            of data points should equal the number of edges. \n
        @param xMidPtLoc: <em> list or array of integers or floats </em> \n
            The x mid point locations for each bin.  If there are N bins, there
            should be N midpoint x locations. \n
        @param yMidPtLoc: <em> list or array of integers or floats </em> \n
            The y mid point locations for each bin.  If there are N bins, there
            should be N midpoint y locations. \n
        @param midPtData: <em> list or array of integers or floats </em> \n
            The mid point values for each bin.  If there are N bins, there
            should be N midpoint z values. \n
        @param uncertainty: <em> list or array of integers or floats </em> \n
            The mid point uncertainties for each bin.  If there are N bins,
            there should be N midpoint uncertainties. \n
        """

        Histogram.__init__(self, name, xEdges, data, xMidPtLoc, midPtData,
                           uncertainty)

        ## @var yEdges: <em> list or array of integers or floats </em>
        # The upper and lower edges for each bin.  If there are N bins, there
        # will be 2N+1 edges.
        self.yEdges = yEdges
        ## @var midPtY: <em> list or array of integers or floats </em>
        # The mid point locations for each bin.  If there are N bins, there
        # should be N midpoint y locations.
        self.midPtY = yMidPtLoc

    def __repr__(self):
        """!
        2D Histogram print function.

        @param self: <em> histogram2d pointer </em> \n
            The histogram2d pointer. \n
        """
        return "Histogram2D({}, {}, {}, {}, {})".format(self.xEdges,
                                                        self.yEdges,
                                                        self.data,
                                                        self.midPtX,
                                                        self.midPtY,
                                                        self.midPtData,
                                                        self.sigma)

    def __str__(self):
        """!
        Human readable 2D histogram print function.

        @param self: <em> histogram2d pointer </em> \n
            The histogram2d pointer. \n
        """

        header = ["\nHistogram2D:"]
        header += ["X        Y       Z"]
        header = "\n".join(header)+"\n"
        tmp = ""
        for i in range(0, len(self.xEdges)):
            for j in range(0, len(self.yEdges)):
                tmp += "{0:<7}{1:<7}{2}\n".format(self.xEdges[i],
                                                  self.yEdges[j],
                                                  self.data[i, j])
        header = header + tmp
        header += "\nMid Point Data\n"
        header += "X        Y       Z      Sigma\n"
        tmp = ""
        for i in range(0, len(self.midPtX)):
            for j in range(0, len(self.midPtX)):
                tmp += "{0:<7} {1} {2}\n".format(self.midPtX[i],
                                                 self.midPtY[j],
                                                 self.midPtData[i, j],
                                                 self.sigma[i, j])
        header = header + tmp
        return header

    def build_2dHisto(self, xEdges, yEdges, data, uncert=None, xEdgeLoc='low',
                      yEdgeLoc='low', name=''):
        """!
        Builds a histogram object from input tabulated data.

        @param self: <em> histogram pointer </em> \n
            The histogram pointer. \n
        @param xEdges: <em> list or array of floats </em> \n
            The x axis bin edges. \n
        @param yEdges: <em> list or array of floats </em> \n
            The y axis bin edges. \n
        @param data: <em> list or array of floats </em> \n
            The data corresponding to the bin structure. \n
        @param uncert: <em> list or array of floats </em> \n
            The uncertainty corresponding to the data. \n
        @param xEdgeLoc: \e string \n
            Indicator for the location of the x boundary edges.  Options
            are "low", "mid", or "up" \n
        @param yEdgeLoc: \e string \n
            Indicator for the location of the y boundary edges.  Options
            are "low", "mid", or "up" \n
        @param name: \e string \n
            An identifier for the histogram. \n
        """

        assert xEdgeLoc == "low" or xEdgeLoc == "mid" or xEdgeLoc == "up", \
         "Valid specifications for the edge location are 'low', 'mid', or 'up'."
        assert yEdgeLoc == "low" or yEdgeLoc == "mid" or yEdgeLoc == "up", \
         "Valid specifications for the edge location are 'low', 'mid', or 'up'."

        # Initialize variables
        self.label = name
        self.sigma = uncert
        self.midPtX = []
        self.midPtY = []

        # Check for expected data consistency
        if xEdgeLoc == "low":
            if len(xEdges) == len(data):
                try:
                    xEdges.append(xEdges[-1]+(xEdges[-1]-xEdges[-2]))
                except TypeError:
                    xEdges = xEdges.tolist()
                    xEdges.append(xEdges[-1]+(xEdges[-1]-xEdges[-2]))
        if yEdgeLoc == "low":
            if len(yEdges) == len(data[0]):
                try:
                    yEdges.append(yEdges[-1]+(yEdges[-1]-yEdges[-2]))
                except TypeError:
                    yEdges = yEdges.tolist()
                    yEdges.append(xEdges[-1]+(yEdges[-1]-yEdges[-2]))

        if xEdgeLoc == "up":
            if len(xEdges) == len(data):
                try:
                    xEdges.insert(0, 0)
                except TypeError:
                    xEdges = xEdges.tolist()
                    xEdges.insert(0, 0)
        if xEdgeLoc == "up":
            if len(yEdges) == len(data):
                try:
                    yEdges.insert(0, 0)
                except TypeError:
                    yEdges = yEdges.tolist()
                    yEdges.insert(0, 0)

        if xEdgeLoc == "mid":
            tmp = []
            width = xEdges[1]-xEdges[0]
            tmp.append(xEdges[0]-0.5*width)
            for i in range(1, len(xEdges)):
                width = xEdges[i]-xEdges[i-1]
                tmp.append(xEdges[i]-0.5*width)
            tmp.append(xEdges[-1]+0.5*width)
            xEdges = tmp

        if yEdgeLoc == "mid":
            tmp = []
            width = yEdges[1]-yEdges[0]
            tmp.append(yEdges[0]-0.5*width)
            for i in range(1, len(yEdges)):
                width = yEdges[i]-yEdges[i-1]
                tmp.append(yEdges[i]-0.5*width)
            tmp.append(yEdges[-1]+0.5*width)
            yEdges = tmp

        # Set Mid point values
        for i in range(1, len(xEdges)):
            self.midPtX = (xEdges[i] - xEdges[i-1])/2.
        for i in range(1, len(yEdges)):
            self.midPtY = (yEdges[i] - yEdges[i-1])/2.

        # Build histogram data
        self.data = data
        self.midPtData = data
        self.xEdges = xEdges
        self.yEdges = yEdges

    def plot2D(self, **kwargs):
        """!
        Plots a 2D histogram.

        @param self: <em> histogram pointer </em> \n
            The histogram pointer. \n
        @param kwargs: <em> optional plotting inputs </em> \n
            An optional list of additional plot options.  This is wrapped in
            kwargs because 2.7 doesn't support args and keyword specified
            arguements.  The options are listed as kwargs parameters below \n
        @param logX: <em> kwargs boolean </em> \n
            Flag to use a log scale on the x axis \n
        @param logY: <em> kwargs boolean </em> \n
            Flag to use a log scale on the y axis \n
        @param logZ: <em> kwargs boolean </em> \n
            Flag to use a log scale on the z axis \n
        @param title: <em> kwargs string </em> \n
            An optional specification for the plot title. \n
        @param xLabel: <em> kwargs string </em> \n
            An optional specification for the x axis label. \n
        @param yLabel: <em> kwargs string </em> \n
            An optional specification for the y axis label. \n
        @param savePath: <em> kwargs string </em> \n
            An optional specification for the save location. \n
        @param xMin: <em> kwargs integer or float </em> \n
            An optional specification for the minimum X axis value. \n
        @param xMax: <em> kwargs integer or float </em> \n
            An optional specification for the maximum X axis value. \n
        @param yMin: <em> kwargs integer or float </em> \n
            An optional specification for the minimum Y axis value. \n
        @param yMax: <em> kwargs integer or float </em> \n
            An optional specification for the maximum Y axis value. \n
        @param zMin: <em> kwargs integer or float </em> \n
            An optional specification for the minimum Z axis value. \n
        @param zMax: <em> kwargs integer or float </em> \n
            An optional specification for the maximum Z axis value. \n
        @param zIntervals: <em> kwargs integer or float </em> \n
            An optional specification for the number of intervals to be
            used for the Z axis colorscale. \n
        """

        # Set defaults if not specified since 2.7 sucks
        if 'logX' not in kwargs.keys():
            kwargs['logX'] = False
        if 'logY' not in kwargs.keys():
            kwargs['logY'] = False
        if 'logZ' not in kwargs.keys():
            kwargs['logZ'] = False
        if 'title' not in kwargs.keys():
            kwargs['title'] = ''
        if 'xLabel' not in kwargs.keys():
            kwargs['xLabel'] = ''
        if 'yLabel' not in kwargs.keys():
            kwargs['yLabel'] = ''
        if 'savePath' not in kwargs.keys():
            kwargs['savePath'] = ''
        if 'xMin' not in kwargs.keys():
            kwargs['xMin'] = min(self.xEdges)
        if 'xMax' not in kwargs.keys():
            kwargs['xMax'] = max(self.xEdges)
        if 'yMin' not in kwargs.keys():
            kwargs['yMin'] = min(self.yEdges)
        if 'yMax' not in kwargs.keys():
            kwargs['yMax'] = max(self.yEdges)
        if 'zMin' not in kwargs.keys():
            wargs['zMin'] = 1E-6
        if 'zMax' not in kwargs.keys():
            kwargs['zMax'] = np.max(self.data)
        if 'zIntervals' not in kwargs.keys():
            kwargs['zIntervals'] = 8

        # Set up figure
        fig = plt.figure(figsize=(11, 8))
        ax1 = fig.add_axes([0.1, 0.1, 0.7, 0.85])
        axcolor = fig.add_axes([0.85, 0.15, 0.05, 0.70])
        ax1.set_title('{}'.format(kwargs['title']), fontsize=30, weight="bold")
        ax1.set_xlabel('{}'.format(kwargs['xLabel']), fontsize=22,
                       weight='bold')
        ax1.set_ylabel('{}'.format(kwargs['yLabel']), fontsize=22,
                       weight="bold")

        # Set axes and scale
        scale = [kwargs['xMin'], kwargs['xMax'], kwargs['yMin'],
                  kwargs['yMax']]
        ax1.axis(scale)
        if kwargs['logX']:
            ax1.set_xscale('log')
        if kwargs['logY']:
            ax1.set_yscale('log')
        if kwargs['logZ']:
            t = np.logspace(np.floor(np.log10(np.abs(kwargs['zMin']))),
                            np.floor(np.log10(np.abs(kwargs['zMax']))),
                            kwargs['zIntervals'])
        else:
            t = np.linspace(kwargs['zMin'], kwargs['zMax'],
                            kwargs['zIntervals'])

        im = ax1.imshow(self.data, interpolation='none', origin='lower',
                      extent=scale, norm=LogNorm(vmin=kwargs['zMin'],
                                                 vmax=kwargs['zMax']))
        fig.colorbar(im, cax=axcolor, ticks=t, format='$%.2e$')
        plt.show()

        if kwargs['savePath'] != '':
            fig.savefig(kwargs['savePath'], bbox_inches='tight')
