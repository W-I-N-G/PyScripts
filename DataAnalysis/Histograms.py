"""!
@file DataAnalysis/Histograms.py
@package DataAnalysis

@defgroup Histograms Histograms

@brief A histogram class and related support functions to handle histograms of
       data

@author James Bevins

@date 7Mar17
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

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
        Constructor to build the calibParams class. If path is specified, the
        object attributes are populated.

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
        ## @var yValues: <em> list or array of integers or floats </em>
        # The y values at the upper and lower edges for each bin.  The number
        # of yvalues should equal the number of edges.
        self.yValues = data
        ## @var midPtX: <em> list or array of integers or floats </em>
        # The mid point locations for each bin.  If there are N bins, there
        # should be N midpoint x locations.
        self.midPtX = midPtLoc
        ## @var midPtY: <em> list or array of integers or floats </em>
        # The mid point values for each bin.  If there are N bins, there
        # should be N midpoint y values.
        self.midPtY = midPtData
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
                                                             self.yValues,
                                                             self.midPtX,
                                                             self.midPtY,
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
            tmp += "{0:<7}{1}\n".format(self.xEdges[i], self.yValues[i])
        header = header + tmp
        header += "\nMid Point Data\n"
        header += "X        Y       Sigma\n"
        tmp = ""
        for i in range(0, len(self.midPtX)):
            tmp += "{0:<7} {1} {2}\n".format(self.midPtX[i], self.midPtY[i],
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
            Indicator for the location of the energy boundary edges.  Options
            are "low", "mid", or "up" \n
        """

        assert edgeLoc == "low" or edgeLoc == "mid" or edgeLoc == "up", \
         "Valid specifications for the edge location are 'low', 'mid', or 'up'."

        # Initialize variables
        self.label = name
        self.sigma = uncert
        self.xEdges = []
        self.yValues = []
        self.midPtX = []
        self.midPtY = []

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

        # Build histogram data
        for i in range(1, len(data)+1):
            self.xEdges.append(edges[i-1])
            self.yValues.append(data[i-1])
            self.xEdges.append(edges[i])
            self.yValues.append(data[i-1])
            self.midPtX.append((edges[i-1]+edges[i])/2)
            self.midPtY.append(data[i-1])

        self.xEdges.append(edges[-1])
        self.yValues.append(0.0)

    def plot(self, *args, **kwargs):
        """!
        Builds a histogram object from input tabulated data.

        @param self: <em> histogram pointer </em> \n
            The histogram pointer. \n
        @param args: \e histograms \n
            An optional list of additional histograms to plot. \n
        @param kwargs: <em> optional plotting inputs </em> \n
            An optional list of additional plot options.  This is wrapped in
            kwargs because 2.7 doesn't support args and keyword specified
            arguements.  The options are listed as kwargs parameters below \n
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

        # Allow use of Tex sybols
        plt.rc('text', usetex=True)

        # Set up figure
        #fig = plt.figure()
        fig = plt.figure(figsize=(11, 8))
        ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.85])

        # Preset data set format scheme
        linewidth = [2.5]
        linestyle = ['-', '-', ':', '-.']
        dashes = [[10, 0.1], [5, 2, 5, 2], [0.5, 0.5], [10, 2.5, 1, 2.5]]
        ax1.set_prop_cycle(color=['k', 'k'])

        # Set axes
        ax1.axis([0, max(self.xEdges)+1,
                  0.5*min(y for y in self.yValues if y > 0),
                  1.5*max(self.yValues)])
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
        ax1.plot(self.xEdges, self.yValues, linewidth=linewidth[0],
                 linestyle=linestyle[0], dashes=dashes[0], marker=None,
                 label=self.label)
        if self.sigma != None:
            ax1.errorbar(self.midPtX, self.midPtY, yerr=self.sigma, marker=None,
                         linestyle='None', capsize=4, capthick=1.5)

        # Plot additional histograms, if specified
        num = 1
        for arg in args:
            ax1.plot(arg.xEdges, arg.yValues, linewidth=linewidth[0],
                     linestyle=linestyle[num%4], dashes=dashes[num%4],
                     marker=None, label=arg.label)
            if arg.sigma != []:
                ax1.errorbar(arg.midPtX, arg.midPtY, yerr=arg.sigma,
                             marker=None, linestyle='None', capsize=4,
                             capthick=1.5)
            num += 1

        # Add and locate legend
        if self.label != '':
            plt.legend(borderaxespad=0.75, loc=1, fontsize=16, handlelength=5,
                       borderpad=0.5, labelspacing=0.75, fancybox=True,
                       framealpha=0.5, numpoints=1)

        plt.show()
