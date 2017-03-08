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

#------------------------------------------------------------------------------#
class histogram(object):
    """!
    @ingroup Histograms
    The class creates a histogram object representing the data.  This defines
    the bin edges and midpoints necessary to manipulate and plot the histogram.
    """

    ##
    def __init__(self, edges=[], data=[], midPtLoc=[], midPtData=[],
                 uncertainty=[]):
        """!
        Constructor to build the calibParams class. If path is specified, the
        object attributes are populated.

        @param self: <em> object pointer </em> \n
            The object pointer. \n
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
        return "ADVANTG Settings({}, {}, {}, {}, {})".format(self.xEdges,
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

    def buildHisto(self, edges, data, uncert=[], edgeLoc="low"):
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
              
    def plotHisto(self, **kwargs):
        """!
        Builds a histogram object from input tabulated data.

        @param self: <em> histogram pointer </em> \n
            The histogram pointer. \n
        @param args: <em> list of histograms </em> \n
            An optional list of additional histograms to plot where the
            keyword will be used for the data set label. \n
        """   
        
        # Allow use of Tex sybols
        plt.rc('text', usetex=True)

        # Set up figure
        #fig = plt.figure()
        fig = plt.figure(figsize=(10,6))
        #ax1 = fig.add_subplot(111)
        ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.85])

        # Preset data set format scheme
        s=10
        linewidth=['2.5']
        marker=['.','o','v','^','<','>','1','2','3','4','8','s','p','*','h','H','+','x','d','D']
        linestyle=['-','--',':','-.']
        dashes=[[10, 0.1],[10, 5, 10, 5],[10,2.5,1,2.5]]
        minorLocator = MultipleLocator(1)

        # Set Line color cycle
        ax1.set_color_cycle(['k', 'k', 'k', 'k'])

        # Set axes
        ax1.axis([0, 25, 0.0001, 1.5*max(heprowHisto.yValues)])
        #ax1.set_xscale('log')
        ax1.set_yscale('log')

        # Set axes labels and plot title.
        ax1.set_title('\\textbf{16MeV D Breakup on Ta}', fontsize=18, weight="bold")    
        ax1.set_xlabel('\\textbf{Energy [MeV]}', fontsize=18, weight="bold")
        ax1.set_ylabel('\\textbf{Neutron PDF}', fontsize=18, weight="bold")
        ax1.tick_params(axis='both', which='major', labelsize=18, width=2)
        ax1.tick_params(axis='both', which='minor', width=2)
        ax1.xaxis.set_minor_locator(minorLocator)

        # Add data set to plot
        ax1.errorbar(heprowHisto.midPtX, heprowHisto.midPtY, yerr=heprowHisto.sigma, marker=None, linestyle='None')
        ax1.plot(heprowHisto.xEdges, heprowHisto.yValues, linewidth=linewidth[0], linestyle=linestyle[0], 
                 marker=None,label="HEPROW Unfolded", dashes=dashes[0]) 
        ax1.errorbar(nsdHisto.midPtX, nsdHisto.midPtY, yerr=nsdHisto.sigma, marker=None, linestyle='None')
        ax1.plot(nsdHisto.xEdges, nsdHisto.yValues, linewidth=linewidth[0], linestyle=linestyle[1], 
                 marker=None,label="NSD Unfolded", dashes=dashes[1]) 
        ax1.plot(meuldersHisto.xEdges, meuldersHisto.yValues, linewidth=linewidth[0], linestyle=linestyle[2], 
                 label="Meulders") 
        #ax1.plot(loneHisto.xEdges, loneHisto.yValues, linewidth=linewidth[0], linestyle=linestyle[3], 
        #         label="Lone", dashes=dashes[2]) 


        # Add and locate legend
        leg = ax1.legend()
        plt.legend(borderaxespad=0.75, loc=1, fontsize=16, handlelength=5, borderpad=0.5,\
                    labelspacing=0.75, fancybox=True, framealpha=0.5, numpoints=1);

        plt.show()
