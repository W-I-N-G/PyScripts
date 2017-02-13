#######################################################################################################
#
# Module : Plotting.py
#
# Contains : Plotting and support routines
#
# Author : James Bevins
#
# Last Modified: 24Jan17
#
#######################################################################################################

import matplotlib.pyplot as plt

#---------------------------------------------------------------------------------------#    
class Histogram:
    """
    Creates a object representing the data required to plot histograms from binned data
   
    Attributes
    ==========
    xEdges : list of integers or floats
        The x edge locations for each bin
        [Deafult: []]
    yValues : list of integers or floats
        The y data values for each bin
        [Deafult: []]
    midPtX : list of integers or floats
        The mid point location for each bin
        [Deafult: []]
    midPtY : list of integers or floats
        The mid point value for each bin
        [Deafult: []]
    uncertainty : list of integers or floats
        The uncertainty for each bin
        [Deafult: []]
    """
        
        
    def __init__(self,xEdges=[],yValues=[],midPtX=[],midPtY=[],uncertainty=[]):  
        self.xEdges=xEdges
        self.yValues=yValues
        self.midPtX=midPtX
        self.midPtY=midPtY
        self.sigma=uncertainty   
        
    def __repr__(self):
        return "ADVANTG Settings({}, {}, {}, {}, {})".format(self.xEdges,self.yValues,self.midPtX,self.midPtY,self.sigma)
    
    def __str__(self):
        header = ["\nHistogram:"]
        header += ["X        Y"]
        header = "\n".join(header)+"\n"
        tmp=""
        for i in range(0,len(self.xEdges)):
            tmp+="{0:<7}{1}\n".format(self.xEdges[i],self.yValues[i])
        header = header + tmp
        header += "\nMid Point Data\n"
        header += "X        Y       Sigma\n"
        tmp=""
        for i in range(0,len(self.midPtX)):
            tmp+="{0:<7} {1} {2}\n".format(self.midPtX[i],self.midPtY[i],self.sigma[i])
        header = header + tmp
        return header

    def buildHisto(self,edges,data,uncert=[],edgeLoc="low"):
        """
        Builds a histogram object from input tabulated data. 

        Parameters
        ==========
        edges : list or array of floats
            The lower or upper bin energies.  This list should have a size that is one greater than the size of the data.
        data: list or array of floats
            The data corresponding to the bin structure. 

        Optional
        ========
        uncert: list or array of floats
            The uncertainty corresponding to the data.  
            [Default=[]]
        edgeLoc : string
            Indicator for the location of the energy boundary edges.  Options are "low", "mid", or "up"  
            [Deafult = "low"]
        """ 

        assert edgeLoc=="low" or edgeLoc=="mid" or edgeLoc=="up", \
          "Valid specifications for the edge location are 'low', 'mid', or 'up'."

        # Initialize variables
        histoEdges=[]
        histoData=[]
        histoMidX=[]
        histoMidY=[]

        # Check for expected data consistency
        if edgeLoc=="low":
            if len(edges)==len(data):
                try: 
                    edges.append(edges[-1]+(edges[-1]-edges[-2]))
                except TypeError:
                    edges=edges.tolist()
                    edges.append(edges[-1]+(edges[-1]-edges[-2]))
        if edgeLoc=="up":
            if len(edges)==len(data):
                try: 
                    edges.insert(0,0)
                except TypeError:
                    edges=edges.tolist()
                    edges.insert(0,0)
        if edgeLoc=="mid":
            tmp=[]
            width=edges[1]-edges[0]
            tmp.append(edges[0]-0.5*width)
            for i in range(0,len(edges)):
                tmp.append(edges[i]+0.5*width)
            edges=tmp

        # Build histogram data
        for i in range(1,len(data)+1):
            histoEdges.append(edges[i-1])
            histoData.append(data[i-1])
            histoEdges.append(edges[i])
            histoData.append(data[i-1])
            histoMidX.append((edges[i-1]+edges[i])/2)
            histoMidY.append(data[i-1])

        histoEdges.append(edges[-1])
        histoData.append(0.0)

        # Store data in histogram object
        self.xEdges=histoEdges
        self.yValues=histoData
        self.midPtX=histoMidX
        self.midPtY=histoMidY
        self.sigma=uncert 