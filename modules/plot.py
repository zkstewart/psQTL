import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pycirclize import Circos
from matplotlib.patches import Patch

from .gff3 import GFF3

SAMPLE_AESTHETICS = [["#000000", "dotted"], ["#002D7E", "dashed"], ["#ECE45A", "dashdot"]]
COVERAGE_COLOURS = ["#004488", "#ddaa33"] # set aside to ensure contrast of colours
NUM_SAMPLE_LINES = len(SAMPLE_AESTHETICS) # for validation

def bin_values(values, start, end, binSize, binThreshold):
    '''
    Parameters:
        values -- a list of lists containing three values: [position, contigID, statistical value]
        start -- an integer value indicating the start of the region
        end -- an integer value indicating the end of the region
        binSize -- an integer value indicating the size of the bins/windows
        binThreshold -- an integer value indicating the threshold for counting a variant
                        within a bin/window
    Returns:
        histo -- a numpy array with length equal to the number of bins/windows
                 within the given region boundaries, with each value indicating
                 the number of variants within that bin/window that met the
                 binThreshold
    '''
    histo = np.array([
        0 for _ in range(math.ceil((end - start + 1) / binSize))
    ])
    for pos, _, stat in values:
        if abs(stat) >= binThreshold:
            binIndex = (pos-start) // binSize
            histo[binIndex] += 1
    return histo

def WMA(s, period):
    """
    See https://stackoverflow.com/questions/74518386/improving-weighted-moving-average-performance
    
    Parameters:
        s -- a numpy array of values to smooth
        period -- an integer value indicating the number of previous values to consider
                  during weighted moving average calculation
    Returns:
        sw -- a pandas Series of the smoothed values
    """
    w = np.arange(period)+1
    w_s = w.sum()
    
    try:
        swv = np.lib.stride_tricks.sliding_window_view(s.flatten(), window_shape=period)
    except ValueError:
        "Less data points than period size causes this error"
        return None
    sw = (swv * w).sum(axis=1) / w_s
    
    # Need to now return it as a normal series
    sw = np.concatenate((np.full(period - 1, np.nan), sw))
    try:
        sw[0:period] = sw[period] # set first n=period values to be same as first smoothed value
    except:
        "len(sw)==1 causes this error"
        return None
    return pd.Series(sw)

class Plot:
    RESULT_TYPES = ["depth", "call"]
    MEASUREMENT_TYPES = ["ed", "splsda"]
    PLOT_TYPES = ["line", "scatter", "histogram", "coverage", "genes"]
    STANDARD_DIMENSION = 5
    
    def __init__(self, resultTypes, measurementTypes, plotTypes, regions,
                 wmaSize=5, binSize=100000, binThreshold=0.4,
                 width=None, height=None):
        # Mandatory parameters
        self.resultTypes = resultTypes
        self.measurementTypes = measurementTypes
        self.plotTypes = plotTypes
        self.regions = regions
        
        # Optional parameters
        self.wmaSize = wmaSize
        self.binSize = binSize
        self.binThreshold = binThreshold
        self.width = width
        self.height = height
        
        # Figure-related parameters (not to be set by user)
        self.fig = None
        self.axs = None
        self.rowNum = None
    
    def start_plotting(self):
        raise NotImplementedError("start_plotting() must be implemented in subclasses")
    
    def savefig(self, outputFileName):
        raise NotImplementedError("savefig() must be implemented in subclasses")
    
    def set_col_labels(self, labels):
        raise NotImplementedError("set_col_labels() must be implemented in subclasses")
    
    def set_row_labels(self, labels):
        raise NotImplementedError("set_row_labels() must be implemented in subclasses")
    
    @property
    def resultTypes(self):
        return self._resultTypes
    
    @resultTypes.setter
    def resultTypes(self, valueList):
        if not isinstance(valueList, list):
            raise TypeError("resultTypes must be a list")
        if len(valueList) < 1 or len(valueList) > 2:
            raise ValueError("resultTypes must be a list of length 1 or 2")
        for value in valueList:
            if value not in Plot.RESULT_TYPES:
                raise ValueError(f"resultTypes must be one of {Plot.RESULT_TYPES}")
        
        self._resultTypes = valueList
    
    @property
    def measurementTypes(self):
        return self._measurementTypes
    
    @measurementTypes.setter
    def measurementTypes(self, valueList):
        if not isinstance(valueList, list):
            raise TypeError("measurementTypes must be a list")
        if len(valueList) < 1 or len(valueList) > 2:
            raise ValueError("measurementTypes must be a list of length 1 or 2")
        for value in valueList:
            if value not in Plot.MEASUREMENT_TYPES:
                raise ValueError(f"measurementTypes must be one of {Plot.MEASUREMENT_TYPES}")
        
        self._measurementTypes = valueList
    
    @property
    def plotTypes(self):
        return self._plotTypes
    
    @plotTypes.setter
    def plotTypes(self, valueList):
        if not isinstance(valueList, list):
            raise TypeError("plotTypes must be a list")
        if len(valueList) < 1 or len(valueList) > len(Plot.PLOT_TYPES):
            raise ValueError("plotTypes must be a list of length 1 or 2")
        for value in valueList:
            if value not in Plot.PLOT_TYPES:
                raise ValueError(f"plotTypes must be one of {Plot.PLOT_TYPES}")
        
        self._plotTypes = valueList
    
    @property
    def regions(self):
        return self._regions
    
    @regions.setter
    def regions(self, value):
        "Validation should have occurred prior to setting this property"
        if not isinstance(value, list):
            raise TypeError("regions must be a list")
        if len(value) < 1:
            raise ValueError("regions must be a list of length >= 1")
        
        self._regions = value
    
    @property
    def wmaSize(self):
        return self._wmaSize
    
    @wmaSize.setter
    def wmaSize(self, value):
        if not isinstance(value, int):
            raise TypeError("wmaSize must be an integer")
        if value < 1:
            raise ValueError(f"wmaSize must be >= 1")
        
        self._wmaSize = value
    
    @property
    def binSize(self):
        return self._binSize
    
    @binSize.setter
    def binSize(self, value):
        if not isinstance(value, int):
            raise TypeError("binSize must be an integer")
        if value < 2:
            raise ValueError(f"binSize must be >= 2")
        
        self._binSize = value
    
    @property
    def binThreshold(self):
        return self._binThreshold
    
    @binThreshold.setter
    def binThreshold(self, value):
        if not isinstance(value, float) and not isinstance(value, int):
            raise TypeError("binThreshold must be a float or int")
        if value < 0:
            raise ValueError(f"binThreshold must be >= 0")
        
        self._binThreshold = value
    
    @property
    def ncol(self):
        return len(self.regions)
    
    @property
    def nrow(self):
        oneRowTypes = [ x for x in ["coverage", "genes"] if x in self.plotTypes ]
        otherPlotTypes = [ x for x in self.plotTypes if x not in oneRowTypes ]
        
        return (len(self.resultTypes) * len(self.measurementTypes) * \
            (
                len(otherPlotTypes) if not all(["scatter" in otherPlotTypes, "line" in otherPlotTypes ]) 
                else len(otherPlotTypes) - 1 # scatter and line are in the same row
            )) + len(oneRowTypes)
    
    @property
    def width(self):
        if self._width is None:
            return Plot.STANDARD_DIMENSION * self.ncol
        return self._width
    
    @width.setter
    def width(self, value):
        if value == None:
            self._width = None
            return
        else:
            if not isinstance(value, int):
                raise TypeError("width must be an integer")
            if value < 1:
                raise ValueError(f"width must be >= 1")
            
            self._width = value
    
    @property
    def height(self):
        if self._height is None:
            return Plot.STANDARD_DIMENSION * self.nrow
        return self._height
    
    @height.setter
    def height(self, value):
        if value == None:
            self._height = None
            return
        else:
            if not isinstance(value, int):
                raise TypeError("height must be an integer")
            if value < 1:
                raise ValueError(f"height must be >= 1")

            self._height = value
    
    def scatter(self, windowedNCLS, contigID, start, end):
        '''
        Returns data suited for scatter plotting of WindowedNCLS values.
        
        Parameters:
            windowedNCLS -- a WindowedNCLS object with statistical values
                            queryable by contigID and start/end positions
            contigID -- a string indicating the contig ID
            start -- an integer indicating the start position of the region
            end -- an integer indicating the end position of the region
        Returns:
            x -- a numpy array of the x values (positions)
            y -- a numpy array of the y values (statistical values)
        '''
        regionValues = windowedNCLS.find_overlap(contigID, start, end)
        x, y = [], []
        for pos, _, ed in regionValues:
            x.append(pos)
            y.append(ed)
        x = np.array(x)
        y = np.array(y)
        
        return x, y
    
    def line(self, windowedNCLS, contigID, start, end):
        '''
        Returns data suited for line plotting of WindowedNCLS values.
        
        Parameters:
            windowedNCLS -- a WindowedNCLS object with statistical values
                            queryable by contigID and start/end positions
            contigID -- a string indicating the contig ID
            start -- an integer indicating the start position of the region
            end -- an integer indicating the end position of the region
        Returns:
            x -- a numpy array of the x values (positions)
            smoothedY -- a pandas Series of the smoothed y values (statical
                         value) values OR the original y values if smoothing
                         was not possible (i.e., not enough data points)
        '''
        x, y = self.scatter(windowedNCLS, contigID, start, end)
        smoothedY = WMA(y, self.wmaSize)
        if smoothedY is None:
            print(f"WARNING: region '{contigID, start, end}' has too few data points to apply WMA smoothing")
            smoothedY = y
        return x, smoothedY
    
    def histogram(self, windowedNCLS, contigID, start, end):
        '''
        Returns data suited for histogram plots of binned ED values.
        
        Parameters:
            windowedNCLS -- a WindowedNCLS object with statistical values
                            queryable by contigID and start/end positions
            contigID -- a string indicating the contig ID
            start -- an integer indicating the start position of the region
            end -- an integer indicating the end position of the region
        Returns:
            x -- a numpy array of the x values (positions)
            y -- a numpy array of the y values (binned statistical values)
        '''
        y = bin_values(windowedNCLS.find_overlap(contigID, start, end),
                        start, end, self.binSize, self.binThreshold)
        x = np.array([ (i * self.binSize) + start for i in range(len(y)) ])
        return x, y
    
    def genes(self, gff3Obj, contigID, start, end):
        '''
        Returns data suited for histogram plots of binned ED values.
        
        Parameters:
            gff3Obj -- a GFF3 object with ncls methods to locate geneFeatures
                       within given regions
            contigID -- a string indicating the contig ID
            start -- an integer indicating the start position of the region
            end -- an integer indicating the end position of the region
        Returns:
            mrnaFeatures -- a list of mRNA features within the region
        '''
        # Get longest isoform for each gene in this region
        geneFeatures = gff3Obj.ncls_finder(start, end, "contig", contigID)
        mrnaFeatures = [
            GFF3.longest_isoform(geneFeature)
            for geneFeature in geneFeatures
            if hasattr(geneFeature, "mRNA")
        ]
        return mrnaFeatures
    
    def coverage(self, depthNCLSDict, samples, contigID, start, end):
        '''
        Returns data suited for coverage plotting of depthNCLSDict values.
        
        Parameters:
            depthNCLSDict -- a dictionary with structure like:
                         {
                             "bulk1": {
                                 "sample1": WindowedNCLS,
                                 "sample2": WindowedNCLS,
                                 ...
                            },
                             "bulk2": { ... }
                         }
            samples -- a list of sample names to plot individually
            contigID -- a string indicating the contig ID
            start -- an integer indicating the start position of the region
            end -- an integer indicating the end position of the region
        '''
        coverageData = {}
        
        # Get the data for each bulk
        for bulk in ["bulk1", "bulk2"]:
            sampleDict = depthNCLSDict[bulk]
            coverageData[bulk] = {}
            
            # Get the average and upper/lower quantiles of each bulk's coverage values
            bulkValues = []
            for sample, windowedNCLS in sampleDict.items():
                # Skip any samples being individually plotted
                if sample in samples:
                    continue
                
                # Get values within this region
                overlappingBins = list(windowedNCLS.find_overlap(contigID, start, end))
                y = [ stat for _, _, stat in overlappingBins ]
                bulkValues.append(y)
            bulkValues = np.array(bulkValues)
            
            # Skip if there are no values
            if bulkValues.size == 0:
                coverageData[bulk] = None
                continue
            
            # Get the median and Q1/Q3 quantiles
            q1, median, q3 = np.percentile(bulkValues, [25, 50, 75], axis=0)
            
            # Figure out the x values
            x = []
            for index, (windowStart, windowEnd, _) in enumerate(overlappingBins):
                x.append(windowStart + (windowEnd - windowStart)/2)
            
            # Store the data
            coverageData[bulk] = {
                "x": x,
                "q1": q1,
                "median": median,
                "q3": q3
            }
        
        # Get the data for individual samples
        for sampleIndex, sample in enumerate(samples):
            # Get values within this region
            sampleWindowedNCLS = depthNCLSDict["bulk1"][sample] \
                if sample in depthNCLSDict["bulk1"] \
                else depthNCLSDict["bulk2"][sample]
            y = [ stat for _, _, stat in sampleWindowedNCLS.find_overlap(contigID, start, end) ]
            
            # Store the data
            coverageData[sample] = y # x can be reused from the bulk data
        
        # Return the data
        return coverageData
    
    def __repr__(self):
        return f"Plot(resultTypes={self.resultTypes}, measurementTypes={self.measurementTypes}, plotTypes={self.plotTypes})"

class HorizontalPlot(Plot):
    YLIM_HEADSPACE = 0.1 # proportion of ylim to add to the top of the plot
    SPACING = 0.1 # padding for gene plots
    
    def __init__(self, resultTypes, measurementTypes, plotTypes, regions,
                 wmaSize=5, binSize=100000, binThreshold=0.4,
                 width=None, height=None):
        super().__init__(resultTypes, measurementTypes, plotTypes, regions, 
                         wmaSize, binSize, binThreshold,
                         width, height)
    
    def start_plotting(self):
        '''
        Initialises a matplotlib figure and axes for plotting. Method is not called
        immediately to allow for customisation of optional attributes especially
        the figure width and height if not set during object initialisation.
        '''
        self.fig, self.axs = plt.subplots(nrows=self.nrow, ncols=self.ncol,
                                          figsize=(self.width, self.height))
        self.axs = np.reshape(self.axs, (self.nrow, self.ncol)) # ensure shape is as expected
        self.fig.tight_layout()
        self.rowNum = -1 # to keep track of the current row number
    
    def savefig(self, outputFileName):
        self.fig.savefig(outputFileName, bbox_inches="tight")
    
    def set_col_labels(self, labels):
        '''
        Sets the column labels for the plot.
        
        Parameters:
            labels -- a list of strings indicating the labels for each column
        '''
        if self.axs is None:
            raise ValueError("Call .start_plotting() before setting column labels")
        if len(labels) != self.ncol:
            raise ValueError(f"Number of labels ({len(labels)}) does not match number of columns ({self.ncol})")
        
        for ax, label in zip(self.axs[0], labels):
            ax.set_title(label, fontweight="bold")
    
    def set_row_labels(self, labels):
        '''
        Sets the row labels for the plot.
        
        Parameters:
            labels -- a list of strings indicating the labels for each row
        '''
        if self.axs is None:
            raise ValueError("Call .start_plotting() before setting column labels")
        if len(labels) != self.nrow:
            raise ValueError(f"Number of labels ({len(labels)}) does not match number of rows ({self.nrow})")
        
        for ax, label in zip(self.axs[:,0], labels):
            ax.set_ylabel(label)
    
    def plot_linescatter(self, scatterNCLS, lineNCLS,
                         linewidth=1, dotsize=3):
        '''
        Plots the data for a line or scatter plot.
        
        Parameters:
            scatterNCLS -- a WindowedNCLS object with statistical values
                            queryable by contigID and start/end positions;
                            used for scatter plots
            lineNCLS -- a WindowedNCLS object with statistical values
                        queryable by contigID and start/end positions;
                        used for line plots
        '''
        if self.fig is None:
            self.start_plotting()
        self.rowNum += 1 # increment row number for plotting
        plotScaleBar = self.rowNum+1 == self.nrow # set up scale bar if this is the last row
        
        maxY = 0 # to set y limits at end
        for colNum, (contigID, start, end, reverse) in enumerate(self.regions):
            # Set xlim
            self.axs[self.rowNum, colNum].set_xlim(start, end)
            
            # Turn off ytick labels if not the first column
            if colNum > 0:
                self.axs[self.rowNum, colNum].set_yticklabels([])
            
            # Plot scatter (if applicable)
            if "scatter" in self.plotTypes and contigID in scatterNCLS.contigs:
                x, y = self.scatter(scatterNCLS, contigID, start, end)
                self.axs[self.rowNum, colNum].scatter(x, y, color="red", s=dotsize,
                                                      alpha=0.5, zorder=0)
                if y.size != 0:
                    maxY = max(maxY, max(y))
            
            # Plot line (if applicable)
            if "line" in self.plotTypes and contigID in lineNCLS.contigs:
                x, smoothedY = self.line(lineNCLS, contigID, start, end)
                self.axs[self.rowNum, colNum].plot(x, smoothedY, linewidth=linewidth,
                                                   zorder=1)
                if smoothedY.size != 0:
                    maxY = max(maxY, max(smoothedY))
            
            # Set up scale bar if this is the last row
            if plotScaleBar == True:
                self.scalebar(colNum, start, end)
            
            # Otherwise, turn off x labels
            else:
                self.axs[self.rowNum, colNum].set_xticklabels([])
                self.axs[self.rowNum, colNum].locator_params(axis='x', nbins=4) # use less ticks; avoid clutter
            
            # Reverse x axis if necessary
            if reverse:
                self.axs[self.rowNum, colNum].invert_xaxis()
        
        # Set y limits based on the maximum Y value
        for colNum in range(len(self.regions)):
            if maxY == 0:
                self.axs[self.rowNum, colNum].set_ylim(0, 1)
            else:
                self.axs[self.rowNum, colNum].set_ylim(0, maxY + (maxY * HorizontalPlot.YLIM_HEADSPACE))
    
    def plot_histogram(self, windowedNCLS):
        '''
        Plots the data for a histogram.
        
        Parameters:
            windowedNCLS -- a WindowedNCLS object with statistical values
                            queryable by contigID and start/end positions
        '''
        if self.fig is None:
            self.start_plotting()
        self.rowNum += 1 # increment row number for plotting
        plotScaleBar = self.rowNum+1 == self.nrow # set up scale bar if this is the last row
        
        # Get the maximum Y value across all regions
        maxY = 0
        for contigID, start, end, reverse in self.regions:
            if contigID in windowedNCLS.contigs:
                x, y = self.histogram(windowedNCLS, contigID, start, end)
                if y.size != 0:
                    maxY = max(maxY, max(y))
        
        # Plot each region
        for colNum, (contigID, start, end, reverse) in enumerate(self.regions):
            # Set limits
            self.axs[self.rowNum, colNum].set_xlim(start, end)
            if maxY == 0:
                self.axs[self.rowNum, colNum].set_ylim(0, 1)
            else:
                self.axs[self.rowNum, colNum].set_ylim(0, maxY + (maxY * HorizontalPlot.YLIM_HEADSPACE))
            
            # Turn off ytick labels if not the first column
            if colNum > 0:
                self.axs[self.rowNum, colNum].set_yticklabels([])
            
            # Plot bars
            if contigID in windowedNCLS.contigs:
                x, y = self.histogram(windowedNCLS, contigID, start, end)
                self.axs[self.rowNum, colNum].bar(x, y, align="edge", width=self.binSize)
            
            # Set up scale bar if this is the last row
            if plotScaleBar == True:
                self.scalebar(colNum, start, end)
            
            # Otherwise, turn off x labels
            else:
                self.axs[self.rowNum, colNum].set_xticklabels([])
                self.axs[self.rowNum, colNum].locator_params(axis='x', nbins=4) # use less ticks; avoid clutter
            
            # Reverse x axis if necessary
            if reverse:
                self.axs[self.rowNum, colNum].invert_xaxis()
    
    def plot_genes(self, gff3Obj):
        '''
        Plots the data for gene models.
        
        Parameters:
            gff3Obj -- a GFF3 class object from gff3.py in this repository
        '''
        if self.fig is None:
            self.start_plotting()
        self.rowNum += 1 # increment row number for plotting
        plotScaleBar = self.rowNum+1 == self.nrow # set up scale bar if this is the last row
        
        self.fig.canvas.draw() # need to draw the figure to get the renderer
        alreadyWarned = False
        
        for colNum, (contigID, start, end, reverse) in enumerate(self.regions):        
            # Get longest isoform for each gene in this region
            mrnaFeatures = self.genes(gff3Obj, contigID, start, end)
            
            # Turn off y labels
            self.axs[self.rowNum, colNum].set_yticklabels([])
            self.axs[self.rowNum, colNum].set_yticks([])
            
            # Set xlim to allow renderer to calculate text box positions correctly
            self.axs[self.rowNum, colNum].set_xlim(start, end)
            
            # Plot each gene in non-overlapping lanes
            lanes = []
            for mrnaFeature in mrnaFeatures:
                if not hasattr(mrnaFeature, "exon"):
                    if alreadyWarned == False:
                        print(f"WARNING: Gene '{mrnaFeature.ID}' has no exons; skipping...")
                        print("(This warning will only be shown once)")
                        alreadyWarned = True
                    continue
                
                # Resolve overlaps
                placed = False
                laneNum = 0
                for lane in lanes:
                    if (mrnaFeature.start - ARROW_SIZE) > lane[-1]: # ARROW_SIZE not used on first iter, gets set later
                        placed = True
                        break
                    laneNum += 1
                if not placed:
                    lane = [] # start a new lane
                    lanes.append(lane)
                
                # Plot intron line
                "Intron line, when at the bottommost layer, can just be the length of the gene"
                self.axs[self.rowNum, colNum].plot([max(mrnaFeature.start, start), min(mrnaFeature.end, end)], # truncated to region
                    [laneNum + HorizontalPlot.SPACING + (1-HorizontalPlot.SPACING)/2]*2, # y values for start and end
                    color="black", linewidth=1, zorder=0)
                
                # Plot gene directionality
                if mrnaFeature.strand in ["+", "-"]:
                    if mrnaFeature.strand == "+":
                        lastExon = max(mrnaFeature.exon, key=lambda x: x.end)
                        arrowAnnot = self.axs[self.rowNum, colNum].annotate("▷",
                            xy=(lastExon.end, laneNum + HorizontalPlot.SPACING + (1-HorizontalPlot.SPACING)/2), xycoords="data",
                            fontsize=8, color="grey", zorder=1, ha="left", va="center")
                    else:
                        lastExon = min(mrnaFeature.exon, key=lambda x: x.start)
                        arrowAnnot = self.axs[self.rowNum, colNum].annotate("◁",
                            xy=(lastExon.start, laneNum + HorizontalPlot.SPACING + (1-HorizontalPlot.SPACING)/2), xycoords="data",
                            fontsize=8, color="grey", zorder=1, ha="right", va="center")
                
                # Plot exon boxes
                exonCoords = [
                    (
                        max(exonFeature.start, start),
                        min(exonFeature.end, end)
                    )
                    for exonFeature in mrnaFeature.exon
                ]
                self.axs[self.rowNum, colNum].broken_barh(
                    [
                        (exonStart, exonEnd - exonStart)
                        for exonStart, exonEnd in exonCoords
                        if exonStart < exonEnd # this occurs if the exon exists outside of the specified region
                    ],
                    (laneNum+HorizontalPlot.SPACING, 1-HorizontalPlot.SPACING),
                    facecolors="dodgerblue", edgecolor="black",
                    zorder=2 # above gene directionality
                )
                
                # Plot CDS boxes
                if hasattr(mrnaFeature, "CDS"):
                    cdsCoords = [
                        (
                            max(cdsFeature.start, start),
                            min(cdsFeature.end, end)
                        )
                        for cdsFeature in mrnaFeature.CDS
                    ]
                    self.axs[self.rowNum, colNum].broken_barh(
                        [
                            (cdsStart, cdsEnd - cdsStart)
                            for cdsStart, cdsEnd  in cdsCoords
                            if cdsStart < cdsEnd # this occurs if the exon exists outside of the specified region
                        ],
                        (laneNum+HorizontalPlot.SPACING, 1-HorizontalPlot.SPACING),
                        facecolors="coral", edgecolor="black",
                        zorder=3 # above exon boxes and gene directionality
                    )
                
                # Get the size of the arrow
                transf = self.axs[self.rowNum, colNum].transData.inverted()
                arrowBbox = arrowAnnot.get_window_extent(renderer = self.fig.canvas.renderer)
                arrowDataCoords = arrowBbox.transformed(transf)
                ARROW_SIZE = arrowDataCoords.x1 - arrowDataCoords.x0
                
                # Plot gene name
                geneName = mrnaFeature.ID
                namePosition = mrnaFeature.end + ARROW_SIZE if (mrnaFeature.end + ARROW_SIZE) < end \
                    else end + ARROW_SIZE # prevent text from going off the plot
                
                textBox = self.axs[self.rowNum, colNum].text(
                    namePosition, 
                    laneNum + HorizontalPlot.SPACING + (1-HorizontalPlot.SPACING)/2, # y position
                    geneName, # text
                    horizontalalignment="left", verticalalignment="center", # alignment
                    fontsize=8,
                    zorder=4 # above everything else
                )
                bb = textBox.get_window_extent(renderer = self.fig.canvas.renderer)
                bb_datacoords = bb.transformed(transf)
                
                # Store the rightmost x value for this lane
                lane.append(bb_datacoords.x1)
                
                # Reposition gene name and arrow if necessary
                if reverse:
                    # Reposition the name
                    textBox.remove()
                    
                    namePosition = mrnaFeature.start - ARROW_SIZE if (mrnaFeature.start - ARROW_SIZE) > start \
                        else start - ARROW_SIZE # prevent text from going off the plot
                    textBox = self.axs[self.rowNum, colNum].text(
                        namePosition,
                        laneNum + HorizontalPlot.SPACING + (1-HorizontalPlot.SPACING)/2, # y position
                        geneName, # text
                        horizontalalignment="left", verticalalignment="center", # alignment
                        fontsize=8,
                        zorder=4 # above everything else
                    )
                    
                    # Reposition the arrow
                    arrowAnnot.remove()
                    if mrnaFeature.strand == "+":
                        lastExon = max(mrnaFeature.exon, key=lambda x: x.end)
                        arrowAnnot = self.axs[self.rowNum, colNum].annotate("◁",
                            xy=(lastExon.end, laneNum + HorizontalPlot.SPACING + (1-HorizontalPlot.SPACING)/2), xycoords="data",
                            fontsize=8, color="grey", zorder=1, ha="right", va="center")
                    else:
                        lastExon = min(mrnaFeature.exon, key=lambda x: x.start)
                        arrowAnnot = self.axs[self.rowNum, colNum].annotate("▷",
                            xy=(lastExon.start, laneNum + HorizontalPlot.SPACING + (1-HorizontalPlot.SPACING)/2), xycoords="data",
                            fontsize=8, color="grey", zorder=1, ha="left", va="center")
            
            # Set up scale bar if this is the last row
            if plotScaleBar == True:
                self.scalebar(colNum, start, end)
            # Otherwise, turn off x labels
            else:
                self.axs[self.rowNum, colNum].set_xticklabels([])
                self.axs[self.rowNum, colNum].locator_params(axis='x', nbins=4) # use less ticks; avoid clutter
            
            # Reverse x axis if necessary
            if reverse:
                self.axs[self.rowNum, colNum].invert_xaxis()
        
        # Set ylim to the maximum number of lanes
        if len(lanes) == 0:
            self.axs[self.rowNum, colNum].set_ylim(0, 1)
        else:
            self.axs[self.rowNum, colNum].set_ylim(0, len(lanes)+HorizontalPlot.SPACING)
    
    def plot_coverage(self, depthNCLSDict, samples,
                      linewidth=1):
        '''
        Parameters:
            depthNCLSDict -- a dictionary with structure like:
                            {
                                "bulk1": {
                                    "sample1": EDNCLS,
                                    "sample2": EDNCLS,
                                    ...
                                },
                                "bulk2": { ... }
                            }
            samples -- a list of strings indicating the sample names to plot individual
                       lines for
            linewidth -- OPTIONAL; an integer value indicating the width of the line plot (default=1)
        '''
        if self.fig is None:
            self.start_plotting()
        self.rowNum += 1 # increment row number for plotting
        plotScaleBar = self.rowNum+1 == self.nrow # set up scale bar if this is the last row
        
        # Plot each region
        minY = 0 # to set y limits at end
        maxY = 0
        for colNum, (contigID, start, end, reverse) in enumerate(self.regions):
            coverageData = self.coverage(depthNCLSDict, samples, contigID, start, end) # keys: bulk1, bulk2, [*samples]
            
            # Plot each bulk
            for bulk in ["bulk1", "bulk2"]:
                bulkData = coverageData[bulk] # keys: x, q1, median, q3
                
                # Skip if there are no values
                if bulkData == None:
                    print(f"WARNING: No coverage values found for region '{contigID, start, end}' for '{bulk}'")
                    continue
                
                # Extend tails for better visualisation
                x = np.concatenate(([start], bulkData["x"], [end]))
                median = np.concatenate((
                    [bulkData["median"][0]], bulkData["median"], [bulkData["median"][-1]]
                ))
                q1 = np.concatenate((
                    [bulkData["q1"][0]], bulkData["q1"], [bulkData["q1"][-1]]
                ))
                q3 = np.concatenate((
                    [bulkData["q3"][0]], bulkData["q3"], [bulkData["q3"][-1]]
                ))
                
                # Plot median and Q1/Q3 lines
                self.axs[self.rowNum, colNum].plot(
                    x, median, linewidth=linewidth,
                    color=COVERAGE_COLOURS[0] if bulk == "bulk1" else COVERAGE_COLOURS[1])
                
                self.axs[self.rowNum, colNum].fill_between(
                    x, q1, q3, alpha = 0.5,
                    color=COVERAGE_COLOURS[0] if bulk == "bulk1" else COVERAGE_COLOURS[1],
                    label="_nolegend_")
                
                # Get the maximum Y value for this region
                maxY = max(maxY, np.percentile(median, 90)) # 90th percentile to trim outliers
            
            # Plot individual samples
            for sampleIndex, sample in enumerate(samples):
                "x can be reused from the bulk plot"
                y = coverageData[sample]
                
                # Extend tails for better visualisation
                y = np.concatenate(([y[0]], y, [y[-1]]))
                
                # Get the line colour and type
                lineColour, lineType = SAMPLE_AESTHETICS[sampleIndex]
                
                # Plot the line
                self.axs[self.rowNum, colNum].plot(x, y, color=lineColour,
                                              linestyle=lineType,
                                              linewidth=linewidth)
            
            # Set xlim
            self.axs[self.rowNum, colNum].set_xlim(start, end)
            
            # Turn off ytick labels if not the first column
            if colNum > 0:
                self.axs[self.rowNum, colNum].set_yticklabels([])
            
            # Set up scale bar if this is the last row
            if plotScaleBar == True:
                self.scalebar(colNum, start, end)
            # Otherwise, turn off x labels
            else:
                self.axs[self.rowNum, colNum].set_xticklabels([])
                self.axs[self.rowNum, colNum].locator_params(axis='x', nbins=4) # use less ticks; avoid clutter
            
            # Reverse x axis if necessary
            if reverse:
                self.axs[self.rowNum, colNum].invert_xaxis()
        
        # Set y limits
        for colNum in range(len(self.regions)):
            if np.floor(minY) == np.ceil(maxY):
                self.axs[self.rowNum, colNum].set_ylim(0, 1)
            else:
                self.axs[self.rowNum, colNum].set_ylim(np.floor(minY), np.ceil(maxY))
        
        # Set legend
        legendLabels = ["bulk1", "bulk2"] + samples # samples can be []
        self.axs[self.rowNum, colNum].legend(legendLabels, # colNum is the last column
                                        loc="center left",
                                        bbox_to_anchor=(1, 0.5),
                                        ncol=1)
    
    def scalebar(self, colNum, start, end):
        '''
        Start and end values are unused now but may be used in the future
        again.
        
        Parameters:
            colNum -- an integer value indicating the column index to plot to
            start -- [UNUSED]; an integer value indicating the start of the region
            end -- [UNUSED]; an integer value indicating the end of the region
        '''
        self.axs[self.rowNum, colNum].locator_params(axis='x', nbins=4) # use less ticks; avoid clutter
        self.axs[self.rowNum, colNum].set_xlabel("Chromosome position", weight="bold") 
    
    def __repr__(self):
        return f"HorizontalPlot(resultTypes={self.resultTypes}, measurementTypes={self.measurementTypes}, plotTypes={self.plotTypes})"

class CircosPlot(Plot):
    START_POSITION = 95
    TRACK_GAP = 2
    OUTER_HEIGHT = 0.3
    CENTRE_SPACE = 30
    COLOURS = ["#cc6677", "#332288", "#ddcc77", # Paul Tol muted colour palette
               "#117733", "#88ccee", "#882255",
               "#44aa99", "#999933", "#aa4499",
               "#dddddd"] # 10 colours as 10 rows are the max that psQTL_post can produce
    INTERVALS = {
        1: [1, "bp"], # 1 bp
        10: [1, "bp"], # 10 bp
        100: [1, "bp"], # 100 bp
        1000: [1000, "Kb"], # 1 Kb
        10000: [1000, "Kb"], # 10 Kb
        100000: [1000, "Kb"], # 100 Kb
        1000000: [1000000, "Mb"], # 1 Mb
        10000000: [1000000, "Mb"], # 10 Mb
        100000000: [1000000, "Mb"], # 100 Mb
        1000000000: [1000000000, "Gb"], # 1 Gb
        10000000000: [1000000000, "Gb"], # 10 Gb
        100000000000: [1000000000, "Gb"], # 100 Gb
        1000000000000: [1000000000000, "Tb"], # 1 Tb
        10000000000000: [1000000000000, "Tb"], # 10 Tb
        100000000000000: [1000000000000, "Tb"], # 100 Tb
        1000000000000000: [1000000000000000, "Pb"] # that's got to future-proof it for a while
    }
    NUM_MAJOR_TICKS = 5 # number of major ticks to aim for on the scale bar
    
    def __init__(self, resultTypes, measurementTypes, plotTypes, regions,
                 wmaSize=5, binSize=100000, binThreshold=0.4,
                 width=None, height=None):
        super().__init__(resultTypes, measurementTypes, plotTypes, regions, 
                         wmaSize, binSize, binThreshold,
                         width, height)
        
        # Figure-related parameters (not to be set by user)
        self.axs = None
        self.rowNum = None
    
    def start_plotting(self):
        seqid2size = {
            f"{contigID}:{start}-{end}" if not reverse else f"{contigID}:{end}-{start}": (start, end)
            for contigID, start, end, reverse in self.regions
        }
        self.circos = Circos(seqid2size, space = 0 if len(seqid2size) == 1 else 2)
        self.handles = [] # to store legend handles (which are row labels)
        
        self.trackHeight = (100 - (
            (CircosPlot.TRACK_GAP * self.nrow) +
            (100 - CircosPlot.START_POSITION) +
            CircosPlot.CENTRE_SPACE)) / self.nrow # height of each track
        
        # Establish axes for each region/sector
        "Enables us to use a similar interface to the HorizontalPlot class for plotting to specific [row,col] indices"
        self.axs = []
        for colNum, sector in enumerate(self.circos.sectors):
            currentPosition = CircosPlot.START_POSITION
            
            # Set up outer track
            outer_track = sector.add_track((currentPosition-CircosPlot.OUTER_HEIGHT, currentPosition))
            outer_track.axis(fc="black")
            
            # Determine scale bar intervals
            major_interval = 10 ** (len(str(sector.size)) - 1)
            divisor, unit = CircosPlot.INTERVALS[major_interval]
            minor_interval = int(major_interval / 10)
            major_tick_interval = math.ceil((sector.size / CircosPlot.NUM_MAJOR_TICKS) / major_interval) * major_interval
            
            # Set up the scale bar
            if sector.size > minor_interval:
                outer_track.xticks_by_interval(major_tick_interval, label_formatter=lambda v: f"{v / divisor:.0f} {unit}")
                outer_track.xticks_by_interval(minor_interval, tick_length=1, show_label=False)
            currentPosition -= (CircosPlot.OUTER_HEIGHT + CircosPlot.TRACK_GAP)
            
            # Set all inner row/tracks
            column = []
            for rowNum in range(self.nrow):
                track = sector.add_track((currentPosition-self.trackHeight, currentPosition))
                track.axis()
                column.append(track)
                currentPosition -= (self.trackHeight + CircosPlot.TRACK_GAP)
            
            # Store the column of tracks
            self.axs.append(column)
        
        self.axs = np.column_stack(self.axs)
        self.rowNum = -1 # to keep track of the current row/track number
    
    def savefig(self, outputFileName):
        fig = self.circos.plotfig()
        _ = self.circos.ax.legend(
            handles=self.handles,
            bbox_to_anchor=(0.5, 0.5),
            loc="center",
            ncols=2, ## TBD: check how this looks
        )
        fig.savefig(outputFileName, dpi=300)
    
    def set_col_labels(self, labels):
        '''
        Sets the column labels for the plot.
        
        Parameters:
            labels -- a list of strings indicating the labels for each column
        '''
        if self.axs is None:
            raise ValueError("Call .start_plotting() before setting column labels")
        if len(labels) != self.ncol:
            raise ValueError(f"Number of labels ({len(labels)}) does not match number of columns ({self.ncol})")
        
        for label, sector in zip(labels, self.circos.sectors): ## TBD: make sure ordering is stable
            sector.text(label, size=10)
    
    def set_row_labels(self, labels):
        '''
        Sets the row labels for the plot.
        
        Parameters:
            labels -- a list of strings indicating the labels for each row
        '''
        if self.axs is None:
            raise ValueError("Call .start_plotting() before setting column labels")
        if len(labels) != self.nrow:
            raise ValueError(f"Number of labels ({len(labels)}) does not match number of rows ({self.nrow})")
        
        self.handles = [
            Patch(color=CircosPlot.COLOURS[i], label=x)
            for i, x in enumerate(labels)
        ]
    
    def plot_linescatter(self, scatterNCLS, lineNCLS,
                         linewidth=1, dotsize=3):
        '''
        Plots the data for a line or scatter plot.
        
        Parameters:
            scatterNCLS -- a WindowedNCLS object with statistical values
                            queryable by contigID and start/end positions;
                            used for scatter plots
            lineNCLS -- a WindowedNCLS object with statistical values
                        queryable by contigID and start/end positions;
                        used for line plots
        '''
        if self.axs is None:
            self.start_plotting()
        self.rowNum += 1 # increment row number for plotting
        
        # Derive y limits from the maximum Y value across all regions
        maxY = 0 # to set y limits at end
        for colNum, (contigID, start, end, reverse) in enumerate(self.regions):
            if "scatter" in self.plotTypes and contigID in scatterNCLS.contigs:
                x, y = self.scatter(scatterNCLS, contigID, start, end)
                if y.size != 0:
                    maxY = max(maxY, max(y))
            if "line" in self.plotTypes and contigID in lineNCLS.contigs:
                x, smoothedY = self.line(lineNCLS, contigID, start, end)
                if smoothedY.size != 0:
                    maxY = max(maxY, max(smoothedY))
        
        # Plot each region
        for colNum, (contigID, start, end, reverse) in enumerate(self.regions):
            # Plot scatter (if applicable)
            if "scatter" in self.plotTypes and contigID in scatterNCLS.contigs:
                x, y = self.scatter(scatterNCLS, contigID, start, end)
                #x = CircosPlot.adjustX(x, start, reverse) ## TBD: implement this
                self.axs[self.rowNum, colNum].scatter(np.clip(x, start, end), y, vmax=maxY,
                                                      color=CircosPlot.COLOURS[self.rowNum],
                                                      s=dotsize, alpha=0.5,
                                                      zorder=0)
            
            # Plot line (if applicable)
            if "line" in self.plotTypes and contigID in lineNCLS.contigs:
                x, smoothedY = self.line(lineNCLS, contigID, start, end)
                smoothedY = smoothedY.to_numpy()
                #x = CircosPlot.adjustX(x, start, reverse) ## TBD: implement this
                self.axs[self.rowNum, colNum].line(np.clip(x, start, end), smoothedY, vmax=maxY,
                                                   color=CircosPlot.COLOURS[self.rowNum],
                                                   linewidth=linewidth,
                                                   zorder=1)
    
    def plot_histogram(self, windowedNCLS):
        '''
        Plots the data for a histogram.
        
        Parameters:
            windowedNCLS -- a WindowedNCLS object with statistical values
                            queryable by contigID and start/end positions
        '''
        if self.axs is None:
            self.start_plotting()
        self.rowNum += 1 # increment row number for plotting
        
        # Get the maximum Y value across all regions
        maxY = 0
        for contigID, start, end, reverse in self.regions:
            if contigID in windowedNCLS.contigs:
                x, y = self.histogram(windowedNCLS, contigID, start, end)
                if y.size != 0:
                    maxY = max(maxY, max(y))
        
        # Plot each region
        for colNum, (contigID, start, end, reverse) in enumerate(self.regions):
            if contigID in windowedNCLS.contigs:
                x, y = self.histogram(windowedNCLS, contigID, start, end)
                self.axs[self.rowNum, colNum].bar(np.clip(x, start, end), y, vmax=maxY,
                                                  width=self.binSize,
                                                  color=CircosPlot.COLOURS[self.rowNum],
                                                  align="edge")
    
    def plot_genes(self, gff3Obj):
        '''
        Plots the data for gene models.
        
        Parameters:
            gff3Obj -- a Gff().get_seqid2features(feature_tyoe=None) instance
                       from pyCirclize
        '''
        if self.axs is None:
            self.start_plotting()
        self.rowNum += 1 # increment row number for plotting
        
        for colNum, (contigID, start, end, reverse) in enumerate(self.regions):
            for feature in gff3Obj[contigID]:
                if feature.type == "gene":
                    featureStart = feature.location.start.real
                    featureEnd = feature.location.end.real
                    # Check if the feature is within the region
                    if featureEnd > start and featureStart < end:
                        self.axs[self.rowNum, colNum].genomic_features(
                            [feature], plotstyle="arrow",
                            fc=CircosPlot.COLOURS[self.rowNum])
    
    def plot_coverage(self, depthNCLSDict, samples,
                      linewidth=1):
        '''
        Parameters:
            depthNCLSDict -- a dictionary with structure like:
                            {
                                "bulk1": {
                                    "sample1": EDNCLS,
                                    "sample2": EDNCLS,
                                    ...
                                },
                                "bulk2": { ... }
                            }
            samples -- a list of strings indicating the sample names to plot individual
                       lines for
            linewidth -- OPTIONAL; an integer value indicating the width of the line plot (default=1)
        '''
        if self.axs is None:
            self.start_plotting()
        self.rowNum += 1 # increment row number for plotting
        
        # Get the maximum Y value across all regions
        maxY = 0
        for contigID, start, end, reverse in self.regions:
            coverageData = self.coverage(depthNCLSDict, samples, contigID, start, end)
            for bulk in ["bulk1", "bulk2"]:
                bulkData = coverageData[bulk]
                if bulkData != None:
                    maxY = max(maxY, np.percentile(bulkData["median"], 90)) # 90th percentile to trim outliers
        if maxY == 0:
            maxY = 1 # this can occur if all values are 0
        
        # Plot each region
        for colNum, (contigID, start, end, reverse) in enumerate(self.regions):
            coverageData = self.coverage(depthNCLSDict, samples, contigID, start, end) # keys: bulk1, bulk2, [*samples]
            
            # Plot each bulk
            extendStart = False
            extendEnd = False
            for bulk in ["bulk1", "bulk2"]:
                bulkData = coverageData[bulk] # keys: x, q1, median, q3
                
                # Skip if there are no values
                if bulkData == None:
                    print(f"WARNING: No coverage values found for region '{contigID, start, end}' for '{bulk}'")
                    continue
                
                # Extend tails for better visualisation
                if bulkData["x"][0] < start:
                    bulkData["x"][0] = start
                else:
                    bulkData["x"] = np.concatenate(([start], bulkData["x"]))
                    bulkData["median"] = np.concatenate(([bulkData["median"][0]], bulkData["median"]))
                    bulkData["q1"] = np.concatenate(([bulkData["q1"][0]], bulkData["q1"]))
                    bulkData["q3"] = np.concatenate(([bulkData["q3"][0]], bulkData["q3"]))
                    extendStart = True
                
                if bulkData["x"][-1] > end:
                    bulkData["x"][-1] = end
                else:
                    bulkData["x"] = np.concatenate((bulkData["x"], [end]))
                    bulkData["median"] = np.concatenate((bulkData["median"], [bulkData["median"][-1]]))
                    bulkData["q1"] = np.concatenate((bulkData["q1"], [bulkData["q1"][-1]]))
                    bulkData["q3"] = np.concatenate((bulkData["q3"], [bulkData["q3"][-1]]))
                    extendEnd = True
                
                # Cap values at maxY
                bulkData["median"] = np.clip(bulkData["median"], 0, maxY)
                bulkData["q1"] = np.clip(bulkData["q1"], 0, maxY)
                bulkData["q3"] = np.clip(bulkData["q3"], 0, maxY)
                
                # Plot median and Q1/Q3 lines
                self.axs[self.rowNum, colNum].fill_between(
                    bulkData["x"], bulkData["q1"], bulkData["q3"], vmax=maxY,
                    alpha = 0.5, color=COVERAGE_COLOURS[0] if bulk == "bulk1" else COVERAGE_COLOURS[1],
                    label="_nolegend_", zorder=0)
                
                self.axs[self.rowNum, colNum].line(
                    bulkData["x"], bulkData["median"], vmax=maxY,
                    color=COVERAGE_COLOURS[0] if bulk == "bulk1" else COVERAGE_COLOURS[1],
                    linewidth=linewidth, zorder=1)
            
            # Plot individual samples
            for sampleIndex, sample in enumerate(samples):
                "x can be reused from the bulk plot"
                y = coverageData[sample]
                
                # Extend tails for better visualisation
                if extendStart:
                    y = np.concatenate(([y[0]], y))
                if extendEnd:
                    y = np.concatenate((y, [y[-1]]))
                
                # Cap values at maxY
                y = np.clip(y, 0, maxY)
                
                # Get the line colour and type
                lineColour, lineType = SAMPLE_AESTHETICS[sampleIndex]
                
                # Plot the line
                self.axs[self.rowNum, colNum].line(
                    bulkData["x"], y, vmax=maxY, 
                    color=lineColour, linestyle=lineType, linewidth=linewidth)
    
    def __repr__(self):
        return f"CircosPlot(resultTypes={self.resultTypes}, measurementTypes={self.measurementTypes}, plotTypes={self.plotTypes})"
