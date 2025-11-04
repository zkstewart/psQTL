import os, sys, math, re
import numpy as np
import matplotlib.pyplot as plt

from Bio.SeqFeature import SeqFeature, SimpleLocation, ExactPosition
from pycirclize import Circos
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from gff3 import GFF3Graph
from reporting import is_overlapping

SAMPLE_AESTHETICS = [["#000000", "dotted"], ["#002D7E", "dashed"], ["#ECE45A", "dashdot"]]
LINESCATTER_COLOURS = [["#2166ac", "#b2182b"], # blue to red, for ED measurements
                       ["#762a83", "#e7b745"]]  # purple to yellow, for BA measurements
INTEGRATED_AESTHETICS = ["#1b7837", "D", "Integrated SNP/CNV"] # green with diamond marker for integrated sPLSDA
SPLSDA_DOTSIZE = 15 # size of the dots for sPLSDA selected SNPs/CNVs
COVERAGE_COLOURS = ["#004488", "#ddaa33"] # set aside to ensure contrast of colours
GENE_COLOURS = ["coral", "dodgerblue"]
HIGHLIGHT_COLOUR = "#f9e659" # orange-yellow for highlights
NUM_SAMPLE_LINES = len(SAMPLE_AESTHETICS) # for validation

SMOOTHING_BUFFER = 100000 # buffer to apply to the start and end positions for smoothing

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
    return sw

def minmax_norm(x, minValue, maxValue):
    '''
    Normalize an integer x to the range [0, 1] based on the provided minValue and maxValue.
    
    Parameters:
        x -- an integer to normalise
        minValue -- the minimum value
        maxValue -- the maximum value
    Returns:
        x_norm -- a numpy array of normalized values
    '''
    if minValue == maxValue:
        return 0
    return (x - minValue) / (maxValue - minValue)

class Plot:
    RESULT_TYPES = ["depth", "call"]
    MEASUREMENT_TYPES = ["ed", "splsda"]
    PLOT_TYPES = ["line", "scatter", "histogram", "coverage", "genes"]
    STANDARD_DIMENSION = 5
    def __init__(self, regions, highlights=None,
                 callED=None, depthED=None,
                 callSPLSDA=None, depthSPLSDA=None, integratedSPLSDA=None,
                 coverageNCLSDict=None, coverageSamples=None,
                 annotationGFF3=None,
                 power=1, wmaSize=5, width=None, height=None):
        '''
        Parameters:
            regions -- a list of dictionaries, each dict structured like:
                       {
                            "contig": contigID, # string
                            "start": start, # int
                            "end": end, # int
                            "reverse": reverse] # bool
                        }
            highlights -- values equivalently formatted to 'regions' but to denote
                          regions that should be highlighted in the plot with
                          a coloured opaque background
            callED -- a WindowedNCLS object with Euclidean Distance (ED) values
                      obtained from the call method
            depthED -- a WindowedNCLS object with ED values obtained from the depth method
            callSPLSDA -- a list or tuple of two WindowedNCLS objects indicating 1)
                          selected SNPs and 2) the Balanced Accuracy (BA) for sPLSDA
                          in windows.
            depthSPLSDA -- a list or tuple of two WindowedNCLS objects indicating 1)
                           selected SNPs and 2) the BA for sPLSDA in windows.
            integratedSPLSDA -- a WindowedNCLS object with integrated sPLSDA values
                                indicating the selected SNPs and CNVs when assessing
                                selected features from 'call' and 'depth simultaneously.
            coverageNCLSDict -- a dictionary with structure like:
                                {
                                    "group1": {
                                        "sample1": WindowedNCLS,
                                        "sample2": WindowedNCLS,
                                        ...
                                    },
                                    "group2": { ... }
                                }
            coverageSamples -- a list of sample names available in coverageNCLSDict
                               that should be plotted individually OR None to not
                               plot individual samples
            annotationGFF3 -- a GFF3 object with gene annotations; type may vary
            power -- an integer indicating the power that ED values were raised to;
                     provide this for aesthetic reasons, does not affect the data
            wmaSize -- an integer indicating the size of the weighted moving average
                       (WMA) window to apply to line plots; must be >= 1
            width -- an integer indicating the width of the plot in inches; if None,
                     defaults to Plot.STANDARD_DIMENSION * number of regions
            height -- an integer indicating the height of the plot in inches; if None,
                      defaults to Plot.STANDARD_DIMENSION * number of rows
        '''
        # Mandatory parameters
        self.regions = regions
        
        # Plot data
        self.callED = callED
        self.depthED = depthED
        self.callSPLSDA = callSPLSDA
        self.depthSPLSDA = depthSPLSDA
        self.integratedSPLSDA = integratedSPLSDA
        self.coverageNCLSDict = coverageNCLSDict
        self.annotationGFF3 = annotationGFF3
        
        # Aesthetic parameters
        self.highlights = highlights
        self.power = power
        self.wmaSize = wmaSize
        self.width = width
        self.height = height
        self.coverageSamples = coverageSamples
        
        # Defaults not set during object init
        self.showGeneNames = True
        
        # Figure-related parameters (not to be set by user)
        self.fig = None
        self.axs = None
        self.rowNum = None
    
    def plot(self):
        raise NotImplementedError("plot() must be implemented in subclasses")
    
    @property
    def callED(self):
        return self._callED
    
    @callED.setter
    def callED(self, value):
        if value is None:
            self._callED = None
            return
        if (not hasattr(value, "isWindowedNCLS") or not value.isWindowedNCLS) and (not hasattr(value, "isRangeNCLS") or not value.isRangeNCLS):
            raise TypeError("callED must be a WindowedNCLS or RangeNCLS object")
        self._callED = value
    
    @property
    def depthED(self):
        return self._depthED
    
    @depthED.setter
    def depthED(self, value):
        if value is None:
            self._depthED = None
            return
        if (not hasattr(value, "isWindowedNCLS") or not value.isWindowedNCLS) and (not hasattr(value, "isRangeNCLS") or not value.isRangeNCLS):
            raise TypeError("depthED must be a WindowedNCLS or RangeNCLS object")
        self._depthED = value
    
    @property
    def callSPLSDA(self):
        return self._callSPLSDA
    
    @callSPLSDA.setter
    def callSPLSDA(self, value):
        if value is None:
            self._callSPLSDA = None
            return
        if not isinstance(value, list) and not isinstance(value, tuple):
            raise TypeError("callSPLSDA must be a list or tuple of WindowedNCLS objects")
        if not len(value) == 2:
            raise ValueError("callSPLSDA must be a list or tuple of length 2")
        for i, val in enumerate(value):
            if (not hasattr(val, "isWindowedNCLS") or not val.isWindowedNCLS) and (not hasattr(val, "isRangeNCLS") or not val.isRangeNCLS):
                raise TypeError(f"callSPLSDA[{i}] must be a WindowedNCLS or RangeNCLS object")
        self._callSPLSDA = value
    
    @property
    def depthSPLSDA(self):
        return self._depthSPLSDA
    
    @depthSPLSDA.setter
    def depthSPLSDA(self, value):
        if value is None:
            self._depthSPLSDA = None
            return
        if not isinstance(value, list) and not isinstance(value, tuple):
            raise TypeError("depthSPLSDA must be a list or tuple of WindowedNCLS objects")
        if not len(value) == 2:
            raise ValueError("depthSPLSDA must be a list or tuple of length 2")
        for i, val in enumerate(value):
            if (not hasattr(val, "isWindowedNCLS") or not val.isWindowedNCLS) and (not hasattr(val, "isRangeNCLS") or not val.isRangeNCLS):
                raise TypeError(f"depthSPLSDA[{i}] must be a WindowedNCLS or RangeNCLS object")
        self._depthSPLSDA = value
    
    @property
    def integratedSPLSDA(self):
        return self._integratedSPLSDA
    
    @integratedSPLSDA.setter
    def integratedSPLSDA(self, value):
        if value is None:
            self._integratedSPLSDA = None
            return
        if not isinstance(value, list) and not isinstance(value, tuple):
            raise TypeError("integratedSPLSDA must be a list or tuple of WindowedNCLS objects")
        if not len(value) == 2:
            raise ValueError("integratedSPLSDA must be a list or tuple of length 2")
        for i, val in enumerate(value):
            if (not hasattr(val, "isWindowedNCLS") or not val.isWindowedNCLS) and (not hasattr(val, "isRangeNCLS") or not val.isRangeNCLS):
                raise TypeError(f"integratedSPLSDA[{i}] must be a WindowedNCLS or RangeNCLS object")
        self._integratedSPLSDA = value
    
    @property
    def coverageNCLSDict(self):
        return self._coverageNCLSDict
    
    @coverageNCLSDict.setter
    def coverageNCLSDict(self, value):
        if value is None:
            self._coverageNCLSDict = None
            return
        if not isinstance(value, dict):
            raise TypeError("coverageNCLSDict must be a dictionary")
        if not set(value.keys()) == {"group1", "group2"}:
            raise ValueError("coverageNCLSDict must have keys 'group1' and 'group2'")
        for group, sampleDict in value.items():
            for sampleID, sampleValue in sampleDict.items():
                if (not hasattr(sampleValue, "isWindowedNCLS") or not sampleValue.isWindowedNCLS) and (not hasattr(sampleValue, "isRangeNCLS") or not sampleValue.isRangeNCLS):
                    raise TypeError("coverageNCLSDict must index WindowedNCLS or RangeNCLS objects")
        self._coverageNCLSDict = value
    
    @property
    def coverageSamples(self):
        return self._coverageSamples
    
    @coverageSamples.setter
    def coverageSamples(self, value):
        "Internally convert None to an empty list"
        if value is None:
            self._coverageSamples = []
            return
        if not isinstance(value, list):
            raise TypeError("coverageSamples must be a list")
        if len(value) > NUM_SAMPLE_LINES:
            raise ValueError(f"coverageSamples must have at most {NUM_SAMPLE_LINES} samples")
        for sample in value:
            found = False
            for key, v in self.coverageNCLSDict.items():
                if sample in v:
                    found = True
                    break
            if not found:
                raise ValueError(f"coverageSamples contains sample '{sample}' not found in coverageNCLSDict")
        
        self._coverageSamples = value
    
    @property
    def annotationGFF3(self):
        return self._annotationGFF3
    
    @annotationGFF3.setter
    def annotationGFF3(self, value):
        if value is None:
            self._annotationGFF3 = None
            return
        if not hasattr(value, "isGFF3Graph") or not value.isGFF3Graph:
            raise TypeError("annotationGFF3 must be a GFF3Graph object")
        self._annotationGFF3 = value
    
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
    def highlights(self):
        return self._highlights
    
    @highlights.setter
    def highlights(self, value):
        "Validation should have occurred prior to setting this property"
        if value is None:
            value = []
        if not isinstance(value, list):
            raise TypeError("highlights must be a list")
        
        self._highlights = value
    
    @property
    def power(self):
        return self._power
    
    @power.setter
    def power(self, value):
        if not isinstance(value, int):
            raise TypeError("power must be an integer")
        if value < 1:
            raise ValueError(f"power must be >= 1")
        
        self._power = value
    
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
    def showGeneNames(self):
        return self._showGeneNames
    
    @showGeneNames.setter
    def showGeneNames(self, value):
        if not isinstance(value, bool):
            raise TypeError("showGeneNames must be boolean")
        
        self._showGeneNames = value
    
    @property
    def ncol(self):
        return len(self.regions)
    
    @property
    def nrow(self):
        return sum([
            1 if self.callED is not None else 0,
            1 if self.depthED is not None else 0,
            1 if self.callSPLSDA is not None else 0,
            1 if self.depthSPLSDA is not None else 0,
            1 if self.coverageNCLSDict is not None else 0,
            1 if self.annotationGFF3 is not None else 0
        ])
    
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
    
    @staticmethod
    def match_scatter_to_line(x, y, positions):
        '''
        Return the y value at a given position on the line defined by x and y.
        Assumes that x is sorted, and that any given position will be within the
        range of x.
        
        Parameters:
            x -- a numpy array of x values (positions); must be sorted!
            y -- a numpy array of y values (statistical values); must match x length
            position -- a list or numpy array of X positions to match to the line
        Returns:
            matchedY -- a numpy array y values corresponding to the Y values
                        we expect to see at the given positions
        '''
        if len(x) != len(y):
            raise ValueError("x and y must have the same length")
        if len(x) < 2:
            raise ValueError("x and y must have at least two values to interpolate")
        
        matchedY = []
        for position in positions:
            # Locate the nearest index in x to the position
            xindex = np.searchsorted(x, position)
            xvalue = x[xindex]
            
            # If the position is outside the range of x, raise an error
            if xvalue < position or xindex == len(x):
                raise ValueError(f"Position {position} is outside the range of x values")
            
            # If the position is exactly at an x value, store the corresponding y value
            if xvalue == position:
                matchedY.append(y[xindex])
                continue
            
            # If the position is between two x values, interpolate to find the y value
            xprev, xnext = x[xindex - 1], x[xindex]
            yprev, ynext = y[xindex - 1], y[xindex]
            matchedY.append(Plot.interpolate(xprev, xnext, yprev, ynext, position))
        return matchedY
    
    @staticmethod
    def interpolate(xprev, xnext, yprev, ynext, position):
        '''
        Calculate the slope between two points (xprev, yprev) and (xnext, ynext)
        and return the interpolated y value at a given position.
        
        Parameters:
            xprev -- the x value of the previous point
            xnext -- the x value of the next point
            yprev -- the y value of the previous point
            ynext -- the y value of the next point
            position -- the x value at which to interpolate the y value
        Returns:
            y -- the interpolated y value at the given position
        '''
        slope = (ynext - yprev) / (xnext - xprev)
        return yprev + (slope * (position - xprev))
    
    def scatter(self, windowedNCLS, contigID, start, end, clipped=True):
        '''
        Returns data suited for scatter plotting of WindowedNCLS values.
        
        Parameters:
            windowedNCLS -- a WindowedNCLS object with statistical values
                            queryable by contigID and start/end positions
            contigID -- a string indicating the contig ID
            start -- an integer indicating the start position of the region
            end -- an integer indicating the end position of the region
            clipped -- (OPTIONAL) a boolean indicating whether to clip the
                       start and end positions to the nearest available
                       positions within the start->end range. Relevant
                       when the NCLS has been generated with a window size
                       greater than 0. Default is True.
        Returns:
            x -- a numpy array of the x values (positions)
            y -- a numpy array of the y values (statistical values)
        '''
        regionValues = windowedNCLS.find_overlap(contigID, start, end)
        x, y = [], []
        for pos, _, ed in regionValues:
            if clipped and pos < start: # can occur if windowSize > 1; this is treated as our first value
                pos = start
            if clipped and pos > end: # this should not happen
                continue
            x.append(pos)
            y.append(ed)
        x = np.array(x)
        y = np.array(y)
        
        return x, y
    
    def line(self, windowedNCLS, contigID, start, end, applyWMA=True, buffer=0):
        '''
        Returns data suited for line plotting of WindowedNCLS values. Extends or
        truncates the line to match the start and end positions, applying
        weighted moving average (WMA) smoothing if requested.
        
        Parameters:
            windowedNCLS -- a WindowedNCLS object with statistical values
                            queryable by contigID and start/end positions
            contigID -- a string indicating the contig ID
            start -- an integer indicating the start position of the region
            end -- an integer indicating the end position of the region
            applyWMA -- (OPTIONAL) a boolean indicating whether to apply WMA
                        smoothing to the line; default is True
            buffer -- (OPTIONAL) an integer indicating a buffer length to
                      apply to the start and end positions for smoothing
                      and avoiding line truncation; default is 0
        Returns:
            x -- a numpy array of the x values (positions)
            smoothedY -- a pandas Series of the smoothed y values (statical
                         value) values OR the original y values if smoothing
                         was not possible (i.e., not enough data points)
        '''
        x, y = self.scatter(windowedNCLS, contigID, start-buffer, end+buffer, clipped=False)
        if applyWMA:
            smoothedY = WMA(y, self.wmaSize)
            if smoothedY is None:
                print(f"WARNING: region '{contigID, start, end}' has too few data points to apply WMA smoothing")
                smoothedY = y
        else:
            smoothedY = y
        
        # Handle the case where there are no data points
        if len(x) == 0 or len(smoothedY) == 0:
            print(f"WARNING: region '{contigID, start, end}' has no data points to plot as a line")
            return np.array([]), np.array([]) # implicitly handles len(x) != len(smoothedY)
        
        # Locate the start and end indices in the x array
        xstart = 0
        while (xstart+1) < len(x) and x[xstart] < start:
            xstart += 1
        
        xend = len(x) - 1
        while xend > 0 and x[xend] > end:
            xend -= 1
        
        # Prepend values
        if x[xstart] != start: # if we need to extend the line backwards
            xprepend = [start]
            if xstart > 0:
                xprev, xnext = x[xstart - 1], x[xstart]
                yprev, ynext = y[xstart - 1], y[xstart]
                yprepend = [Plot.interpolate(xprev, xnext, yprev, ynext, start)]
            else:
                yprepend = [smoothedY[xstart]] # just project the first value
        else: # if x[xstart] == start; x[xstart] can never be < start because of how np.searchsorted works
            xprepend = []
            yprepend = []
        
        # Append values
        if x[xend] < end: # if we need to extend the line forwards
            xappend = [end]
            if xend+1 < len(x):
                xprev, xnext = x[xend], x[xend + 1]
                yprev, ynext = y[xend], y[xend + 1]
                yappend = [Plot.interpolate(xprev, xnext, yprev, ynext, end)]
            else:
                yappend = [smoothedY[xend]] # just project the last value
        else: # if x[xend] == end; x[xend] can never be > end because of the while loop above
            xappend = []
            yappend = []
        
        # Identify how many points we want to obtain from the x/y arrays
        arrayLength = 0
        for xValue in x:
            if start <= xValue <= end:
                arrayLength += 1
        
        # Generate the new x and smoothedY arrays
        x = np.concatenate((xprepend, x[xstart:xstart+arrayLength], xappend))
        smoothedY = np.concatenate((yprepend, smoothedY[xstart:xstart+arrayLength], yappend))
        
        return x, smoothedY
    
    def histogram(self, windowedNCLS, contigID, start, end):
        '''
        NOTE: histograms have been removed from the front end code as they
        are complicated to configure correctly, hard to explain, and of
        limited use. They are still here for reference and may be
        reintroduced in the future.
        
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
            GFF3Graph.longest_isoform(geneFeature)
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
                             "group1": {
                                 "sample1": WindowedNCLS,
                                 "sample2": WindowedNCLS,
                                 ...
                            },
                             "group2": { ... }
                         }
            samples -- a list of sample names to plot individually; provide
                       an empty list to not plot individual samples
            contigID -- a string indicating the contig ID
            start -- an integer indicating the start position of the region
            end -- an integer indicating the end position of the region
        '''
        coverageData = {}
        
        # Get the data for each group
        for group in ["group1", "group2"]:
            sampleDict = depthNCLSDict[group]
            coverageData[group] = {}
            
            # Get the average and upper/lower quantiles of each group's coverage values
            groupValues = []
            for sample, windowedNCLS in sampleDict.items():
                # Skip any samples being individually plotted
                if sample in samples:
                    continue
                
                # Get values within this region
                overlappingBins = list(windowedNCLS.find_overlap(contigID, start, end))
                y = [ stat for _, _, stat in overlappingBins ]
                groupValues.append(y)
            groupValues = np.array(groupValues)
            
            # Skip if there are no values
            if groupValues.size == 0:
                coverageData[group] = None
                continue
            
            # Get the median and Q1/Q3 quantiles
            q1, median, q3 = np.percentile(groupValues, [25, 50, 75], axis=0)
            
            # Figure out the x values
            x = [ windowStart for windowStart, _, _ in overlappingBins ]
            
            # Store the data
            coverageData[group] = {
                "x": x,
                "q1": q1,
                "median": median,
                "q3": q3
            }
        
        # Get the data for individual samples
        for sampleIndex, sample in enumerate(samples):
            # Get values within this region
            sampleWindowedNCLS = depthNCLSDict["group1"][sample] \
                if sample in depthNCLSDict["group1"] \
                else depthNCLSDict["group2"][sample]
            y = [ stat for _, _, stat in sampleWindowedNCLS.find_overlap(contigID, start, end) ]
            
            # Store the data
            coverageData[sample] = y # x can be reused from the group data
        
        # Return the data
        return coverageData
    
    def __repr__(self):
        return f"Plot(<inherit me please>)"

class HorizontalPlot(Plot):
    YLIM_HEADSPACE = 0.1 # proportion of ylim to add to the top of the plot
    SPACING = 0.1 # padding for gene plots
    
    def __init__(self, regions, highlights=None,
                 callED=None, depthED=None,
                 callSPLSDA=None, depthSPLSDA=None, integratedSPLSDA=None,
                 coverageNCLSDict=None, coverageSamples=None,
                 annotationGFF3=None,
                 power=1, wmaSize=5, width=None, height=None):
        super().__init__(regions, highlights, callED, depthED, callSPLSDA, depthSPLSDA, integratedSPLSDA,
                         coverageNCLSDict, coverageSamples, annotationGFF3,
                         power, wmaSize, width, height)
    
    def plot(self, plotTypes, outputFileName):
        '''
        Initialises a matplotlib figure and axes for plotting. Method is not called
        immediately to allow for customisation of optional attributes especially
        the figure width and height if not set during object initialisation.
        
        Parameters:
            plotTypes -- a list of strings indicating the types of plots to create
            outputFileName -- a string indicating the file name to save the plot to
        '''
        # Validate plot types
        if (plotTypes == None) or (not isinstance(plotTypes, list)) or (len(plotTypes) == 0):
            raise ValueError("plotTypes must be a non-empty list of plot types")
        for plotType in plotTypes:
            if plotType not in Plot.PLOT_TYPES:
                raise ValueError(f"Invalid plot type '{plotType}'; must be one of {Plot.PLOT_TYPES}")
        if len(set(plotTypes)) != len(plotTypes):
            raise ValueError("plotTypes must not contain duplicate values")
        
        # Initialise the axes
        self.fig, self.axs = plt.subplots(nrows=self.nrow, ncols=self.ncol,
                                          figsize=(self.width, self.height))
        self.axs = np.reshape(self.axs, (self.nrow, self.ncol)) # ensure shape is as expected
        self.fig.tight_layout()
        
        # Establish column labels
        self.colLabels = [f"{regionDict['contig']}:{regionDict['end']}-{regionDict['start']}" if regionDict['reverse'] == True
                          else f"{regionDict['contig']}:{regionDict['start']}-{regionDict['end']}" if regionDict['full'] == False
                          else regionDict['contig'] # indicate only the contig if we are plotting the full contig
                          for regionDict in self.regions]
        for ax, label in zip(self.axs[0], self.colLabels):
            ax.set_title(label, fontweight="bold")
        
        # Init values for storing data during iteration
        self.rowLabels = []
        self.rowNum = -1 # to keep track of the current row number
        
        # Build the plot
        if "line" in plotTypes or "scatter" in plotTypes:
            if self.callED != None:
                self.plot_linescatter(self.callED if "scatter" in plotTypes else None,
                                      self.callED if "line" in plotTypes else None,
                                      LINESCATTER_COLOURS[0],
                                      applyWMA=True if self.wmaSize > 1 else False,
                                      lineLabel=f"WMA $ED^{self.power}$" if self.wmaSize > 1 else f"$ED^{self.power}$",
                                      scatterLabel=f"SNP $ED^{self.power}$")
                self.rowLabels.append(f"SNP $ED^{self.power}$")
            if self.callSPLSDA != None:
                self.plot_linescatter(self.callSPLSDA[0],
                                      self.callSPLSDA[1],
                                      LINESCATTER_COLOURS[1],
                                      applyWMA=True if self.wmaSize > 1 else False,
                                      lineLabel=f"WMA $BA$" if self.wmaSize > 1 else "$BA$",
                                      scatterLabel=f"Selected SNP",
                                      scatterShape="D", scatterFollowsLine=True, # Y values different to ED or BA
                                      lineZorder=0, scatterZorder=1, # make selected SNPs more visible
                                      dotsize=SPLSDA_DOTSIZE, dotAlpha=1, # make selected SNPs more visible
                                      integratedNCLS=self.integratedSPLSDA[0] if self.integratedSPLSDA != None else None)
                self.rowLabels.append(f"SNP $BA$")
            if self.depthED != None:
                self.plot_linescatter(self.depthED if "scatter" in plotTypes else None,
                                      self.depthED if "line" in plotTypes else None,
                                      LINESCATTER_COLOURS[0],
                                      lineLabel=f"WMA $ED^{self.power}$" if self.wmaSize > 1 else f"$ED^{self.power}$",
                                      scatterLabel=f"CNV $ED^{self.power}$")
                self.rowLabels.append(f"CNV $ED^{self.power}$")
            if self.depthSPLSDA != None:
                self.plot_linescatter(self.depthSPLSDA[0],
                                      self.depthSPLSDA[1],
                                      LINESCATTER_COLOURS[1],
                                      applyWMA=True if self.wmaSize > 1 else False,
                                      lineLabel=f"WMA $BA$" if self.wmaSize > 1 else "$BA$",
                                      scatterLabel=f"Selected CNV",
                                      scatterShape="D", scatterFollowsLine=True,
                                      lineZorder=0, scatterZorder=1, # make selected CNVs more visible
                                      dotsize=SPLSDA_DOTSIZE, dotAlpha=1, # make selected CNVs more visible
                                      integratedNCLS=self.integratedSPLSDA[1] if self.integratedSPLSDA != None else None)
                self.rowLabels.append(f"CNV $BA$")
        if self.coverageNCLSDict != None:
            self.plot_coverage(self.coverageNCLSDict, self.coverageSamples)
            self.rowLabels.append("Median-normalised coverage")
        if self.annotationGFF3 != None:
            self.plot_genes(self.annotationGFF3)
            self.rowLabels.append("Gene annotations")
        
        # Add background colour to highlights regions
        if self.highlights != []:
            # For each region, ...
            for colNum, regionDict in enumerate(self.regions):
                # ... find if any highlights overlap with it
                for highlightDict in self.highlights:
                    if regionDict["contig"] == highlightDict["contig"] and \
                    is_overlapping(regionDict["start"], regionDict["end"], highlightDict["start"], highlightDict["end"]):
                        # If they do, set the background colour on all rows in this column
                        for rowNum in range(self.nrow):
                            self.axs[rowNum, colNum].axvspan(highlightDict["start"], highlightDict["end"],
                                                             color=HIGHLIGHT_COLOUR, alpha=0.5, zorder=-1)
        
        # Set row labels
        for ax, label in zip(self.axs[:,0], self.rowLabels):
            ax.set_ylabel(label)
        
        # Save the figure
        self.fig.savefig(outputFileName, bbox_inches="tight")
    
    def plot_linescatter(self, scatterNCLS, lineNCLS, colours, applyWMA=True,
                         lineLabel="WMA", scatterLabel="Values", scatterShape="o",
                         lineZorder=1, scatterZorder=0, linewidth=1,
                         dotsize=3, dotAlpha=0.5,
                         integratedNCLS=None, scatterFollowsLine=False):
        '''
        Plots the data for a line or scatter plot.
        
        Parameters:
            scatterNCLS -- a WindowedNCLS object with statistical values
                            queryable by contigID and start/end positions;
                            used for scatter plots OR None if not plotting
                            scatter points
            lineNCLS -- a WindowedNCLS object with statistical values
                        queryable by contigID and start/end positions;
                        used for line plots OR None if not plotting
                        a line
            colours -- a list of two colours to use for the line (first) and scatter
                       points (second)
            applyWMA -- (OPTIONAL) a boolean indicating whether to apply WMA
                        smoothing to the line; default is True
            lineLabel -- (OPTIONAL) a string indicating the label for the line;
                          default is "WMA"
            scatterLabel -- (OPTIONAL) a string indicating the label for the scatter
                            points; default is "Values"
            scatterShape -- (OPTIONAL) a string indicating the shape of the scatter
                            points; default is "o" for standard dots
            lineZorder -- (OPTIONAL) an integer indicating the z-order of the line;
                          default is 1
            scatterZorder -- (OPTIONAL) an integer indicating the z-order of the scatter
                             points; default is 0
            linewidth -- (OPTIONAL) an integer indicating the width of the line;
                         default is 1
            dotsize -- (OPTIONAL) an integer indicating the size of the scatter
                        points; default is 3
            dotAlpha -- (OPTIONAL) a float indicating the alpha value for the scatter
                        points; default is 0.5
            integratedNCLS -- (OPTIONAL) a WindowedNCLS object with integrated sPLSDA
                              values indicating the selected SNPs and CNVs when assessing
                              selected features from 'call' and 'depth' simultaneously;
                              used for scatter plots OR None if not plotting integrated
                              values
            scatterFollowsLine -- (OPTIONAL) a boolean indicating whether the scatter y values
                                should follow the line y values; default is False.
        '''
        # Validate that one of either scatterNCLS or lineNCLS is provided
        if scatterNCLS is None and lineNCLS is None:
            raise ValueError("At least one of scatterNCLS or lineNCLS must be provided")
        
        # Validate that colours is a list of two colours
        if colours is None or len(colours) != 2:
            raise ValueError("colours must be a list of two colours")
        
        self.rowNum += 1 # increment row number for plotting
        plotScaleBar = self.rowNum+1 == self.nrow # set up scale bar if this is the last row
        
        maxY = 0 # to set y limits at end
        for colNum, regionDict in enumerate(self.regions):
            # Decompose dictionary into variables
            contigID = regionDict["contig"]
            start = regionDict["start"]
            end = regionDict["end"]
            reverse = regionDict["reverse"]
            
            # Set xlim
            self.axs[self.rowNum, colNum].set_xlim(start, end)
            
            # Turn off ytick labels if not the first column
            if colNum > 0:
                self.axs[self.rowNum, colNum].set_yticklabels([])
            
            # Plot line (if applicable)
            if lineNCLS != None and contigID in lineNCLS.contigs:
                lineX, smoothedY = self.line(lineNCLS, contigID, start, end,
                                             applyWMA=applyWMA, buffer=SMOOTHING_BUFFER)
                if smoothedY.size != 0:
                    self.axs[self.rowNum, colNum].plot(lineX, smoothedY, color=colours[0],
                                                       linewidth=linewidth,
                                                       zorder=lineZorder,
                                                       label=lineLabel)
                    maxY = max(maxY, max(smoothedY))
            
            # Plot scatter (if applicable)
            if scatterNCLS != None and contigID in scatterNCLS.contigs:
                x, y = self.scatter(scatterNCLS, contigID, start, end)
                if y.size != 0:
                    if scatterFollowsLine:
                        if lineNCLS == None:
                            raise ValueError("scatterFollowsLine is True but no lineNCLS provided; cannot interpolate scatter points")
                        y = Plot.match_scatter_to_line(lineX, smoothedY, x)
                    self.axs[self.rowNum, colNum].scatter(x, y, color=colours[1], alpha=dotAlpha,
                                                          s=dotsize, marker=scatterShape,
                                                          zorder=scatterZorder,
                                                          label=scatterLabel)
                    maxY = max(maxY, max(y))
            
            # Plot integrated sPLSDA values (if applicable)
            if integratedNCLS is not None and contigID in integratedNCLS.contigs:
                x, y = self.scatter(integratedNCLS, contigID, start, end)
                if y.size != 0:
                    if scatterFollowsLine:
                        if lineNCLS == None:
                            raise ValueError("scatterFollowsLine is True but no lineNCLS provided; cannot interpolate scatter points")
                        y = Plot.match_scatter_to_line(lineX, smoothedY, x)
                    self.axs[self.rowNum, colNum].scatter(x, y, color=INTEGRATED_AESTHETICS[0], alpha=1, # alpha always 1
                                                          s=dotsize+1, marker=INTEGRATED_AESTHETICS[1],
                                                          zorder=2, # integrated points are on top
                                                          label=INTEGRATED_AESTHETICS[2])
            
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
        
        # Set legend
        self.axs[self.rowNum, colNum].legend(
            loc="center left", bbox_to_anchor=(1, 0.5), ncol=1
        )
    
    def plot_histogram(self, windowedNCLS):
        '''
        Plots the data for a histogram.
        
        Parameters:
            windowedNCLS -- a WindowedNCLS object with statistical values
                            queryable by contigID and start/end positions
        '''
        self.rowNum += 1 # increment row number for plotting
        plotScaleBar = self.rowNum+1 == self.nrow # set up scale bar if this is the last row
        
        # Get the maximum Y value across all regions
        maxY = 0
        for regionDict in self.regions:
            if regionDict["contig"] in windowedNCLS.contigs:
                x, y = self.histogram(windowedNCLS, regionDict["contig"], regionDict["start"], regionDict["end"])
                if y.size != 0:
                    maxY = max(maxY, max(y))
        
        # Plot each region
        for colNum, regionDict in enumerate(self.regions):
            # Decompose dictionary into variables
            contigID = regionDict["contig"]
            start = regionDict["start"]
            end = regionDict["end"]
            reverse = regionDict["reverse"]
            
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
        self.rowNum += 1 # increment row number for plotting
        plotScaleBar = self.rowNum+1 == self.nrow # set up scale bar if this is the last row
        
        self.fig.canvas.draw() # need to draw the figure to get the renderer
        alreadyWarned = False
        
        for colNum, regionDict in enumerate(self.regions):
            # Decompose dictionary into variables
            contigID = regionDict["contig"]
            start = regionDict["start"]
            end = regionDict["end"]
            reverse = regionDict["reverse"]
             
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
                    facecolors=GENE_COLOURS[1], edgecolor="black", linewidth=0.5,
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
                        facecolors=GENE_COLOURS[0], edgecolor="black", linewidth=0.5,
                        zorder=3 # above exon boxes and gene directionality
                    )
                
                # Get the size of the arrow
                transf = self.axs[self.rowNum, colNum].transData.inverted()
                arrowBbox = arrowAnnot.get_window_extent(renderer = self.fig.canvas.renderer)
                arrowDataCoords = arrowBbox.transformed(transf)
                ARROW_SIZE = arrowDataCoords.x1 - arrowDataCoords.x0
                
                # Plot gene name
                if self.showGeneNames:
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
                if self.showGeneNames:
                    lane.append(bb_datacoords.x1)
                else:
                    lane.append(mrnaFeature.end + ARROW_SIZE)
                
                # Reposition gene name and arrow if necessary
                if reverse:
                    # Reposition the name
                    if self.showGeneNames:
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
                                "group1": {
                                    "sample1": EDNCLS,
                                    "sample2": EDNCLS,
                                    ...
                                },
                                "group2": { ... }
                            }
            samples -- a list of strings indicating the sample names to plot individual
                       lines for
            linewidth -- OPTIONAL; an integer value indicating the width of the line plot (default=1)
        '''
        self.rowNum += 1 # increment row number for plotting
        plotScaleBar = self.rowNum+1 == self.nrow # set up scale bar if this is the last row
        
        # Plot each region
        minY = 0 # to set y limits at end
        maxY = 0
        for colNum, regionDict in enumerate(self.regions):
            # Decompose dictionary into variables
            contigID = regionDict["contig"]
            start = regionDict["start"]
            end = regionDict["end"]
            reverse = regionDict["reverse"]
            
            # Get coverage data for this region
            coverageData = self.coverage(depthNCLSDict, samples, contigID, start, end) # keys: group1, group2, [*samples]
            
            # Plot each group
            for group in ["group1", "group2"]:
                groupData = coverageData[group] # keys: x, q1, median, q3
                
                # Skip if there are no values
                if groupData == None:
                    print(f"WARNING: No coverage values found for region '{contigID, start, end}' for '{group}'")
                    continue
                
                # Extend tails for better visualisation
                x = np.concatenate(([start], groupData["x"], [end]))
                median = np.concatenate((
                    [groupData["median"][0]], groupData["median"], [groupData["median"][-1]]
                ))
                q1 = np.concatenate((
                    [groupData["q1"][0]], groupData["q1"], [groupData["q1"][-1]]
                ))
                q3 = np.concatenate((
                    [groupData["q3"][0]], groupData["q3"], [groupData["q3"][-1]]
                ))
                
                # Plot median and Q1/Q3 lines
                self.axs[self.rowNum, colNum].plot(
                    x, median, linewidth=linewidth,
                    color=COVERAGE_COLOURS[0] if group == "group1" else COVERAGE_COLOURS[1])
                
                self.axs[self.rowNum, colNum].fill_between(
                    x, q1, q3, alpha = 0.5,
                    color=COVERAGE_COLOURS[0] if group == "group1" else COVERAGE_COLOURS[1],
                    label="_nolegend_")
                
                # Get the maximum Y value for this region
                maxY = max(maxY, np.percentile(median, 90)) # 90th percentile to trim outliers
            
            # Plot individual samples
            for sampleIndex, sample in enumerate(samples):
                "x can be reused from the group plot"
                y = coverageData[sample]
                
                # Skip if there are no values
                if y is None or len(y) == 0:
                    continue
                
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
            if np.floor(minY) == np.ceil(maxY): # minY is always 0, so this is only true if maxY is also 0
                self.axs[self.rowNum, colNum].set_ylim(0, 1)
            else:
                self.axs[self.rowNum, colNum].set_ylim(np.floor(minY), np.ceil(maxY))
        
        # Set legend
        legendLabels = ["group1", "group2"] + samples # samples can be []
        self.axs[self.rowNum, colNum].legend( # colNum is the last column
            legendLabels, loc="center left", bbox_to_anchor=(1, 0.5), ncol=1)
    
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
        return (f"HorizontalPlot(regions={self.regions}" ) + \
               (", callED=" + "True" if self.callED is not None else "False") + \
               (", depthED=" + "True" if self.depthED is not None else "False") + \
               (", callSPLSDA=" + "True" if self.callSPLSDA is not None else "False") + \
               (", depthSPLSDA=" + "True" if self.depthSPLSDA is not None else "False") + \
               (", coverage=" + "True" if self.coverageNCLSDict is not None else "False") + \
               (", genes=" + "True" if self.annotationGFF3 is not None else "False") + \
               (f", power={self.power}, wmaSize={self.wmaSize}, width={self.width}, height={self.height})")

class CircosPlot(Plot):
    # Aesthetic constants
    START_POSITION = 95
    TRACK_GAP = 3
    OUTER_HEIGHT = 0.3
    CENTRE_SPACE = 20
    AXIS_SPACE = 10
    ARC_DISTANCE = 1.5
    
    TOP=1.05
    LEFT=-0.1
    BOTTOM=-0.05
    RIGHT=1.1
    LEGEND_POSITIONS = [
        (LEFT, TOP, "upper left"), # top left
        (RIGHT, TOP, "upper right"), # top right
        (LEFT, BOTTOM, "lower left") # bottom left
    ]
    ROW_LABEL_POSITION = (RIGHT, BOTTOM, "lower right") # bottom right; reserved for row labels
    
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
    DEGREES_PER_X_TICK = 45 # the span of degrees for each x tick on the scale bar
    
    NUM_Y_TICKS = 3 # number of y ticks to aim for on each track
    STANDARD_DIMENSION = 8
    
    ASCII_LETTERS_START = ord("a") # ASCII value of 'a'
    
    def __init__(self, regions, highlights=None,
                 callED=None, depthED=None,
                 callSPLSDA=None, depthSPLSDA=None, integratedSPLSDA=None,
                 coverageNCLSDict=None, coverageSamples=None,
                 annotationGFF3=None,
                 power=1, wmaSize=5, width=None, height=None):
        super().__init__(regions, highlights, callED, depthED, callSPLSDA, depthSPLSDA, integratedSPLSDA,
                         coverageNCLSDict, coverageSamples, annotationGFF3,
                         power, wmaSize, width, height)
        
        # Figure-related parameters (not to be set by user)
        self.axs = None
        self.rowNum = None
    
    @staticmethod
    def number_formatter(value, divisor, unit):
        '''
        Parameters:
            value -- a float value to format
            divisor -- an integer value to divide the value by
            unit -- a string indicating the unit to append to the formatted value
        Returns:
            formattedNumber -- a string representing the formatted value
        '''
        dividedValue = value / divisor
        formattedValue = str(dividedValue).rstrip("0") # remove trailing zeros
        if formattedValue.endswith("."):
            formattedValue = formattedValue[:-1] # remove trailing decimal if it exists
        
        return f"{formattedValue} {unit}"
    
    @staticmethod
    def truncate(value, decimals):
        '''
        Truncates a float value to a specified number of decimal places. This prevents
        rounding up higher than maxY which triggers an error with pycirclize.
        
        Parameters:
            value -- a float value to truncate
            decimals -- an integer value indicating the number of decimal places to truncate to
        Returns:
            truncatedValue -- a float value truncated to the specified number of decimal places
        '''
        integer, decimal, scientific = re.match(r"^(\d+)(\.\d+)?(e[+-]?\d+)?$", str(value)).groups()
        if decimal is None:
            decimal = ""
        else:
            firstNonZero = re.search(r"[^0.]", decimal)
            firstNonZero = 0 if firstNonZero is None else firstNonZero.start()
            decimalLength = max(decimals, firstNonZero)
            decimal = decimal[:decimalLength+1] # +1 to include the decimal point
        if scientific is None:
            scientific = ""
        return float(f"{integer}{decimal}{scientific}")
    
    @property
    def width(self):
        if self._width is None:
            return CircosPlot.STANDARD_DIMENSION
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
            return CircosPlot.STANDARD_DIMENSION
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
    
    @property
    def axisSpace(self):
        if self._axisSpace is None:
            return CircosPlot.AXIS_SPACE
        return self._axisSpace
    
    @axisSpace.setter
    def axisSpace(self, value):
        if value == None:
            self._axisSpace = None
            return
        else:
            if not isinstance(value, int):
                raise TypeError("axisSpace must be an integer")
            if value < 0:
                raise ValueError(f"axisSpace must be >= 0")
            
            self._axisSpace = value
    
    @property
    def plottedLength(self):
        "Returns the total length of all regions"
        return sum(regionDict["end"] - regionDict["start"] for regionDict in self.regions) # end > start since reversed plotting isn't allowed
    
    def _set_axs(self):
        "Enables us to use a similar interface to the HorizontalPlot class for plotting to specific [row,col] indices"
        self.axs = []
        self.linescatterHandles = []
        self.coverageHandles = []
        self.rowLegend = []
        
        for colNum, sector in enumerate(self.circos.sectors):
            currentPosition = CircosPlot.START_POSITION
                        
            # Set up outer track
            outer_track = sector.add_track((currentPosition-CircosPlot.OUTER_HEIGHT, currentPosition))
            outer_track.axis(fc="black")
            
            # Determine the divisor and unit for the scale bar
            midwayPoint = int((sector.start + sector.end - 1) / 2) # midpoint of the sector
            divisor, unit = CircosPlot.INTERVALS[10 ** (len(str(midwayPoint)) - 1)] # divisor is based on the unit of the region
            
            # Determine scale bar intervals
            majorInterval = 10 ** (len(str(sector.size-1)) - 1) # majorInterval is based on the length of the region
            minorInterval = int(majorInterval / 10)
            if minorInterval < 1:
                minorInterval = 1
            
            # Determine the approximate amount of degrees this sector will take up
            sectorDegrees = (sector.size / self.plottedLength) * 360
            numTicksInThisSector = math.ceil(sectorDegrees / CircosPlot.DEGREES_PER_X_TICK) + 2 # first and last ticks are unlabelled
            if numTicksInThisSector < 3:
                numTicksInThisSector = 3 # need at least 3 ticks to offset first and last ticks not being labelled
            
            # Determine where the major ticks will be placed
            majorTickInterval = math.ceil((sector.size / numTicksInThisSector) / majorInterval) * majorInterval
            
            ongoingCount = 1
            while (sector.size / majorTickInterval < numTicksInThisSector) and (majorTickInterval > majorInterval): # coerce the number of ticks to be at least numTicksInThisSector
                majorTickInterval = math.ceil((sector.size / (numTicksInThisSector+ongoingCount)) / majorInterval) * majorInterval
                ongoingCount += 1
            
            # Set up the scale bar
            if sector.size > minorInterval:
                firstTick = sector.start - (sector.start % majorTickInterval) # first tick may preceed sector.start
                xTicks = [
                    x if x != 0 else 1 # first tick at 1 needs special handling
                    for x in range(firstTick, sector.end+1, majorTickInterval)
                    if x == 0 or (x >= sector.start and x <= sector.end) # first tick at 1 needs special handling
                ] # ensure we don't go past the end of the sector
                xLabels = [
                    None if i == 0 or i == len(xTicks)-1
                    else CircosPlot.number_formatter(x, divisor, unit)
                    for i, x in enumerate(xTicks)
                ]
                outer_track.xticks(xTicks, labels=xLabels)
            
            outer_track.xticks_by_interval(minorInterval, tick_length=1, show_label=False)
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
        
        self.axs = np.column_stack(self.axs) # gives [nrow, ncol] shape
    
    def plot(self, plotTypes, outputFileName):
        '''
        Initialises a circos plot with axes for plotting. Method is not called
        immediately to allow for customisation of optional attributes especially
        the figure width and height if not set during object initialisation.
        
        Parameters:
            plotTypes -- a list of strings indicating the types of plots to create
            outputFileName -- a string indicating the file name to save the plot to
        '''
        # Validate plot types
        if (plotTypes == None) or (not isinstance(plotTypes, list)) or (len(plotTypes) == 0):
            raise ValueError("plotTypes must be a non-empty list of plot types")
        for plotType in plotTypes:
            if plotType not in Plot.PLOT_TYPES:
                raise ValueError(f"Invalid plot type '{plotType}'; must be one of {Plot.PLOT_TYPES}")
        if len(set(plotTypes)) != len(plotTypes):
            raise ValueError("plotTypes must not contain duplicate values")
        
        # Initialise the circos figure object
        seqid2size = {
            f"{regionDict['contig']}:{regionDict['start']}-{regionDict['end']}": (regionDict['start'], regionDict['end']) # reversed plotting is not supported in Circos
            for regionDict in self.regions
        }
        self.circos = Circos(seqid2size, space = 0 if len(seqid2size) == 1 else 2,
                             end=360-self.axisSpace) # leave space for y ticks
        
        # Establish axes for each region/column and track/row
        self._set_axs()
        
        # Establish column labels
        self.colLabels = [f"{regionDict['contig']}:{regionDict['end']}-{regionDict['start']}" if regionDict['reverse'] == True
                          else f"{regionDict['contig']}:{regionDict['start']}-{regionDict['end']}" if regionDict['full'] == False
                          else regionDict['contig'] # indicate only the contig if we are plotting the full contig
                          for regionDict in self.regions]
        for label, sector in zip(self.colLabels, self.circos.sectors):
            sector.text(label, size=10)
        
        # Init values for storing data during iteration
        self.rowLabels = []
        self.rowNum = -1 # to keep track of the current row number
        
        # Build the plot
        if "line" in plotTypes or "scatter" in plotTypes:
            if self.callED != None:
                self.plot_linescatter(self.callED if "scatter" in plotTypes else None,
                                      self.callED if "line" in plotTypes else None,
                                      LINESCATTER_COLOURS[0],
                                      applyWMA=True if self.wmaSize > 1 else False,
                                      lineLabel=f"WMA $ED^{self.power}$" if self.wmaSize > 1 else f"$ED^{self.power}$",
                                      scatterLabel=f"SNP/CNV $ED^{self.power}$")
                self.rowLabels.append(f"SNP $ED^{self.power}$")
            if self.callSPLSDA != None:
                self.plot_linescatter(self.callSPLSDA[0],
                                      self.callSPLSDA[1],
                                      LINESCATTER_COLOURS[1],
                                      applyWMA=True if self.wmaSize > 1 else False,
                                      lineLabel=f"WMA $BA$" if self.wmaSize > 1 else f"$BA$",
                                      scatterLabel=f"Selected SNP/CNV",
                                      scatterShape="D", scatterFollowsLine=True, # Y values different to ED or BA
                                      lineZorder=0, scatterZorder=1, # make selected SNPs more visible
                                      dotsize=SPLSDA_DOTSIZE, dotAlpha=1, # make selected SNPs more visible
                                      integratedNCLS=self.integratedSPLSDA[0] if self.integratedSPLSDA != None else None)
                self.rowLabels.append(f"SNP $BA$")
            if self.depthED != None:
                self.plot_linescatter(self.depthED if "scatter" in plotTypes else None,
                                      self.depthED if "line" in plotTypes else None,
                                      LINESCATTER_COLOURS[0],
                                      applyWMA=True if self.wmaSize > 1 else False,
                                      lineLabel=f"WMA $ED^{self.power}$" if self.wmaSize > 1 else f"$ED^{self.power}$",
                                      scatterLabel=f"SNP/CNV $ED^{self.power}$")
                self.rowLabels.append(f"CNV $ED^{self.power}$")
            if self.depthSPLSDA != None:
                self.plot_linescatter(self.depthSPLSDA[0],
                                      self.depthSPLSDA[1],
                                      LINESCATTER_COLOURS[1],
                                      applyWMA=True if self.wmaSize > 1 else False,
                                      lineLabel=f"WMA $BA$" if self.wmaSize > 1 else f"$BA$",
                                      scatterLabel=f"Selected SNP/CNV",
                                      scatterShape="D", scatterFollowsLine=True,
                                      lineZorder=0, scatterZorder=1, # make selected CNVs more visible
                                      dotsize=SPLSDA_DOTSIZE, dotAlpha=1, # make selected CNVs more visible
                                      integratedNCLS=self.integratedSPLSDA[1] if self.integratedSPLSDA != None else None)
                self.rowLabels.append(f"CNV $BA$")
        if self.coverageNCLSDict != None:
            self.plot_coverage(self.coverageNCLSDict, self.coverageSamples)
            self.rowLabels.append("Coverage")
        if self.annotationGFF3 != None:
            self.plot_genes(self.annotationGFF3)
            self.rowLabels.append("Genes")
        
        # Add background colour to highlights regions
        if self.highlights != []:
            # For each region, ...
            for colNum, regionDict in enumerate(self.regions):
                # ... find if any highlights overlap with it
                for highlightDict in self.highlights:
                    if regionDict["contig"] == highlightDict["contig"] and \
                    is_overlapping(regionDict["start"], regionDict["end"], highlightDict["start"], highlightDict["end"]):
                        # If they do, set the background colour for this region/sector
                        sector = self.circos.sectors[colNum]
                        sector.rect(highlightDict["start"], highlightDict["end"],
                                    r_lim=(CircosPlot.CENTRE_SPACE, CircosPlot.START_POSITION),
                                    color=HIGHLIGHT_COLOUR, alpha=0.5, zorder=-1)
        
        # Set panel annotations
        minRcenter = self.axs[self.nrow-1, 0].r_center
        maxRcenter = self.axs[0, 0].r_center
        
        for i in range(self.nrow):
            # Get information for this track
            track = self.axs[i, 0]
            r = track.r_center
            
            # Get the min-max normalised scaling for the arc distance
            "We scale the arc distance (shorter at edges, longer at centre) since it makes the result look better."
            scalingFactor = (1 - minmax_norm(track.r_center, minRcenter, maxRcenter))
            if self.nrow == 1:
                scalingFactor = 1 # add the padding if there is only one track
            
            # Calculate the angle for the text label
            """Given the known radius distance from the centre of the plot, we want to 
            position the text label at a consistent arc distance away from the leftmost
            edge of the track. We do this by calculating the angle of the straight line
            with the set radius length that will pass through the point where we have
            the desired arc distance.
            """
            theta0 = math.radians(360 - self.axisSpace)
            arcAngle = (CircosPlot.ARC_DISTANCE*scalingFactor) / r
            thetaFinal = theta0 + arcAngle
            angleDegrees = math.degrees(thetaFinal)
            
            # Add the text label to the track
            self.circos.text(chr(CircosPlot.ASCII_LETTERS_START+i), r=track.r_center,
                             deg=angleDegrees)
        
        # Create the figure object
        "Necessary prior to legend addition"
        fig = self.circos.plotfig(figsize=(self.width, self.height))
        
        # Set line/scatter legend
        numLegends = 0
        if len(self.linescatterHandles) > 0:
            # Remove duplicates from the handles list
            linescatterHandles = []
            for colour, label, shape in self.linescatterHandles:
                if (colour, label, shape) not in linescatterHandles:
                    linescatterHandles.append((colour, label, shape))
            
            # Convert to Line2D handles for legend
            self.linescatterHandles = [
                Line2D([], [], color=colour, label=label, marker=shape, ms=5, ls="None")
                if shape != "line" else
                Line2D([], [], color=colour, label=label, linewidth=1)
                for colour, label, shape in linescatterHandles
            ]
            
            # Add the line/scatter legend to the circos plot
            hasED = self.callED != None or self.depthED != None
            hasSPLSDA = self.callSPLSDA != None or self.depthSPLSDA != None
            legendTitle = ("ED" if hasED else "") + \
                          ("/" if (hasED and hasSPLSDA) else "") + \
                          ("BA" if hasSPLSDA else "")
            linescatterLegend = self.circos.ax.legend(
                handles=self.linescatterHandles,
                bbox_to_anchor=CircosPlot.LEGEND_POSITIONS[numLegends][0:2],
                loc=CircosPlot.LEGEND_POSITIONS[numLegends][2],
                fontsize=8,
                title=legendTitle,
                handlelength=2
            )
            self.circos.ax.add_artist(linescatterLegend)
            numLegends += 1
        
        # Set coverage legend
        if len(self.coverageHandles) > 0:
            # Remove duplicates from the handles list
            coverageHandles = []
            for colour, label, shape in self.coverageHandles:
                if (colour, label, shape) not in coverageHandles:
                    coverageHandles.append((colour, label, shape))
            
            # Convert to Line2D handles for legend
            self.coverageHandles = [
                Line2D([], [], color=colour, label=label, linewidth=1, ls=shape)
                for colour, label, shape in coverageHandles
            ]
            
            # Add the coverage legend to the circos plot
            legendTitle = "Median-norm\nCoverage"
            coverageLegend = self.circos.ax.legend(
                handles=self.coverageHandles,
                bbox_to_anchor=CircosPlot.LEGEND_POSITIONS[numLegends][0:2],
                loc=CircosPlot.LEGEND_POSITIONS[numLegends][2],
                fontsize=8,
                title=legendTitle,
                handlelength=2
            )
            self.circos.ax.add_artist(coverageLegend)
            numLegends += 1
        
        # Set panel legend
        self.rowHandles = [
            Patch(label=f"{chr(CircosPlot.ASCII_LETTERS_START+i)}: {label}") # produces 'a: SNP ED^4', etc.
            for i, label in enumerate(self.rowLabels)
        ]
        legendTitle = "Panels"
        rowLabelLegend = self.circos.ax.legend(
            handles=self.rowHandles,
            bbox_to_anchor=CircosPlot.ROW_LABEL_POSITION[0:2],
            loc=CircosPlot.ROW_LABEL_POSITION[2],
            fontsize=6,
            title_fontsize=8,
            title=legendTitle,
            handlelength=0, handletextpad=0, # hide the handles
            alignment="left" if len(self.rowHandles) == 1 else "center",
            ncols=2
        )
        self.circos.ax.add_artist(rowLabelLegend)
        
        # Save the figure
        fig.savefig(outputFileName, dpi=300, pad_inches=1.5)
    
    @property
    def trackHeight(self):
        return (100 - (
            (CircosPlot.TRACK_GAP * self.nrow) +
            (100 - CircosPlot.START_POSITION) +
            CircosPlot.CENTRE_SPACE)
        ) / self.nrow # height of each track
    
    def _format_y_ticks(self, maxY):        
        yticks = np.linspace(0, maxY, CircosPlot.NUM_Y_TICKS)
        if maxY < 1:
            yticks = [ CircosPlot.truncate(x, 2) for i, x in enumerate(yticks) ]
        else:
            yticks = [ int(x) for x in yticks ]
        ylabels = [ str(x) for x in yticks ]
        return yticks, ylabels
    
    def plot_linescatter(self, scatterNCLS, lineNCLS, colours, applyWMA=True,
                         lineLabel="WMA", scatterLabel="Values", scatterShape="o",
                         lineZorder=1, scatterZorder=0, linewidth=1,
                         dotsize=3, dotAlpha=0.5,
                         integratedNCLS=None, scatterFollowsLine=False):
        '''
        Plots the data for a line or scatter plot.
        
        Parameters:
            scatterNCLS -- a WindowedNCLS object with statistical values
                            queryable by contigID and start/end positions;
                            used for scatter plots OR None if not plotting
                            scatter points
            lineNCLS -- a WindowedNCLS object with statistical values
                        queryable by contigID and start/end positions;
                        used for line plots OR None if not plotting
                        a line
            colours -- a list of two colours to use for the line (first) and scatter
                       points (second)
            applyWMA -- (OPTIONAL) a boolean indicating whether to apply WMA
                        smoothing to the line; default is True
            lineLabel -- (OPTIONAL) a string indicating the label for the line;
                          default is "WMA"
            scatterLabel -- (OPTIONAL) a string indicating the label for the scatter
                            points; default is "Values"
            scatterShape -- (OPTIONAL) a string indicating the shape of the scatter
                            points; default is "o" for standard dots
            lineZorder -- (OPTIONAL) an integer indicating the z-order of the line;
                          default is 1
            scatterZorder -- (OPTIONAL) an integer indicating the z-order of the scatter
                             points; default is 0
            linewidth -- (OPTIONAL) an integer indicating the width of the line;
                         default is 1
            dotsize -- (OPTIONAL) an integer indicating the size of the scatter
                        points; default is 3
            dotAlpha -- (OPTIONAL) a float indicating the alpha value for the scatter
                        points; default is 0.5
            integratedNCLS -- (OPTIONAL) a WindowedNCLS object with integrated sPLSDA
                              values indicating the selected SNPs and CNVs when assessing
                              selected features from 'call' and 'depth' simultaneously;
                              used for scatter plots OR None if not plotting integrated
                              values
            scatterFollowsLine -- (OPTIONAL) a boolean indicating whether the scatter y values
                                should follow the line y values; default is False.
        '''
        # Validate that one of either scatterNCLS or lineNCLS is provided
        if scatterNCLS is None and lineNCLS is None:
            raise ValueError("At least one of scatterNCLS or lineNCLS must be provided")
        
        # Validate that colours is a list of two colours
        if colours is None or len(colours) != 2:
            raise ValueError("colours must be a list of two colours")
        
        self.rowNum += 1 # increment row number for plotting
        
        # Derive y limits from the maximum Y value across all regions
        maxY = 0 # to set y limits at end
        for colNum, regionDict in enumerate(self.regions):
            if scatterNCLS != None and regionDict["contig"] in scatterNCLS.contigs:
                x, y = self.scatter(scatterNCLS, regionDict["contig"], regionDict["start"], regionDict["end"])
                if y.size != 0 and not scatterFollowsLine:
                    maxY = max(maxY, max(y))
            if lineNCLS != None and regionDict["contig"] in lineNCLS.contigs:
                x, smoothedY = self.line(lineNCLS, regionDict["contig"], regionDict["start"], regionDict["end"],
                                         applyWMA=applyWMA, buffer=SMOOTHING_BUFFER)
                if smoothedY.size != 0:
                    maxY = max(maxY, max(smoothedY))
        
        # Add y-axis
        yticks, ylabels = self._format_y_ticks(maxY)
        self.axs[self.rowNum, 0].yticks(yticks, ylabels, vmin=0, vmax=maxY, side="left")
        
        # Plot each region
        for colNum, regionDict in enumerate(self.regions):
            # Decompose dictionary into variables
            contigID = regionDict["contig"]
            start = regionDict["start"]
            end = regionDict["end"]
            reverse = regionDict["reverse"]
            
            # Adjust start and end for one-based indexing
            # Plot line (if applicable)
            if lineNCLS != None and contigID in lineNCLS.contigs:
                lineX, smoothedY = self.line(lineNCLS, contigID, start, end,
                                             applyWMA=applyWMA, buffer=SMOOTHING_BUFFER)
                if smoothedY.size != 0:
                    self.axs[self.rowNum, colNum].line(lineX, smoothedY, vmax=maxY,
                                                       color=colours[0],
                                                       linewidth=linewidth,
                                                       zorder=lineZorder,
                                                       label=lineLabel)
                    self.linescatterHandles.append([colours[0], lineLabel, "line"])
            
            # Plot scatter (if applicable)
            if scatterNCLS != None and contigID in scatterNCLS.contigs:
                x, y = self.scatter(scatterNCLS, contigID, start, end)
                if y.size != 0:
                    if scatterFollowsLine:
                        if lineNCLS == None:
                            raise ValueError("scatterFollowsLine is True but no lineNCLS provided; cannot interpolate scatter points")
                        y = Plot.match_scatter_to_line(lineX, smoothedY, x)
                    self.axs[self.rowNum, colNum].scatter(x, y, vmax=maxY,
                                                          color=colours[1],
                                                          s=dotsize, alpha=dotAlpha, marker=scatterShape,
                                                          zorder=scatterZorder,
                                                          label=scatterLabel)
                    self.linescatterHandles.append([colours[1], scatterLabel, scatterShape])
            
            # Plot integrated sPLSDA values (if applicable)
            if integratedNCLS is not None and contigID in integratedNCLS.contigs:
                x, y = self.scatter(integratedNCLS, contigID, start, end)
                if y.size != 0:
                    if scatterFollowsLine:
                        if lineNCLS == None:
                            raise ValueError("scatterFollowsLine is True but no lineNCLS provided; cannot interpolate scatter points")
                        y = Plot.match_scatter_to_line(lineX, smoothedY, x)
                    self.axs[self.rowNum, colNum].scatter(x, y, vmax=maxY,
                                                          color=INTEGRATED_AESTHETICS[0], alpha=1, # always alpha=1
                                                          s=dotsize+1, marker=INTEGRATED_AESTHETICS[1],
                                                          zorder=2, # integrated points are on top
                                                          label=INTEGRATED_AESTHETICS[2])
                    self.linescatterHandles.append([INTEGRATED_AESTHETICS[0], INTEGRATED_AESTHETICS[2], INTEGRATED_AESTHETICS[1]])
    
    def plot_histogram(self, windowedNCLS):
        '''
        Plots the data for a histogram.
        
        Parameters:
            windowedNCLS -- a WindowedNCLS object with statistical values
                            queryable by contigID and start/end positions
        '''
        self.rowNum += 1 # increment row number for plotting
        
        # Get the maximum Y value across all regions
        maxY = 0
        for regionDict in self.regions:
            if regionDict["contig"] in windowedNCLS.contigs:
                x, y = self.histogram(windowedNCLS, regionDict["contig"], regionDict["start"], regionDict["end"])
                if y.size != 0:
                    maxY = max(maxY, max(y))
        
        # Add y-axis
        yticks, ylabels = self._format_y_ticks(maxY)
        self.axs[self.rowNum, 0].yticks(yticks, ylabels, vmin=0, vmax=maxY, side="left")
        
        # Plot each region
        for colNum, regionDict in enumerate(self.regions):
            if regionDict["contig"] in windowedNCLS.contigs:
                x, y = self.histogram(windowedNCLS, regionDict["contig"], regionDict["start"], regionDict["end"])
                self.axs[self.rowNum, colNum].bar(np.clip(x, regionDict["start"], regionDict["end"]), y, vmax=maxY,
                                                  width=self.binSize,
                                                  color="grey",
                                                  align="edge")
    
    def plot_genes(self, gff3Obj):
        '''
        Plots the data for gene models.
        
        Parameters:
            gff3Obj -- a GFF3Graph object containing gene models with ncls indexing
        '''
        self.rowNum += 1 # increment row number for plotting
        
        for colNum, regionDict in enumerate(self.regions):
            # Decompose dictionary into variables
            contigID = regionDict["contig"]
            start = regionDict["start"]
            end = regionDict["end"]
            reverse = regionDict["reverse"]
            
            # Get longest isoform for each gene in this region
            mrnaFeatures = self.genes(gff3Obj, contigID, start, end)
            
            # Plot each feature
            for mrnaFeature in mrnaFeatures:
                # Get start and end positions with truncation if they extend outside the region
                newStart = max(start, mrnaFeature.start)
                newEnd = min(end, mrnaFeature.end)
                
                # Format a feature amenable to pycirclize plotting
                feature = SeqFeature(
                    location=SimpleLocation(
                        ExactPosition(newStart),
                        ExactPosition(newEnd),
                        strand=-1 if mrnaFeature.strand == "-" else 1
                    ),
                    type="mRNA",
                    id=mrnaFeature.ID,
                    qualifiers={"ID": [mrnaFeature.ID]}
                )
                
                # Plot the feature
                self.axs[self.rowNum, colNum].genomic_features(
                    [feature], plotstyle="arrow",
                        fc=GENE_COLOURS[0])
    
    def plot_coverage(self, depthNCLSDict, samples,
                      linewidth=1):
        '''
        Parameters:
            depthNCLSDict -- a dictionary with structure like:
                            {
                                "group1": {
                                    "sample1": EDNCLS,
                                    "sample2": EDNCLS,
                                    ...
                                },
                                "group2": { ... }
                            }
            samples -- a list of strings indicating the sample names to plot individual
                       lines for
            linewidth -- OPTIONAL; an integer value indicating the width of the line plot (default=1)
        '''
        self.rowNum += 1 # increment row number for plotting
        
        # Get the maximum Y value across all regions
        maxY = 0
        for regionDict in self.regions:
            coverageData = self.coverage(depthNCLSDict, samples, regionDict["contig"], regionDict["start"], regionDict["end"])
            for group in ["group1", "group2"]:
                groupData = coverageData[group]
                if groupData != None:
                    maxY = max(maxY, np.percentile(groupData["median"], 90)) # 90th percentile to trim outliers
        if maxY == 0:
            maxY = 1 # this can occur if all values are 0
        
        # Add y-axis
        yticks, ylabels = self._format_y_ticks(maxY)
        self.axs[self.rowNum, 0].yticks(yticks, ylabels, vmin=0, vmax=maxY, side="left")
        
        # Plot each region
        for colNum, regionDict in enumerate(self.regions):
            # Decompose dictionary into variables
            contigID = regionDict["contig"]
            start = regionDict["start"]
            end = regionDict["end"]
            reverse = regionDict["reverse"]
            
            # Get coverage data for this region
            coverageData = self.coverage(depthNCLSDict, samples, contigID, start, end) # keys: group1, group2, [*samples]
            
            # Plot each group
            extendStart = False
            extendEnd = False
            for group in ["group1", "group2"]:
                groupData = coverageData[group] # keys: x, q1, median, q3
                
                # Skip if there are no values
                if groupData == None:
                    print(f"WARNING: No coverage values found for region '{contigID, start, end}' for '{group}'")
                    continue
                
                # Extend tails for better visualisation
                if groupData["x"][0] < start:
                    groupData["x"][0] = start
                else:
                    groupData["x"] = np.concatenate(([start], groupData["x"]))
                    groupData["median"] = np.concatenate(([groupData["median"][0]], groupData["median"]))
                    groupData["q1"] = np.concatenate(([groupData["q1"][0]], groupData["q1"]))
                    groupData["q3"] = np.concatenate(([groupData["q3"][0]], groupData["q3"]))
                    extendStart = True
                
                if groupData["x"][-1] > end:
                    groupData["x"][-1] = end
                else:
                    groupData["x"] = np.concatenate((groupData["x"], [end]))
                    groupData["median"] = np.concatenate((groupData["median"], [groupData["median"][-1]]))
                    groupData["q1"] = np.concatenate((groupData["q1"], [groupData["q1"][-1]]))
                    groupData["q3"] = np.concatenate((groupData["q3"], [groupData["q3"][-1]]))
                    extendEnd = True
                
                # Cap values at maxY
                groupData["median"] = np.clip(groupData["median"], 0, maxY)
                groupData["q1"] = np.clip(groupData["q1"], 0, maxY)
                groupData["q3"] = np.clip(groupData["q3"], 0, maxY)
                
                # Plot median and Q1/Q3 lines
                self.axs[self.rowNum, colNum].fill_between(
                    groupData["x"], groupData["q1"], groupData["q3"], vmax=maxY,
                    alpha = 0.5, color=COVERAGE_COLOURS[0] if group == "group1" else COVERAGE_COLOURS[1],
                    label="_nolegend_", zorder=0)
                self.axs[self.rowNum, colNum].line(
                    groupData["x"], groupData["median"], vmax=maxY,
                    color=COVERAGE_COLOURS[0] if group == "group1" else COVERAGE_COLOURS[1],
                    linewidth=linewidth, zorder=1)
                
                self.coverageHandles.append([COVERAGE_COLOURS[0] if group == "group1" else COVERAGE_COLOURS[1],
                                             "Group 1" if group == "group1" else "Group 2", "solid"])
            
            # Plot individual samples
            for sampleIndex, sample in enumerate(samples):
                "x can be reused from the group plot"
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
                    groupData["x"], y, vmax=maxY, 
                    color=lineColour, linestyle=lineType,
                    linewidth=linewidth, zorder=2)
                self.coverageHandles.append([lineColour, sample, lineType])
    
    def __repr__(self):
        return (f"CircosPlot(regions={self.regions}" ) + \
               (", callED=" + "True" if self.callED is not None else "False") + \
               (", depthED=" + "True" if self.depthED is not None else "False") + \
               (", callSPLSDA=" + "True" if self.callSPLSDA is not None else "False") + \
               (", depthSPLSDA=" + "True" if self.depthSPLSDA is not None else "False") + \
               (", coverage=" + "True" if self.coverageNCLSDict is not None else "False") + \
               (", genes=" + "True" if self.annotationGFF3 is not None else "False") + \
               (f", power={self.power}, wmaSize={self.wmaSize}, width={self.width}, height={self.height})")
