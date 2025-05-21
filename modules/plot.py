import math
import numpy as np
import pandas as pd

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
    
    def __init__(self, resultType, measurementType, plotType, regions):
        # Mandatory parameters
        self.resultType = resultType
        self.measurementType = measurementType
        self.plotType = plotType
        self.regions = regions
        
        # Optional parameters
        self.wmaSize = 5
        self.binSize = 100000
        self.binThreshold = 0.4
        self.width = None
        self.height = None
    
    @property
    def resultType(self):
        return self._resultType
    
    @resultType.setter
    def resultType(self, valueList):
        if not isinstance(valueList, list):
            raise TypeError("resultType must be a list")
        if len(valueList) < 1 or len(valueList) > 2:
            raise ValueError("resultType must be a list of length 1 or 2")
        for value in valueList:
            if value not in Plot.RESULT_TYPES:
                raise ValueError(f"resultType must be one of {Plot.RESULT_TYPES}")
        
        self._resultType = valueList
    
    @property
    def measurementType(self):
        return self._measurementType
    
    @measurementType.setter
    def measurementType(self, valueList):
        if not isinstance(valueList, list):
            raise TypeError("measurementType must be a list")
        if len(valueList) < 1 or len(valueList) > 2:
            raise ValueError("measurementType must be a list of length 1 or 2")
        for value in valueList:
            if value not in Plot.MEASUREMENT_TYPES:
                raise ValueError(f"measurementType must be one of {Plot.MEASUREMENT_TYPES}")
        
        self._measurementType = valueList
    
    @property
    def plotType(self):
        return self._plotType
    
    @plotType.setter
    def plotType(self, valueList):
        if not isinstance(valueList, list):
            raise TypeError("plotType must be a list")
        if len(valueList) < 1 or len(valueList) > len(Plot.PLOT_TYPES):
            raise ValueError("plotType must be a list of length 1 or 2")
        for value in valueList:
            if value not in Plot.PLOT_TYPES:
                raise ValueError(f"plotType must be one of {Plot.PLOT_TYPES}")
        
        self._plotType = valueList
    
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
        return len(self.resultType) * len(self.measurementType) * \
            len(self.plotType)
                if not all(["scatter" in self.plotType, 
                            "line" in self.plotType ]) 
            else len(self.plotType) - 1 # scatter and line are in the same row
    
    @property
    def width(self):
        if self._width is None:
            return Plot.STANDARD_DIMENSION * self.ncol
        return self._width
    
    @width.setter
    def width(self, value):
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
        if not isinstance(value, int):
            raise TypeError("height must be an integer")
        if value < 1:
            raise ValueError(f"height must be >= 1")
        
        self._height = value
    
    def scatter(self, windowedNCLS):
        '''
        Yields data suited for scatter plotting of WindowedNCLS values.
        
        Parameters:
            windowedNCLS -- a WindowedNCLS object with statistical values
                            queryable by contigID and start/end positions
        
        Yields:
            contigID -- a string indicating the contig ID
            start -- an integer indicating the start position of the region
            end -- an integer indicating the end position of the region
            reverse -- a boolean indicating if plotting should be reversed
            x -- a numpy array of the x values (positions)
            y -- a numpy array of the y values (statistical values)
        '''
        for contigID, start, end, reverse in self.regions:
            # Get values within this region
            regionValues = windowedNCLS.find_overlap(contigID, start, end)
            x, y = [], []
            for pos, _, ed in regionValues:
                x.append(pos)
                y.append(ed)
            x = np.array(x)
            y = np.array(y)
            
            # Yield the original values
            yield contigID, start, end, reverse, x, y
    
    def line(self, windowedNCLS):
        '''
        Yields data suited for line plotting of WindowedNCLS values.
        
        Parameters:
            windowedNCLS -- a WindowedNCLS object with statistical values
                            queryable by contigID and start/end positions
        
        Yields:
            contigID -- a string indicating the contig ID
            start -- an integer indicating the start position of the region
            end -- an integer indicating the end position of the region
            reverse -- a boolean indicating if plotting should be reversed
            x -- a numpy array of the x values (positions)
            smoothedY -- a pandas Series of the smoothed y values (statical
                         value) values OR the original y values if smoothing
                         was not possible (i.e., not enough data points)
        '''
        for contigID, start, end, reverse, x, y in self.scatter(windowedNCLS):
            smoothedY = WMA(y, self.wmaSize)
            if smoothedY is None:
                print(f"WARNING: region '{contigID, start, end}' has too few data points to apply WMA smoothing")
                smoothedY = y
            yield contigID, start, end, reverse, x, smoothedY
    
    def histogram(self, windowedNCLS):
        '''
        Yields data suited for histogram plots of binned ED values.
        
        Parameters:
            windowedNCLS -- a WindowedNCLS object with statistical values
                            queryable by contigID and start/end positions
        Yields:
            contigID -- a string indicating the contig ID
            start -- an integer indicating the start position of the region
            end -- an integer indicating the end position of the region
            reverse -- a boolean indicating if plotting should be reversed
            x -- a numpy array of the x values (positions)
            y -- a numpy array of the y values (binned statistical values)
        '''
        for contigID, start, end, reverse in self.regions:
            y = bin_values(windowedNCLS.find_overlap(contigID, start, end),
                           start, end, self.binSize, self.binThreshold)
            x = np.array([ (i * self.binSize) + start for i in range(len(y)) ])
            yield contigID, start, end, reverse, x, y
    
    def genes(self, gff3Obj):
        '''
        Yields data suited for histogram plots of binned ED values.
        
        Parameters:
            gff3Obj -- a GFF3 object with ncls methods to locate geneFeatures
                       within given regions
        Yields:
            contigID -- a string indicating the contig ID
            start -- an integer indicating the start position of the region
            end -- an integer indicating the end position of the region
            reverse -- a boolean indicating if plotting should be reversed
            mrnaFeatures -- a list of mRNA features within the region
        '''
        for contigID, start, end, reverse in self.regions:
            # Get longest isoform for each gene in this region
            geneFeatures = gff3Obj.ncls_finder(start, end, "contig", contigID)
            mrnaFeatures = [
                GFF3.longest_isoform(geneFeature)
                for geneFeature in geneFeatures
                if hasattr(geneFeature, "mRNA")
            ]
            yield contigID, start, end, reverse, mrnaFeatures
    
    def coverage(self, depthNCLSDict, samples):
        for colNum, (contigID, start, end, reverse) in enumerate(regions):
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
            
            # Yield the data
            yield contigID, start, end, reverse, coverageData
    
    def __repr__(self):
        return f"Plot(resultType={self.resultType}, measurementType={self.measurementType}, plotType={self.plotType})"

class HorizontalPlot(Plot):
    def __init__(self, resultType, measurementType, plotType, regions):
        super().__init__(resultType, measurementType, plotType, regions)
    
    def plot(self):
        # Implement the plotting logic here
        pass
    
    def __repr__(self):
        return f"HorizontalPlot(resultType={self.resultType}, measurementType={self.measurementType}, plotType={self.plotType})"

class CircosPlot(Plot):
    def __init__(self, resultType, measurementType, plotType, regions):
        super().__init__(resultType, measurementType, plotType, regions)
    
    def plot(self):
        # Implement the plotting logic here
        pass
