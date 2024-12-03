import os, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.patches import Rectangle
from itertools import product

from .gff3 import GFF3

YLIM_HEADSPACE = 0.1 # proportion of ylim to add to the top of the plot
SAMPLE_AESTHETICS = [["#000000", "dotted"], ["#002D7E", "dashed"], ["#ECE45A", "dashdot"]]
NUM_SAMPLE_LINES = len(SAMPLE_AESTHETICS) # for validation

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
        if stat >= binThreshold:
            binIndex = (pos-start) // binSize
            histo[binIndex] += 1
    return histo

def linescatter(axs, rowNum, edNCLS, regions, wmaSize, line, scatter, 
                power, outputDirectory, plotScalebar, fileSuffix,
                linewidth=1, dotsize=3):
    '''
    Parameters:
        axs -- a list of matplotlib.pyplot Axes objects to plot to
        rowNum -- an integer value indicating the row index to plot to
        edNCLS -- an EDNCLS object
        regions -- a list of lists containing three values: [contigID, start, end]
        wmaSize -- an integer value indicating the number of previous values to consider
                   during weighted moving average calculation
        line -- a boolean value indicating whether to plot a line
        scatter -- a boolean value indicating whether to plot scatter points
        power -- an integer value indicating what power statistical values were raised to
        outputDirectory -- a string indicating the directory to save the output TSV data
        plotScalebar -- a boolean value indicating whether to plot a scalebar on the X axis
        fileSuffix -- a string to append to the output file name
        linewidth -- OPTIONAL; an integer value indicating the width of the line plot (default=1)
        dotsize -- OPTIONAL; an integer value indicating the size of the dots (default=3)
    '''
    if line == False and scatter == False:
        raise ValueError("At least one of 'line' or 'scatter' must be True")
    
    # Plot each region
    maxY = 0 # to set y limits at end
    for colNum, (contigID, start, end) in enumerate(regions):
        # Get values within this region
        regionValues = edNCLS.find_overlap(contigID, start, end)
        x, y = [], []
        for pos, _, ed in regionValues:
            x.append(pos)
            y.append(ed)
        x = np.array(x)
        y = np.array(y)
        
        # Smooth the y values
        if line == True:
            smoothedY = WMA(y, wmaSize)
            if smoothedY is None:
                print(f"WARNING: '{contigID, start, end}' has too few data points to apply WMA smoothing")
                smoothedY = y
        
        # Get the maximum Y value for this region
        if scatter == True:
            maxY = max(maxY, max(y))
        else:
            maxY = max(maxY, max(smoothedY))
        
        # Set xlim
        axs[rowNum, colNum].set_xlim(start, end)
        
        # Turn off ytick labels if not the first column
        if colNum > 0:
            axs[rowNum, colNum].set_yticklabels([])
        
        # Plot scatter (if applicable)
        if scatter == True:
            axs[rowNum, colNum].scatter(x, y, color="red", s=dotsize, alpha=0.5,
                                        zorder=0)
        
        # Plot line (if applicable)
        if line == True:
            axs[rowNum, colNum].plot(x, smoothedY, linewidth=linewidth,
                                     zorder=1)
        
        # Set up scale bar if this is the last row
        if plotScalebar == True:
            scalebar(axs, rowNum, colNum, start, end)
        # Otherwise, turn off x labels
        else:
            axs[rowNum, colNum].set_xticklabels([])
            axs[rowNum, colNum].locator_params(axis='x', nbins=4) # use less ticks; avoid clutter
        
        # Derive our output file name for TSV data
        fileOut = os.path.join(outputDirectory, f"{contigID}.{start}-{end}.{fileSuffix}_line.tsv")
        if os.path.isfile(fileOut):
            print(f"WARNING: Line/scatter plot data for '{contigID, start, end}' already exists as '{fileOut}'; won't overwrite")
            continue
        else:
            with open(fileOut, "w") as fileOutTSV:
                if line == True:
                    fileOutTSV.write("contigID\tposition\ted\tsmoothed_ed\n")
                    for xVal, yVal, smoothedYVal in zip(x, y, smoothedY):
                        fileOutTSV.write(f"{contigID}\t{xVal}\t{yVal}\t{smoothedYVal}\n")
                else:
                    fileOutTSV.write("contigID\tposition\ted\n")
                    for xVal, yVal in zip(x*1000000, y):
                        fileOutTSV.write(f"{contigID}\t{xVal}\t{yVal}\n")
    
    # Set y limits
    for colNum in range(len(regions)):
        axs[rowNum, colNum].set_ylim(0, maxY + (maxY * YLIM_HEADSPACE))

def histogram(axs, rowNum, edNCLS, regions, binSize, binThreshold, power, outputDirectory, plotScalebar):
    '''
    Parameters:
        axs -- a list of matplotlib.pyplot Axes objects to plot to
        rowNum -- an integer value indicating the row index to plot to
        edNCLS -- an EDNCLS object
        regions -- a list of lists containing three values: [contigID, start, end]
        binSize -- an integer value indicating the size of the bins/windows
        binThreshold -- an integer value indicating the threshold for counting a variant
                        within a bin/window
        power -- an integer value indicating what power statistical values were raised to
        outputDirectory -- a string indicating the directory to save the output TSV data
        plotScalebar -- a boolean value indicating whether to plot a scalebar on the X axis
    '''
    # Get the maximum Y value across all regions
    maxY = 0
    for contigID, start, end in regions:
        y = bin_values(edNCLS.find_overlap(contigID, start, end),
                       start, end, binSize, binThreshold)
        maxY = max(maxY, max(y))
    
    # Plot each region
    for colNum, (contigID, start, end) in enumerate(regions):
        y = bin_values(edNCLS.find_overlap(contigID, start, end),
                       start, end, binSize, binThreshold)
        x = np.array([ (i * binSize) + start for i in range(len(y)) ])
        
        # Set limits
        axs[rowNum, colNum].set_xlim(start, end)
        axs[rowNum, colNum].set_ylim(0, maxY + (maxY * YLIM_HEADSPACE))
        
        # Turn off ytick labels if not the first column
        if colNum > 0:
            axs[rowNum, colNum].set_yticklabels([])
        
        # Plot bars
        axs[rowNum, colNum].bar(x, y, align="edge", width=binSize)
        
        # Set up scale bar if this is the last row
        if plotScalebar == True:
            scalebar(axs, rowNum, colNum, start, end)
        # Otherwise, turn off x labels
        else:
            axs[rowNum, colNum].set_xticklabels([])
            axs[rowNum, colNum].locator_params(axis='x', nbins=4) # use less ticks; avoid clutter
        
        # Derive our output file name for TSV data
        fileOut = os.path.join(outputDirectory, f"{contigID}.{start}-{end}.histo.tsv")
        if os.path.isfile(fileOut):
            print(f"WARNING: Histogram plot data for '{contigID, start, end}' already exists as '{fileOut}'; won't overwrite")
            continue
        else:
            with open(fileOut, "w") as fileOutTSV:
                fileOutTSV.write(f"contigID\twindow_num\twindow_start\tnum_variants >= {binThreshold}\n")
                for xVal, yVal in zip(x, y):
                    fileOutTSV.write(f"{contigID}\t{xVal}\t{(xVal * binSize) + start}\t{yVal}\n")

def genes(fig, axs, rowNum, gff3Obj, regions, power, plotScalebar):
    '''
    Parameters:
        fig -- the matplotlib.pyplot Figure object that axes correspond to
        axs -- a list of matplotlib.pyplot Axes objects to plot to
        rowNum -- an integer value indicating the row index to plot to
        gff3Obj -- a GFF3 class object from gff3.py in this repository
        regions -- a list of lists containing three values: [contigID, start, end]
        
        power -- an integer value indicating what power statistical values were raised to
        plotScalebar -- a boolean value indicating whether to plot a scalebar on the X axis
    Returns:
        genePltList -- a list of matplotlib.pyplot objects containing the genes plot data
                       per region
    '''
    fig.canvas.draw() # need to draw the figure to get the renderer
    SPACING = 0.1
    ARROW_PROPORTION = 0.005 # proportion of subplot width to use for arrow size
    
    for colNum, (contigID, start, end) in enumerate(regions):
        ARROW_SIZE = math.ceil(ARROW_PROPORTION * (end - start)) # arrow size is proportional to region width
        
        # Get longest isoform for each gene in this region
        geneFeatures = gff3Obj.ncls_finder(start, end, "contig", contigID)
        mrnaFeatures = [
            GFF3.longest_isoform(geneFeature)
            for geneFeature in geneFeatures
            if hasattr(geneFeature, "mRNA")
        ]
        
        # Turn off y labels
        axs[rowNum, colNum].set_yticklabels([])
        axs[rowNum, colNum].set_yticks([])
        
        # Set xlim to allow renderer to calculate text box positions correctly
        axs[rowNum, colNum].set_xlim(start, end)
        
        # Plot each gene in non-overlapping lanes
        lanes = []
        for mrnaFeature in mrnaFeatures:
            # Resolve overlaps
            placed = False
            laneNum = 0
            for lane in lanes:
                if (mrnaFeature.start - ARROW_SIZE) > lane[-1]:
                    placed = True
                    break
                laneNum += 1
            if not placed:
                lane = [] # start a new lane
                lanes.append(lane)
            
            # Plot intron line
            "Intron line, when at the bottommost layer, can just be the length of the gene"
            axs[rowNum, colNum].plot([max(mrnaFeature.start, start), min(mrnaFeature.end, end)], # truncated to region
                                     [laneNum + SPACING + (1-SPACING)/2]*2, # y values for start and end
                                     color="black", linewidth=1,
                                     zorder=0)
            
            # Plot gene directionality
            if mrnaFeature.strand in ["+", "-"]:
                if mrnaFeature.strand == "+":
                    lastExon = max(mrnaFeature.exon, key=lambda x: x.end)
                    
                    arrowVertices = [
                        (lastExon.end, laneNum + SPACING), # bottom
                        (lastExon.end, laneNum + 1), # top
                        (lastExon.end + ARROW_SIZE, laneNum + SPACING + (1-SPACING)/2) # tip
                    ]
                else:
                    lastExon = min(mrnaFeature.exon, key=lambda x: x.start)
                    
                    arrowVertices = [
                        (lastExon.start, laneNum + SPACING), # bottom
                        (lastExon.start, laneNum + 1), # top
                        (lastExon.start - ARROW_SIZE, laneNum + SPACING + (1-SPACING)/2) # tip
                    ]
                axs[rowNum, colNum].add_patch(
                    plt.Polygon(arrowVertices, closed=True, fill=True,
                                edgecolor="black",
                                facecolor="black",
                                zorder=1) # above exon boxes
                )
            
            # Plot exon boxes
            exonCoords = [
                (
                    max(exonFeature.start, start),
                    min(exonFeature.end, end)
                )
                for exonFeature in mrnaFeature.exon
            ]
            axs[rowNum, colNum].broken_barh(
                [
                    (exonStart, exonEnd - exonStart)
                    for exonStart, exonEnd in exonCoords
                    if exonStart < exonEnd # this occurs if the exon exists outside of the specified region
                ],
                (laneNum+SPACING, 1-SPACING),
                facecolors="dodgerblue", edgecolor="black",
                zorder=2 # above gene directionality
            )
            
            # Plot CDS boxes
            cdsCoords = [
                (
                    max(cdsFeature.start, start),
                    min(cdsFeature.end, end)
                )
                for cdsFeature in mrnaFeature.CDS
            ]
            axs[rowNum, colNum].broken_barh(
                [
                    (cdsStart, cdsEnd - cdsStart)
                    for cdsStart, cdsEnd  in cdsCoords
                    if cdsStart < cdsEnd # this occurs if the exon exists outside of the specified region
                ],
                (laneNum+SPACING, 1-SPACING),
                facecolors="coral", edgecolor="black",
                zorder=3 # above exon boxes and gene directionality
            )
            
            # Plot gene name
            geneName = mrnaFeature.ID
            textBox = axs[rowNum, colNum].text(
                mrnaFeature.end + ARROW_SIZE, # x position
                laneNum + SPACING + (1-SPACING)/2, # y position
                geneName, # text
                horizontalalignment="left", verticalalignment="center", # alignment
                fontsize=8,
                zorder=4 # above everything else
            )
            transf = axs[rowNum, colNum].transData.inverted()
            bb = textBox.get_window_extent(renderer = fig.canvas.renderer)
            bb_datacoords = bb.transformed(transf)
            
            # Store the rightmost x value for this lane
            lane.append(bb_datacoords.x1)
        
        # Set up scale bar if this is the last row
        if plotScalebar == True:
            scalebar(axs, rowNum, colNum, start, end)
        # Otherwise, turn off x labels
        else:
            axs[rowNum, colNum].set_xticklabels([])
            axs[rowNum, colNum].locator_params(axis='x', nbins=4) # use less ticks; avoid clutter
    
    # Set ylim to the maximum number of lanes
    axs[rowNum, colNum].set_ylim(0, len(lanes)+SPACING)

def coverage(axs, rowNum, depthNCLSDict, regions, samples, plotScalebar, linewidth=1):
    '''
    Parameters:
        axs -- a list of matplotlib.pyplot Axes objects to plot to
        rowNum -- an integer value indicating the row index to plot to
        depthNCLSDict -- a dictionary with structure like:
                         {
                             "bulk1": {
                                 "sample1": EDNCLS,
                                 "sample2": EDNCLS,
                                 ...
                            },
                             "bulk2": { ... }
                         }
        regions -- a list of lists containing three values: [contigID, start, end]
        samples -- a list of strings indicating the sample names to plot individual
                   lines for
        plotScalebar -- a boolean value indicating whether to plot a scalebar on the X axis
        linewidth -- OPTIONAL; an integer value indicating the width of the line plot (default=1)
    '''
    # Plot each region
    minY = 0 # to set y limits at end
    maxY = 0
    for colNum, (contigID, start, end) in enumerate(regions):
        # Plot each bulk
        for bulk, sampleDict in depthNCLSDict.items():
            # Get the average and upper/lower quantiles of each bulk's coverage values
            bulkValues = []
            for sample, edNCLS in sampleDict.items():
                # Get values within this region
                overlappingBins = list(edNCLS.find_overlap(contigID, start, end))
                y = [ stat for _, _, stat in overlappingBins ]
                bulkValues.append(y)
                # Get the minimum Y value for this region
                minY = min(minY, min(y))
            
            # Get the median and Q1/Q3 quantiles
            bulkValues = np.array(bulkValues)
            q1, median, q3 = np.percentile(bulkValues, [25, 50, 75], axis=0)
            
            # Figure out the x values
            x = []
            for index, (windowStart, windowEnd, _) in enumerate(overlappingBins):
                x.append(windowStart + (windowEnd - windowStart+1)/2)
            
            # Plot median and Q1/Q3 lines
            axs[rowNum, colNum].plot(x, median, linewidth=linewidth,
                                     color="palevioletred" if bulk == "bulk1" else "mediumseagreen")
            axs[rowNum, colNum].fill_between(x, q1, q3, alpha = 0.5,
                                     color="palevioletred" if bulk == "bulk1" else "mediumseagreen")
            
            # Get the maximum Y value for this region
            maxY = max(maxY, np.percentile(median, 90)) # 90th percentile to trim outliers
        
        # Plot individual samples
        for sampleIndex, sample in enumerate(samples):
            # Get values within this region
            "x can be reused from the bulk plot"
            sampleEdNCLS = depthNCLSDict["bulk1"][sample] \
                           if sample in depthNCLSDict["bulk1"] \
                           else depthNCLSDict["bulk2"][sample]
            y = [ stat for _, _, stat in sampleEdNCLS.find_overlap(contigID, start, end) ]
            
            # Get the line colour and type
            lineColour, lineType = SAMPLE_AESTHETICS[sampleIndex]
            
            # Plot the line
            axs[rowNum, colNum].plot(x, y, color=lineColour,
                                    linestyle=lineType,
                                    linewidth=linewidth)
        
        # Set xlim
        axs[rowNum, colNum].set_xlim(start, end)
        
        # Turn off ytick labels if not the first column
        if colNum > 0:
            axs[rowNum, colNum].set_yticklabels([])
        
        # Set up scale bar if this is the last row
        if plotScalebar == True:
            scalebar(axs, rowNum, colNum, start, end)
        # Otherwise, turn off x labels
        else:
            axs[rowNum, colNum].set_xticklabels([])
            axs[rowNum, colNum].locator_params(axis='x', nbins=4) # use less ticks; avoid clutter
    
    # Set y limits
    for colNum in range(len(regions)):
        axs[rowNum, colNum].set_ylim(np.floor(minY), np.ceil(maxY))

def scalebar(axs, rowNum, colNum, start, end):
    '''
    Start and end values are unused now but may be used in the future
    again.
    
    Parameters:
        axs -- a list of matplotlib.pyplot Axes objects to plot to
        rowNum -- an integer value indicating the row index to plot to
        colNum -- an integer value indicating the column index to plot to
        start -- [UNUSED]; an integer value indicating the start of the region
        end -- [UNUSED]; an integer value indicating the end of the region
    '''
    # Identify scale
    # if end / 1e6 >= 1:
    #     #scale = 1e6
    #     scaleLabel = "Mbp"
    # elif end / 1e3 >= 1:
    #     #scale = 1e3
    #     scaleLabel = "Kbp"
    # else:
    #     #scale = 1
    #     scaleLabel = "bp"
    
    # Set up x axis
    axs[rowNum, colNum].locator_params(axis='x', nbins=4) # use less ticks; avoid clutter
    axs[rowNum, colNum].set_xlabel("Chromosome position")
