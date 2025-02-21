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
        if abs(stat) >= binThreshold:
            binIndex = (pos-start) // binSize
            histo[binIndex] += 1
    return histo

def linescatter(axs, rowNum, edNCLS, regions, wmaSize, line, scatter, 
                outputDirectory, plotScalebar, fileSuffix,
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
        outputDirectory -- a string indicating the directory to save the output TSV data OR
                           None if no output is desired
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
                print(f"WARNING: region '{contigID, start, end}' has too few data points to apply WMA smoothing")
                smoothedY = y
        
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
        
        # Get the maximum Y value for this region
        if scatter == True:
            if y.size != 0:
                maxY = max(maxY, max(y))
        else:
            if smoothedY.size != 0:
                maxY = max(maxY, max(smoothedY))
        
        # Derive our output file name for TSV data
        if outputDirectory != None:
            fileOut = os.path.join(outputDirectory, f"{contigID}.{start}-{end}.{fileSuffix}_line.tsv")
            if os.path.isfile(fileOut):
                print(f"WARNING: Line/scatter plot data for region '{contigID, start, end}' already exists as '{fileOut}'; won't overwrite")
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
        if maxY == 0:
            axs[rowNum, colNum].set_ylim(0, 1)
        else:
            axs[rowNum, colNum].set_ylim(0, maxY + (maxY * YLIM_HEADSPACE))

def histogram(axs, rowNum, edNCLS, regions, binSize, binThreshold, outputDirectory, plotScalebar):
    '''
    Parameters:
        axs -- a list of matplotlib.pyplot Axes objects to plot to
        rowNum -- an integer value indicating the row index to plot to
        edNCLS -- an EDNCLS object
        regions -- a list of lists containing three values: [contigID, start, end]
        binSize -- an integer value indicating the size of the bins/windows
        binThreshold -- an integer value indicating the threshold for counting a variant
                        within a bin/window
        outputDirectory -- a string indicating the directory to save the output TSV data OR
                           None if no output is desired
        plotScalebar -- a boolean value indicating whether to plot a scalebar on the X axis
    '''
    # Get the maximum Y value across all regions
    maxY = 0
    for contigID, start, end in regions:
        y = bin_values(edNCLS.find_overlap(contigID, start, end),
                       start, end, binSize, binThreshold)
        if y.size != 0:
            maxY = max(maxY, max(y))
    
    # Plot each region
    for colNum, (contigID, start, end) in enumerate(regions):
        y = bin_values(edNCLS.find_overlap(contigID, start, end),
                       start, end, binSize, binThreshold)
        x = np.array([ (i * binSize) + start for i in range(len(y)) ])
        
        # Set limits
        axs[rowNum, colNum].set_xlim(start, end)
        if maxY == 0:
            axs[rowNum, colNum].set_ylim(0, 1)
        else:
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
        if outputDirectory != None:
            fileOut = os.path.join(outputDirectory, f"{contigID}.{start}-{end}.histo.tsv")
            if os.path.isfile(fileOut):
                print(f"WARNING: Histogram plot data for region '{contigID, start, end}' already exists as '{fileOut}'; won't overwrite")
                continue
            else:
                with open(fileOut, "w") as fileOutTSV:
                    fileOutTSV.write(f"contigID\twindow_start\twindow_end\tnum_variants >= {binThreshold}\n")
                    for xVal, yVal in zip(x, y):
                        windowEnd = xVal + binSize if xVal + binSize < end else end
                        fileOutTSV.write(f"{contigID}\t{xVal}\t{windowEnd-1}\t{yVal}\n")

def genes(fig, axs, rowNum, gff3Obj, regions, plotScalebar):
    '''
    Parameters:
        fig -- the matplotlib.pyplot Figure object that axes correspond to
        axs -- a list of matplotlib.pyplot Axes objects to plot to
        rowNum -- an integer value indicating the row index to plot to
        gff3Obj -- a GFF3 class object from gff3.py in this repository
        regions -- a list of lists containing three values: [contigID, start, end]
        plotScalebar -- a boolean value indicating whether to plot a scalebar on the X axis
    Returns:
        genePltList -- a list of matplotlib.pyplot objects containing the genes plot data
                       per region
    '''
    fig.canvas.draw() # need to draw the figure to get the renderer
    SPACING = 0.1
    alreadyWarned = False
    
    for colNum, (contigID, start, end) in enumerate(regions):        
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
            axs[rowNum, colNum].plot([max(mrnaFeature.start, start), min(mrnaFeature.end, end)], # truncated to region
                                     [laneNum + SPACING + (1-SPACING)/2]*2, # y values for start and end
                                     color="black", linewidth=1,
                                     zorder=0)
            
            # Plot gene directionality
            if mrnaFeature.strand in ["+", "-"]:
                if mrnaFeature.strand == "+":
                    lastExon = max(mrnaFeature.exon, key=lambda x: x.end)
                    arrowAnnot = axs[rowNum, colNum].annotate("▷", xy=(lastExon.end, laneNum + SPACING + (1-SPACING)/2), xycoords="data",
                                                 #xytext=(lastExon.end, laneNum + SPACING + (1-SPACING)/2), textcoords="data",
                                                 fontsize=8, color="grey", zorder=1,
                                                 ha="left", va="center")
                else:
                    lastExon = min(mrnaFeature.exon, key=lambda x: x.start)
                    arrowAnnot = axs[rowNum, colNum].annotate("◁", xy=(lastExon.start, laneNum + SPACING + (1-SPACING)/2), xycoords="data",
                                                 #xytext=(lastExon.start, laneNum + SPACING + (1-SPACING)/2), textcoords="data",
                                                 fontsize=8, color="grey", zorder=1,
                                                 ha="right", va="center")
            
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
            if hasattr(mrnaFeature, "CDS"):
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
            
            # Get the size of the arrow
            transf = axs[rowNum, colNum].transData.inverted()
            arrowBbox = arrowAnnot.get_window_extent(renderer = fig.canvas.renderer)
            arrowDataCoords = arrowBbox.transformed(transf)
            ARROW_SIZE = arrowDataCoords.x1 - arrowDataCoords.x0
            
            # Plot gene name
            geneName = mrnaFeature.ID
            textBox = axs[rowNum, colNum].text(
                mrnaFeature.end + ARROW_SIZE if (mrnaFeature.end + ARROW_SIZE) < end # x position
                    else end + ARROW_SIZE, # prevent text from going off the plot
                laneNum + SPACING + (1-SPACING)/2, # y position
                geneName, # text
                horizontalalignment="left", verticalalignment="center", # alignment
                fontsize=10,
                zorder=4 # above everything else
            )
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
    if len(lanes) == 0:
        axs[rowNum, colNum].set_ylim(0, 1)
    else:
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
        for bulk in ["bulk1", "bulk2"]:
            sampleDict = depthNCLSDict[bulk]
            
            # Get the average and upper/lower quantiles of each bulk's coverage values
            bulkValues = []
            for sample, edNCLS in sampleDict.items():
                # Skip any samples being individually plotted
                if sample in samples:
                    continue
                
                # Get values within this region
                overlappingBins = list(edNCLS.find_overlap(contigID, start, end))
                y = [ stat for _, _, stat in overlappingBins ]
                bulkValues.append(y)
                # Get the minimum Y value for this region
                if len(y) > 0:
                    minY = min(minY, min(y))
            bulkValues = np.array(bulkValues)
            
            # Skip if there are no values
            if bulkValues.size == 0:
                print(f"WARNING: No coverage values found for region '{contigID, start, end}' for '{bulk}'")
                continue
            
            # Get the median and Q1/Q3 quantiles
            q1, median, q3 = np.percentile(bulkValues, [25, 50, 75], axis=0)
            
            # Figure out the x values
            x = []
            for index, (windowStart, windowEnd, _) in enumerate(overlappingBins):
                x.append(windowStart + (windowEnd - windowStart)/2)
            
            # Extend tails for better visualisation
            x = np.concatenate(([start], x, [end]))
            median = np.concatenate(([median[0]], median, [median[-1]]))
            q1 = np.concatenate(([q1[0]], q1, [q1[-1]]))
            q3 = np.concatenate(([q3[0]], q3, [q3[-1]]))
            
            # Plot median and Q1/Q3 lines
            axs[rowNum, colNum].plot(x, median, linewidth=linewidth,
                                     color="palevioletred" if bulk == "bulk1" else "mediumseagreen")
            
            axs[rowNum, colNum].fill_between(x, q1, q3, alpha = 0.5,
                                     color="palevioletred" if bulk == "bulk1" else "mediumseagreen",
                                     label="_nolegend_")
            
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
            
            # Extend tails for better visualisation
            y = np.concatenate(([y[0]], y, [y[-1]]))
            
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
        if np.floor(minY) == np.ceil(maxY):
            axs[rowNum, colNum].set_ylim(0, 1)
        else:
            axs[rowNum, colNum].set_ylim(np.floor(minY), np.ceil(maxY))
    
    # Set legend
    legendLabels = ["bulk1", "bulk2"] + samples # samples can be []
    axs[rowNum, colNum].legend(legendLabels, # colNum is the last column
                          loc="center left",
                          bbox_to_anchor=(1, 0.5),
                          ncol=1)

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
