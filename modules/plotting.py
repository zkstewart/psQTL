import os, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.patches import Rectangle

from .gff3 import GFF3

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
                power, width, height, outputDirectory,
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
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory to save the output TSV data
        linewidth -- OPTIONAL; an integer value indicating the width of the line plot (default=1)
        dotsize -- OPTIONAL; an integer value indicating the size of the dots (default=3)
    '''
    if line == False and scatter == False:
        raise ValueError("At least one of 'line' or 'scatter' must be True")
    
    # Get the maximum Y value across all regions
    maxY = 0
    for contigID, start, end in regions:
        regionValues = edNCLS.find_overlap(contigID, start, end)
        y = [stat for _, _, stat in regionValues]
        maxY = max(maxY, max(y))
    
    # Plot each region
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
        
        # Set ylim
        axs[rowNum, colNum].set_ylim(0, maxY + 0.1)
        
        # Turn off ytick labels if not the first column
        if colNum > 0:
            axs[rowNum, colNum].set_yticklabels([])
        
        # Plot scatter (if applicable)
        if scatter == True:
            axs[rowNum, colNum].scatter(x, y, color="red", s=dotsize, alpha=0.5, zorder=0)
        
        # Plot line (if applicable)
        if line == True:
            axs[rowNum, colNum].plot(x, smoothedY, zorder=1, linewidth=linewidth)
        
        # Derive our output file name for TSV data
        fileOut = os.path.join(outputDirectory, f"{contigID}.{start}-{end}.line.tsv")
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

def histogram(axs, rowNum, edNCLS, regions, binSize, binThreshold, power, width, height, outputDirectory):
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
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory to save the output TSV data
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
        x = np.arange(len(y))
        
        # Set ylim
        axs[rowNum, colNum].set_ylim(0, maxY + 5)
        
        # Turn off ytick labels if not the first column
        if colNum > 0:
            axs[rowNum, colNum].set_yticklabels([])
        
        # Plot bars
        axs[rowNum, colNum].bar(x, y, zorder=0)
        
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

def genes(axs, rowNum, edNCLS, regions, gff3Obj, power, width, height):
    '''
    Parameters:
        axs -- a list of matplotlib.pyplot Axes objects to plot to
        rowNum -- an integer value indicating the row index to plot to
        edNCLS -- an EDNCLS object
        regions -- a list of lists containing three values: [contigID, start, end]
        gff3Obj -- a GFF3 class object from gff3.py in this repository
        power -- an integer value indicating what power statistical values were raised to
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
    Returns:
        genePltList -- a list of matplotlib.pyplot objects containing the genes plot data
                       per region
    '''
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
        
        # Resolve overlaps into separate lanes
        lanes = []
        for thisMrnaFeature in mrnaFeatures:
            if lanes == []:
                lanes.append([thisMrnaFeature])
                continue
            
            placed = False
            for lane in lanes:
                if (thisMrnaFeature.start - ARROW_SIZE) > lane[-1].end + ARROW_SIZE: # check for possible overlaps either side
                    lane.append(thisMrnaFeature)
                    placed = True
                    break
            if not placed:
                lanes.append([thisMrnaFeature])
        
        # Set limits
        axs[rowNum, colNum].set_xlim(start, end)
        axs[rowNum, colNum].set_ylim(0, len(lanes)+SPACING)
        
        # Turn off y labels
        axs[rowNum, colNum].set_yticklabels([])
        axs[rowNum, colNum].set_yticks([])
        
        # Plot each lane
        for laneNum, lane in enumerate(lanes):
            for mrnaFeature in lane:
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
