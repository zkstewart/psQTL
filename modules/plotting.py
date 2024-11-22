import os
import numpy as np
import matplotlib.pyplot as plt

def linescatter(edNCLS, regions, wmaSize, line, scatter, 
                power, width, height,
                linewidth=1, dotsize=3):
    '''
    Parameters:
        edNCLS -- an EDNCLS object
        regions -- a list of lists containing three values: [contigID, start, end]
        wmaSize -- an integer value indicating the number of previous values to consider
                   during weighted moving average calculation
        line -- a boolean value indicating whether to plot a line
        scatter -- a boolean value indicating whether to plot scatter points
        power -- an integer value indicating what power statistical values were raised to
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        linewidth -- OPTIONAL; an integer value indicating the width of the line plot (default=1)
        dotsize -- OPTIONAL; an integer value indicating the size of the dots (default=3)
    Returns:
        linePltList -- a list of matplotlib.pyplot objects containing the line and/or scatter
                       plot data per region
    '''
    if line == False and scatter == False:
        raise ValueError("At least one of 'line' or 'scatter' must be True")
    
    raise NotImplementedError("This function has not been implemented yet")

def histogram(edNCLS, regions, binSize, binThreshold, power, width, height):
    '''
    Parameters:
        edNCLS -- an EDNCLS object
        regions -- a list of lists containing three values: [contigID, start, end]
        binSize -- an integer value indicating the size of the bins/windows
        binThreshold -- an integer value indicating the threshold for counting a variant
                        within a bin/window
        power -- an integer value indicating what power statistical values were raised to
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
    Returns:
        histoPltList -- a list of matplotlib.pyplot objects containing the histogram plot data
                        per region
    '''
    raise NotImplementedError("This function has not been implemented yet")

def genes(edNCLS, regions, gff3Obj, power, width, height):
    '''
    Parameters:
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
    raise NotImplementedError("This function has not been implemented yet")
