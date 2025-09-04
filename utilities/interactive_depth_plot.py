#! python3
# interactive_depth_plot.py
# Generates an alternative style of plot focusing on coverage values that have
# been median-normalised for clearer depiction of potential CNVs. This will
# produce an interactive plot with plotly that may assist with investigating
# CNVs.

import os, argparse, json, sys
import plotly.graph_objects as go
import numpy as np
from typing import Union

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.parsing import parse_metadata
from modules.depth import parse_bins_as_dict, normalise_coverage_dict, convert_dict_to_depthncls, convert_depth_to_alleles

GROUP_UNIFYING = {
    "bulk1": "group1",
    "bulk 1": "group1",
    "b1": "group1",
    "g1": "group1",
    "group 1": "group1",
    "1": "group1",
    "group1": "group1",
    "bulk2": "group2",
    "bulk 2": "group2",
    "b2": "group2",
    "g2": "group2",
    "group 2": "group2",
    "2": "group2",
    "group2": "group2"
}
COLOURS = {
    "group1": "#89A9F9",
    "group1_median": "#4851F9",
    "group2": "#FFE987",
    "group2_median": "#FF7A32"
}

class DirectoryNotFoundError(Exception):
    pass

def validate_interactive_args(args):
    '''
    Sets:
        depthFileDict --  a dictionary with structure like:
                          {
                              "group1": [[sampleID, depthFile], ...],
                              "group2": [[sampleID, depthFile], ...]
                          }
    '''
    args.workingDirectory = os.path.abspath(args.workingDirectory)
    
    # Validate numeric arguments
    if args.windowSize < 0:
        raise ValueError("-w must be a positive integer")
    if args.ploidy:
        if args.ploidy < 1:
            raise ValueError("--ploidy must be a value of 1 or greater")
    
    # Infer metadata location
    metadataJson = os.path.join(args.workingDirectory, "metadata.json")
    if not os.path.isfile(metadataJson):
        raise FileNotFoundError(f"'{metadataJson}' not found as expected")
    
    # Parse the metadata
    with open(metadataJson, "r") as fileIn:
        jsonData = json.load(fileIn)
    metadataDict = { "group1": jsonData["group1"], "group2": jsonData["group2"] }
    
    # Infer depth file locations from workingDirectory
    depthDir = os.path.join(args.workingDirectory, "depth")
    if not os.path.isdir(depthDir):
        raise DirectoryNotFoundError(f"Expected to find '{depthDir}' based on value given to -d")
    
    args.depthFiles = [
        os.path.join(depthDir, dfile)
        for dfile in os.listdir(depthDir)
        if dfile.endswith(f".binned.{args.windowSize}.tsv")
    ]
    if len(args.depthFiles) == 0:
        raise ValueError(f"Found no files with '.binned.{args.windowSize}.tsv' suffix in '{args.depthDir}'")
    
    # Make sure all depth files have metadata
    metadataValues = [ sampleID for values in metadataDict.values() for sampleID in values ]
    args.depthFileDict = {}
    for depthFile in args.depthFiles:
        base = os.path.basename(depthFile)
        found = False
        for group, values in metadataDict.items():
            args.depthFileDict.setdefault(group, [])
            for v in values:
                if base.startswith(v):
                    args.depthFileDict[group].append([v, depthFile])
                    found = True
                    break
        if not found:
            raise KeyError(f"'{base}' does not have any associated metadata in '{metadataJson}'")
    
    # Validate output file location
    args.outputDirectory = os.path.abspath(args.outputDirectory)
    if not os.path.exists(args.outputDirectory):
        os.makedirs(args.outputDirectory, exist_ok=True)
        print(f"Created '{args.outputDirectory}' as part of argument validation")

def simplify_depth(xy: np.ndarray, tolerance: Union[int, float]) -> np.ndarray:
    '''
    Simplify depth data for plotting by removing data points
    that are redundant or minimally informative.
    
    Algorithm works by obtaining the first and last points, in addition to
    any points that mark a departure from the previously stored point by more
    than the specified tolerance.
    
    Parameters:
        xy -- a 2D numpy array where the first column is the x-axis position
              and the second column is the y-axis value
        tolerance -- an integer or float value indicating the maximum increase
                     or decrease in y-axis value before a change is detected
                     and the point is retained
    Returns:
        simplifiedXY -- a 2D numpy array of the input data with some points removed
                        where deemed to be redundant or minimally informative
    '''
    def process_plateau(simplifiedXY, xPlateau, yPlateau):
        # Process plateaus with only one point
        if len(xPlateau) == 1:
            simplifiedXY.append([xPlateau[0], yPlateau[0]])
        
        # Process plateaus with multiple points
        else:
            # Obtain and store the median value of the plateau
            medianY = np.median(yPlateau)
            
            # Store the plateau's median value
            simplifiedXY.append([xPlateau[0], medianY])
            simplifiedXY.append([x-1, medianY])
    
    # Skip simplification if there are too few points
    if len(xy) < 3:
        return xy
    
    # Note first data points
    simplifiedXY = []
    xPlateau = [xy[0][0]]
    yPlateau = [xy[0][1]]
    
    # Begin iteration through remaining data points
    for i in range(1, len(xy)):
        x, y = xy[i]
        
        # Check for departure from this plateau's tolerance threshold
        if y > yPlateau[0] + tolerance or y < yPlateau[0] - tolerance:
            process_plateau(simplifiedXY, xPlateau, yPlateau)
            
            # Begin a new plateau at the point of departure
            xPlateau = [x]
            yPlateau = [y]
        
        # Otherwise, continue to build the plateau
        else:
            xPlateau.append(x)
            yPlateau.append(y)
    
    # Process the last plateau
    process_plateau(simplifiedXY, xPlateau, yPlateau)
    
    # Return the simplified data as a numpy array
    return np.array(simplifiedXY, dtype=float)

def round_to_pointfive(number):
    return round((number)*2)/2

def main():
    usage = """%(prog)s produces interactive plotly visualisation of median-normalised depths
    for each sample.
    """
    # Establish main parser
    p = argparse.ArgumentParser()
    
    # Set arguments shared by subparsers
    p.add_argument("-d", dest="workingDirectory",
                   required=True,
                   help="Specify the location where the analysis is being performed")
    p.add_argument("-w", dest="windowSize",
                   type=int,
                   required=True,
                   help="""Specify the window size that reads were
                   binned into for CNV calling""")
    p.add_argument("-c", dest="contigs",
                   required=True,
                   nargs="+",
                   help="Specify one or more contig IDs to produce plots for")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Specify the location to write the html output file(s)")
    p.add_argument("--ploidy", dest="ploidy",
                   required=False,
                   help="""Optionally, specify the ploidy level to convert median-normalised
                   depth to putative allele count""",
                   default=None)
    p.add_argument("--outputBed", dest="outputBedFileName",
                   required=False,
                   help="Optionally specify a location to write a BED file with the same data",
                   default=None)
    
    args = p.parse_args()
    validate_interactive_args(args) # sets args.metadataDict, args.depthFileDict
    
    # Parse depth data
    coverageDict = parse_bins_as_dict(args.depthFileDict)
    normalise_coverage_dict(coverageDict)
    depthNCLSDict = convert_dict_to_depthncls(coverageDict, args.windowSize)
    
    # Iterate through all contig(s) and make plot(s)
    for contigID in args.contigs:
        outputFileName = os.path.join(args.outputDirectory, f"{contigID}.coverages.html")
        if os.path.exists(outputFileName):
            print(f"Skipping '{outputFileName}' as it already exists...")
        
        # Produce this contig's figure
        fig = go.Figure()
        for group, sampleDepthNCLSDicts in depthNCLSDict.items():
            group = GROUP_UNIFYING[group]
            
            groupValues = []
            for sampleID, windowedNCLS in sampleDepthNCLSDicts.items():
                # Get X and Y axis data
                overlappingBins = list(windowedNCLS.find_overlap(contigID, 1, 999999999999999))
                x = [ int((start+end)/2) for start, end, _ in overlappingBins ] # pin value at middle of the window
                if args.ploidy:
                    y = [ stat for _, _, stat in overlappingBins ]
                    binDict = { x:y for x,y in zip(x, y)}
                    y = convert_depth_to_alleles(binDict)
                    groupValues.append(y)
                else:
                    y = [ round_to_pointfive(stat) for _, _, stat in overlappingBins ]
                    groupValues.append(y)
                
                # Simplify data to reduce plot size
                "Simplification attempts to reduce superfluous data points"
                xy = np.column_stack((x, y))
                xy = simplify_depth(xy, 0)
                simpleX, simpleY = xy[:, 0], xy[:, 1]
                
                fig.add_trace(go.Scatter(x=simpleX, y=simpleY,
                                         mode="lines",
                                         line=go.scatter.Line(color=COLOURS[group]),
                                         name=sampleID,
                                         opacity=0.25,
                                         legendgroup=group,
                                         showlegend=True))
            
            # Obtain group average
            groupValues = np.array(groupValues)
            if groupValues.size == 0: # skip if there are no values
                continue
            groupMedian = np.percentile(groupValues, 50, axis=0)
            
            # Simplify average values
            xy = np.column_stack((x, groupMedian))
            xy = simplify_depth(xy, 0)
            simpleX, simpleY = xy[:, 0], xy[:, 1]
            
            # Plot it
            fig.add_trace(go.Scatter(x=simpleX, y=simpleY,
                                     mode="lines",
                                     line=go.scatter.Line(color=COLOURS[f"{group}_median"],
                                                          dash="dot"),
                                     name=f"{group}_median",
                                     legendgroup=f"{group}_median",
                                     showlegend=True))
        
        # Add range slider
        fig.update_layout(
            ## Set title
            title_text = f"{contigID} median-normalised depth plot",
            ## Set x-axis slider
            xaxis=dict(
                rangeslider=dict(
                    visible=True
                ),
                type="linear"
            ),
            ## Set y-axis range flexibility
            yaxis=dict(
                autorange=True,
                fixedrange=False
            )
        )
        
        # Save output file
        fig.write_html(outputFileName)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
