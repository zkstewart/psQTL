from .gff3 import GFF3Graph

def is_overlapping(start1, end1, start2, end2):
    '''
    See https://nedbatchelder.com/blog/201310/range_overlap_in_two_compares.html
    
    Parameters:
        start1 -- an integer indicating the start of the first range
        end1 -- an integer indicating the end of the first range
        start2 -- an integer indicating the start of the second range
        end2 -- an integer indicating the end of the second range
    '''
    return end1 >= start2 and end2 >= start1

def report_genes_call(windowedNCLS, gff3Obj, regions, outputFileName, radiusSize=50000):
    '''
    Receives a WindowedNCLS object, a GFF3 object, a list of regions, and an output file name.
    The values in the WindowedNCLS object are expected to be from a 'call' ED analysis with power
    transformation already applied.
    
    Calculates the N1 metric, which is the number of SNPs with ED >= 1.
    
    Writes a summary of the number of SNPs that occur proximally to genes in the region(s)
    and their average ED value to the output file.
    
    Parameters:
        windowedNCLS -- a WindowedNCLS object
        gff3Obj -- a GFF3 class object from gff3.py in this repository
        regions -- a list of dictionaries, each dict structured like:
                       {
                            "contig": contigID, # string
                            "start": start, # int
                            "end": end, # int
                            "reverse": reverse] # bool
                        }
        radiusSize -- OPTIONAL; an integer indicating the radius surrounding a gene that
                      you want to consider as being 'proximal' to a gene; default
                      == 50000 (bp).
    '''
    delimiter = "\t" if outputFileName.endswith(".tsv") else ","
    
    with open(outputFileName, "w") as fileOut:
        # Write header
        fileOut.write(delimiter.join(["gene_id", "contig", "strand", "gene_start", "gene_end",
                                      "N1", "max_ED_within_radius", "num_exon", "exon_mean",
                                      "num_intron", "intron_mean", "num_adjacent", "adjacent_mean"]) + "\n")
        
        # Iterate over each region
        for regionDict in regions:
            # Decompose dictionary into variables
            contigID = regionDict["contig"]
            start = regionDict["start"]
            end = regionDict["end"]
            reverse = regionDict["reverse"]
            
            # Get longest isoform for each gene in this region
            geneFeatures = gff3Obj.ncls_finder(start, end, "contig", contigID)
            mrnaFeatures = [
                GFF3Graph.longest_isoform(geneFeature)
                for geneFeature in geneFeatures
                if hasattr(geneFeature, "mRNA")
            ]
            
            # Associate each gene with their nearest SNPs
            mrnaSnps = { mrnaFeature.ID: [] for mrnaFeature in mrnaFeatures }
            for mrnaFeature in mrnaFeatures:
                # Query for SNPs in proximity to this gene
                markers = []
                for pos, _, ed in windowedNCLS.find_overlap(contigID, mrnaFeature.start-radiusSize, mrnaFeature.end+radiusSize):
                    markers.append((pos, ed))
                
                # Separate SNPs into those that are within the gene and those that are adjacent
                for pos, ed in markers:
                    if mrnaFeature.start <= pos <= mrnaFeature.end:
                        mrnaSnps[mrnaFeature.ID].append((pos, ed, "within"))
                    else:
                        mrnaSnps[mrnaFeature.ID].append((pos, ed, "adjacent"))
            
            # Iterate over each gene in this region
            for mrnaFeature in mrnaFeatures:
                mrnaSnpsList = mrnaSnps[mrnaFeature.ID]
                
                # Skip if no SNPs are found
                if mrnaSnpsList == []:
                    continue
                
                # Skip if this mRNA has no exons
                if not hasattr(mrnaFeature, "exon"):
                    continue
                
                # Get CDS and exon coordinates
                if hasattr(mrnaFeature, "CDS"):
                    cdsCoords = [
                        (cdsFeature.start, cdsFeature.end)
                        for cdsFeature in mrnaFeature.CDS
                    ]
                else:
                    cdsCoords = []
                
                exonCoords = [
                    (exonFeature.start, exonFeature.end)
                    for exonFeature in mrnaFeature.exon
                ]
                
                # Categorise SNP locations
                adjacent, exon, intron = [], [], []
                for pos, ed, location in mrnaSnpsList:
                    # Intergenic
                    if location == "adjacent":
                        adjacent.append((pos, ed))
                    # Within
                    elif any([start <= pos <= end for start, end in exonCoords]):
                        exon.append((pos, ed))
                    else:
                        intron.append((pos, ed))
                
                # Calculate mean ED for each category
                adjacent_mean = sum([ed for _, ed in adjacent]) / len(adjacent) if adjacent != [] else -1
                exon_mean = sum([ed for _, ed in exon]) / len(exon) if exon != [] else -1
                intron_mean = sum([ed for _, ed in intron]) / len(intron) if intron != [] else -1
                
                # Calculate N1 and max ED within radius
                N1 = sum([ 1 for _, ed, _ in mrnaSnpsList if ed >= 1 ])
                max_ED_within_radius = max([ ed for _, ed, _ in mrnaSnpsList ])
                
                # Write to file
                fileOut.write(f"{mrnaFeature.ID}\t{contigID}\t{mrnaFeature.strand}\t{mrnaFeature.start}\t{mrnaFeature.end}\t" +
                              f"{N1}\t{max_ED_within_radius:.4f}\t" +
                              f"{len(exon)}\t{exon_mean:.4f}\t{len(intron)}\t{intron_mean:.4f}\t" +
                              f"{len(adjacent)}\t{adjacent_mean:.4f}\n")

def report_genes_depth(windowedNCLS, gff3Obj, regions, outputFileName, radiusSize=50000):
    '''
    Receives a WindowedNCLS object, a GFF3 object, a list of regions, and an output file name.
    The values in the WindowedNCLS object are expected to be from a 'depth' ED analysis with power
    transformation already applied.
    
    Calculates the N1 metric, which is the number of SNPs with ED >= 1.
    
    Writes a summary of the number of SNPs that occur proximally to genes in the region(s)
    and their average ED value to the output file.
    
    Parameters:
        windowedNCLS -- a WindowedNCLS object
        gff3Obj -- a GFF3 class object from gff3.py in this repository
        regions -- a list of dictionaries, each dict structured like:
                       {
                            "contig": contigID, # string
                            "start": start, # int
                            "end": end, # int
                            "reverse": reverse] # bool
                        }
        radiusSize -- OPTIONAL; an integer indicating the radius surrounding a gene that
                      you want to consider as being 'proximal' to a gene; default
                      == 50000 (bp).
    '''
    delimiter = "\t" if outputFileName.endswith(".tsv") else ","
    
    with open(outputFileName, "w") as fileOut:
        # Write header
        fileOut.write(delimiter.join(["gene_id", "contig", "strand", "gene_start", "gene_end",
                                      "N1", "max_ED_within_radius", "num_overlapping", "overlapping_mean",
                                      "num_adjacent", "adjacent_mean"]) + "\n")
        
        # Iterate over each region
        for regionDict in regions:
            # Decompose dictionary into variables
            contigID = regionDict["contig"]
            start = regionDict["start"]
            end = regionDict["end"]
            reverse = regionDict["reverse"]
            
            # Get longest isoform for each gene in this region
            geneFeatures = gff3Obj.ncls_finder(start, end, "contig", contigID)
            mrnaFeatures = [
                GFF3Graph.longest_isoform(geneFeature)
                for geneFeature in geneFeatures
                if hasattr(geneFeature, "mRNA")
            ]
            
            # Associate each gene with their nearest CNVs
            mrnaSnps = { mrnaFeature.ID: [] for mrnaFeature in mrnaFeatures }
            for mrnaFeature in mrnaFeatures:
                # Query for CNVs in proximity to this gene
                markers = []
                for pos, windowEnd, ed in windowedNCLS.find_overlap(contigID, mrnaFeature.start-radiusSize, mrnaFeature.end+radiusSize):
                    markers.append((pos, windowEnd-1, ed)) # windowEnd gives pos+windowSize, then -1 to account for inclusive ranges
                
                # Separate SNPs into those that are within the gene and those that are adjacent
                for pos, windowEnd, ed in markers:
                    if is_overlapping(mrnaFeature.start, mrnaFeature.end, pos, windowEnd):
                        mrnaSnps[mrnaFeature.ID].append((pos, ed, "overlapping"))
                    else:
                        mrnaSnps[mrnaFeature.ID].append((pos, ed, "adjacent"))
            
            # Iterate over each gene in this region
            for mrnaFeature in mrnaFeatures:
                mrnaSnpsList = mrnaSnps[mrnaFeature.ID]
                
                # Skip if no CNVs are found
                if mrnaSnpsList == []:
                    continue
                
                # Skip if this mRNA has no exons
                if not hasattr(mrnaFeature, "exon"):
                    continue
                
                # Get CDS and exon coordinates
                if hasattr(mrnaFeature, "CDS"):
                    cdsCoords = [
                        (cdsFeature.start, cdsFeature.end)
                        for cdsFeature in mrnaFeature.CDS
                    ]
                else:
                    cdsCoords = []
                
                exonCoords = [
                    (exonFeature.start, exonFeature.end)
                    for exonFeature in mrnaFeature.exon
                ]
                
                # Categorise CNV locations
                adjacent, overlapping = [], []
                for pos, ed, location in mrnaSnpsList:
                    # Intergenic
                    if location == "adjacent":
                        adjacent.append((pos, ed))
                    # Within
                    else:
                        overlapping.append((pos, ed))
                
                # Calculate mean ED for each category
                adjacent_mean = sum([ed for _, ed in adjacent]) / len(adjacent) if adjacent != [] else -1
                overlapping_mean = sum([ed for _, ed in overlapping]) / len(overlapping) if overlapping != [] else -1
                
                # Calculate N1 and max ED within radius
                N1 = sum([ 1 for _, ed, _ in mrnaSnpsList if ed >= 1 ])
                max_ED_within_radius = max([ ed for _, ed, _ in mrnaSnpsList ])
                
                # Write to file
                fileOut.write(f"{mrnaFeature.ID}\t{contigID}\t{mrnaFeature.strand}\t{mrnaFeature.start}\t{mrnaFeature.end}\t" +
                              f"{N1}\t{max_ED_within_radius:.4f}\t" +
                              f"{len(overlapping)}\t{overlapping_mean:.4f}\t" +
                              f"{len(adjacent)}\t{adjacent_mean:.4f}\n")

def report_genes_splsda(windowedNCLSDict, gff3Obj, regions, outputFileName, radiusSize=50000):
    '''
    Receives a windowedNCLSDict object, a GFF3 object, a list of regions, and an output file name.
    The values in the windowedNCLSDict object are expected to be from one or more different
    sPLS-DA analyses (e.g., 'call', 'depth', 'integrated_call', 'integrated_depth').
    
    Writes a summary of the number of features that occur proximally to genes in the region(s)
    and their average ED value to the output file.
    
    Parameters:
        windowedNCLSDict -- one to three WindowedNCLS objects in a dictionary with keys:
            call -- the selected SNPs from sPLS-DA analysis of 'call' variants
            depth -- the selected SNPs from sPLS-DA analysis of 'depth' variants
            integrated_call & integrated_depth -- the selected SNPs from sPLS-DA analysis of 'integrated' variants
        gff3Obj -- a GFF3 class object from gff3.py in this repository
        regions -- a list of dictionaries, each dict structured like:
                       {
                            "contig": contigID, # string
                            "start": start, # int
                            "end": end, # int
                            "reverse": reverse] # bool
                        }
        radiusSize -- OPTIONAL; an integer indicating the radius surrounding a gene that
                      you want to consider as being 'proximal' to a gene; default
                      == 50000 (bp).
    '''
    delimiter = "\t" if outputFileName.endswith(".tsv") else ","
    
    with open(outputFileName, "w") as fileOut:
        # Write header
        fileOut.write(delimiter.join(["gene_id", "contig", "strand", "gene_start", "gene_end",
                                      "num_adjacent_call_selected", "num_overlapping_call_selected",
                                      "num_adjacent_depth_selected", "num_overlapping_depth_selected",
                                      "num_adjacent_integrated_selected", "num_overlapping_integrated_selected"]) + "\n")
        
        # Iterate over each region
        for regionDict in regions:
            # Decompose dictionary into variables
            contigID = regionDict["contig"]
            start = regionDict["start"]
            end = regionDict["end"]
            reverse = regionDict["reverse"]
            
            # Get longest isoform for each gene in this region
            geneFeatures = gff3Obj.ncls_finder(start, end, "contig", contigID)
            mrnaFeatures = [
                GFF3Graph.longest_isoform(geneFeature)
                for geneFeature in geneFeatures
                if hasattr(geneFeature, "mRNA")
            ]
            
            # Associate each gene with their nearest markers
            mrnaSnps = { mrnaFeature.ID: [] for mrnaFeature in mrnaFeatures }
            for mrnaFeature in mrnaFeatures:
                # Query for markers in proximity to this gene
                markers = []
                for datasetKey, windowedNCLS in windowedNCLSDict.items():
                    for pos, windowEnd, splsdaStat in windowedNCLS.find_overlap(contigID, mrnaFeature.start-radiusSize, mrnaFeature.end+radiusSize): # splsdaStat is ignored
                        if "integrated" in datasetKey:
                            markers.append((pos, windowEnd-1, "integrated")) # windowEnd gives pos+windowSize, then -1 to account for inclusive ranges
                        else:
                            markers.append((pos, windowEnd-1, datasetKey))
                
                # Separate markers into those that are within the gene and those that are adjacent
                for pos, windowEnd, datasetKey in markers:
                    if is_overlapping(mrnaFeature.start, mrnaFeature.end, pos, windowEnd):
                        mrnaSnps[mrnaFeature.ID].append((pos, datasetKey, "overlapping"))
                    else:
                        mrnaSnps[mrnaFeature.ID].append((pos, datasetKey, "adjacent"))
            
            # Iterate over each gene in this region
            for mrnaFeature in mrnaFeatures:
                mrnaSnpsList = mrnaSnps[mrnaFeature.ID]
                
                # Skip if no features are found
                if mrnaSnpsList == []:
                    continue
                
                # Skip if this mRNA has no exons
                if not hasattr(mrnaFeature, "exon"):
                    continue
                
                # Get CDS and exon coordinates
                if hasattr(mrnaFeature, "CDS"):
                    cdsCoords = [
                        (cdsFeature.start, cdsFeature.end)
                        for cdsFeature in mrnaFeature.CDS
                    ]
                else:
                    cdsCoords = []
                
                exonCoords = [
                    (exonFeature.start, exonFeature.end)
                    for exonFeature in mrnaFeature.exon
                ]
                
                # Categorise feature locations
                adjacent_call, adjacent_depth, adjacent_integrated = [], [], []
                overlapping_call, overlapping_depth, overlapping_integrated = [], [], []
                for pos, datasetKey, location in mrnaSnpsList:
                    # Intergenic
                    if location == "adjacent":
                        if datasetKey == "call":
                            adjacent_call.append(pos)
                        elif datasetKey == "depth":
                            adjacent_depth.append(pos)
                        elif datasetKey == "integrated":
                            adjacent_integrated.append(pos)
                        else:
                            raise ValueError(f"Unhandled dataset key: {datasetKey}")
                    # Overlapping
                    else:
                        if datasetKey == "call":
                            overlapping_call.append(pos)
                        elif datasetKey == "depth":
                            overlapping_depth.append(pos)
                        elif datasetKey == "integrated":
                            overlapping_integrated.append(pos)
                        else:
                            raise ValueError(f"Unhandled dataset key: {datasetKey}")
                
                # Write to file
                fileOut.write(f"{mrnaFeature.ID}\t{contigID}\t{mrnaFeature.strand}\t{mrnaFeature.start}\t{mrnaFeature.end}\t" +
                              f"{len(adjacent_call)}\t{len(overlapping_call)}\t" +
                              f"{len(adjacent_depth)}\t{len(overlapping_depth)}\t" +
                              f"{len(adjacent_integrated)}\t{len(overlapping_integrated)}\n")

def extract_markers(windowedNCLSObj, contigID, start, end):
    '''
    Helper function to allow report_markers() to extract markers from a WindowedNCLS object
    (e.g., from a 'call' ED analysis) or a dictionary of WindowedNCLS objects as expected
    of sPLS-DA results.
    
    Parameters:
        windowedNCLSObj -- a WindowedNCLS object or a dictionary of WindowedNCLS objects
    Returns:
        markers -- a list of tuples, each containing (position, windowEnd, statistic)
                   where statistic can be the ED value, or the type of feature that was
                   selected
    '''
    # Grab all markers that overlap with the specified region
    markers = {}
    if isinstance(windowedNCLSObj, dict):
        for datasetKey, windowedNCLS in windowedNCLSObj.items():
            for pos, windowEnd, splsdaStat in windowedNCLS.find_overlap(contigID, start, end): # splsdaStat is ignored
                markers.setdefault(pos, [])
                markers[pos].append([windowEnd, datasetKey]) # windowEnd gives pos+windowSize
    else:
        for pos, windowEnd, statistic in windowedNCLSObj.find_overlap(contigID, start, end):
            markers.setdefault(pos, [])
            markers[pos].append([windowEnd, statistic]) # windowEnd gives pos+windowSize
    
    # Remove duplicates
    deduplicated = []
    for pos in markers.keys():
        markersList = markers[pos]
        if len(markersList) > 1:
            # Handle sPLS-DA statistics
            if isinstance(markersList[0][1], str): # sPLS-DA statistics are strings
                markersList.sort(key=lambda x: (0 if "integrated" in x else 1, # sort by integrated > selected
                                                0 if "depth" in x else 1)) # and then by depth > call
            # Handle ED statistics
            else: # ED statistics are floats
                markersList.sort(key=lambda x: x[1], reverse=True) # sort by ED value descending
        deduplicated.append([pos, *markersList[0]]) # reconstitute a tuple format (pos, windowEnd, statistic)
    
    # Convert to dictionary format
    markers = {
        pos: {"statistic": statistic, "windowEnd": windowEnd, "genes": []}
        for pos, windowEnd, statistic in deduplicated
    }
    return markers

def report_markers(windowedNCLS, gff3Obj, regions, outputFileName, radiusSize=50000):
    '''
    Receives a WindowedNCLS object, a GFF3 object, a list of regions, and an output file name.
    The values in the WindowedNCLS object are expected to be from a 'call' ED analysis with power
    transformation already applied.
    
    Calculates the N1 metric, which is the number of SNPs with ED >= 1.
    
    Writes a summary of the number of SNPs that occur proximally to genes in the region(s)
    and their average ED value to the output file.
    
    Parameters:
        windowedNCLS -- a WindowedNCLS object, containing any type of results (call, depth, etc.)
        gff3Obj -- a GFF3 class object from gff3.py in this repository
        regions -- a list of dictionaries, each dict structured like:
                       {
                            "contig": contigID, # string
                            "start": start, # int
                            "end": end, # int
                            "reverse": reverse] # bool
                        }
        radiusSize -- OPTIONAL; an integer indicating the radius surrounding a gene that
                      you want to consider as being 'proximal' to a gene; default
                      == 50000 (bp).
    '''
    delimiter = "\t" if outputFileName.endswith(".tsv") else ","
    
    with open(outputFileName, "w") as fileOut:
        # Write header
        fileOut.write(delimiter.join(["contig", "position", "statistic", "nearest_gene"]) + "\n")
        
        # Iterate over each region
        for regionDict in regions:
            # Decompose dictionary into variables
            contigID = regionDict["contig"]
            start = regionDict["start"]
            end = regionDict["end"]
            reverse = regionDict["reverse"]
            
            # Get all markers that occur within this region
            markers = extract_markers(windowedNCLS, contigID, start, end)
            
            # Associate each marker with its nearest gene
            for pos, posDict in markers.items():
                # Get longest isoform for each gene in this region
                "Minus 1 from windowEnd to account for how WindowedNCLS adds 1 to the end position to allow for inclusive ranges"
                geneFeatures = gff3Obj.ncls_finder(pos-radiusSize, posDict["windowEnd"]-1+radiusSize, "contig", contigID)
                mrnaFeatures = [
                    GFF3Graph.longest_isoform(geneFeature)
                    for geneFeature in geneFeatures
                    if hasattr(geneFeature, "mRNA")
                ]
                
                # Locate genes that overlap with the marker position
                within = [
                    mrnaFeature.ID
                    for mrnaFeature in mrnaFeatures
                    if is_overlapping(mrnaFeature.start, mrnaFeature.end, pos, posDict["windowEnd"]-1)
                ]
                if within != []:
                    for mrnaID in within:
                        posDict["genes"].append((mrnaID, "within"))
                # Handle intergenic markers
                else:
                    # Find most proximal gene
                    adjacent = [
                        (mrnaFeature.ID, min(
                            abs(pos - mrnaFeature.start),
                            abs(pos - mrnaFeature.end),
                            abs(posDict["windowEnd"] - mrnaFeature.start),
                            abs(posDict["windowEnd"] - mrnaFeature.end)
                        ))
                        for mrnaFeature in mrnaFeatures
                    ]
                    adjacent.sort(key=lambda x: x[1])
                    
                    # Associate gene to this marker (allowing ties)
                    if adjacent != []:
                        adjacent = [ (mrnaID, distance) for mrnaID, distance in adjacent if distance == adjacent[0][1] ]
                        for mrnaID, distance in adjacent:
                            if distance <= radiusSize: # filter out SNPs that are too far away
                                posDict["genes"].append((mrnaID, "adjacent"))
            
            # Iterate over each marker in this region
            for pos, posDict in markers.items():
                statistic = f"{posDict['statistic']:.4f}" if isinstance(posDict['statistic'], float) \
                            else posDict['statistic'] # sPLS-DA statistics are strings
                genes = '; '.join([mrnaID for mrnaID, location in posDict['genes']]) if posDict['genes'] != [] else "."
                # Write to file
                fileOut.write(f"{contigID}\t{pos}\t{statistic}\t{genes}\n")
