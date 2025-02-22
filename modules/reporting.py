import math
from .gff3 import GFF3

def report_genes(edNCLS, gff3Obj, regions, outputFileName, radiusSize=50000):
    '''
    Receives an EDNCLS object, a GFF3 object, a list of regions, and an output file name.
    Writes a summary of the number of SNPs that occur proximally to genes in the region(s)
    and their average ED value to the output file.
    
    Parameters:
        edNCLS -- an EDNCLS object
        gff3Obj -- a GFF3 class object from gff3.py in this repository
        regions -- a list of lists containing four values: [contigID, start, end, reverse]
        radiusSize -- OPTIONAL; an integer indicating the radius surrounding a gene that
                      you want to consider as being 'proximal' to a gene; default
                      == 50000 (bp).
    '''
    delimiter = "\t" if outputFileName.endswith(".tsv") else ","
    
    with open(outputFileName, "w") as fileOut:
        # Write header
        fileOut.write(delimiter.join(["gene_id", "contig", "strand", "coords", "left_snps", "left_mean",
                                      "intron_snps", "intron_mean", "utr_snps", "utr_mean",
                                      "cds_snps", "cds_mean", "right_snps", "right_mean"]) + "\n")
        
        # Iterate over each region
        for contigID, start, end, reverse in regions:
            # Get all SNPs that occur within this region
            posList, edList = [], []
            for pos, _, ed in edNCLS.find_overlap(contigID, start, end):
                posList.append(pos)
                edList.append(ed)
            
            # Get longest isoform for each gene in this region
            geneFeatures = gff3Obj.ncls_finder(start, end, "contig", contigID)
            mrnaFeatures = [
                GFF3.longest_isoform(geneFeature)
                for geneFeature in geneFeatures
                if hasattr(geneFeature, "mRNA")
            ]
            
            # Associate each SNP to their nearest gene
            mrnaSnps = { mrnaFeature.ID: [] for mrnaFeature in mrnaFeatures }
            for pos, ed in zip(posList, edList):
                # Handle SNPs that are within a gene
                within = [ mrnaFeature.ID for mrnaFeature in mrnaFeatures if mrnaFeature.start <= pos <= mrnaFeature.end ]
                if within != []:
                    for mrnaID in within:
                        mrnaSnps[mrnaID].append((pos, ed))
                # Handle intergenic SNPs
                else:
                    # Find most proximal gene
                    adjacent = [
                        (mrnaFeature.ID, min(abs(pos - mrnaFeature.start), abs(pos - mrnaFeature.end)))
                        for mrnaFeature in mrnaFeatures
                    ]
                    adjacent.sort(key=lambda x: x[1])
                    
                    # Store SNP in nearest gene (allowing ties)
                    if adjacent != []:
                        adjacent = [ (geneID, distance) for geneID, distance in adjacent if distance == adjacent[0][1] ]
                        for geneID, distance in adjacent:
                            if distance <= radiusSize: # filter out SNPs that are too far away
                                mrnaSnps[geneID].append((pos, ed))
            
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
                left, utr, intron, cds, right = [], [], [], [], []
                for pos, ed in mrnaSnpsList:
                    # Intergenic
                    if pos < mrnaFeature.start:
                        left.append((pos, ed))
                    elif pos > mrnaFeature.end:
                        right.append((pos, ed))
                    # Within
                    elif any([start <= pos <= end for start, end in cdsCoords]):
                        cds.append((pos, ed))
                    elif any([start <= pos <= end for start, end in exonCoords]):
                        utr.append((pos, ed))
                    else:
                        intron.append((pos, ed))
                
                # Calculate mean ED for each category
                left_mean = sum([ed for _, ed in left]) / len(left) if left != [] else 0
                utr_mean = sum([ed for _, ed in utr]) / len(utr) if utr != [] else 0
                intron_mean = sum([ed for _, ed in intron]) / len(intron) if intron != [] else 0
                cds_mean = sum([ed for _, ed in cds]) / len(cds) if cds != [] else 0
                right_mean = sum([ed for _, ed in right]) / len(right) if right != [] else 0
                
                # Write to file
                fileOut.write(f"{mrnaFeature.ID}\t{contigID}\t{mrnaFeature.strand}\t{mrnaFeature.start}-{mrnaFeature.end}\t" +
                              f"{len(left)}\t{left_mean:.4f}\t{len(intron)}\t{intron_mean:.4f}\t" +
                              f"{len(utr)}\t{utr_mean:.4f}\t{len(cds)}\t{cds_mean:.4f}\t" +
                              f"{len(right)}\t{right_mean:.4f}\n")

def report_depth(edNCLS, gff3Obj, regions, outputFileName, radiusSize=50000):
    '''
    Receives an EDNCLS object, a GFF3 object, a list of regions, and an output file name.
    Writes a summary of the number of SNPs that occur proximally to genes in the region(s)
    and their average ED value to the output file.
    
    Parameters:
        edNCLS -- an EDNCLS object
        gff3Obj -- a GFF3 class object from gff3.py in this repository
        regions -- a list of lists containing four values: [contigID, start, end, reverse]
        radiusSize -- OPTIONAL; an integer indicating the radius surrounding a gene that
                      you want to consider as being 'proximal' to a gene; default
                      == 50000 (bp).
    '''
    delimiter = "\t" if outputFileName.endswith(".tsv") else ","
    
    with open(outputFileName, "w") as fileOut:
        # Write header
        fileOut.write(delimiter.join(["gene_id", "contig", "strand", "coords", "proximal_segregating_deletion",
                                      "deletion_location", "deletion_ed"]) + "\n")
        
        # Iterate over each region
        for contigID, start, end, reverse in regions:
            # Get all deletions that occur within this region
            posList, edList = [], []
            for startPos, endPos, ed in edNCLS.find_overlap(contigID, start, end):
                posList.append((startPos, endPos))
                edList.append(ed)
            
            # Get longest isoform for each gene in this region
            geneFeatures = gff3Obj.ncls_finder(start, end, "contig", contigID)
            mrnaFeatures = [
                GFF3.longest_isoform(geneFeature)
                for geneFeature in geneFeatures
                if hasattr(geneFeature, "mRNA")
            ]
            
            # Iterate over each gene to check their proximity to deletions
            for mrnaFeature in mrnaFeatures:
                # Find the most segregating deletion that occurs proximally to this gene
                closest = [None, None, -math.inf]
                for startPos, endPos, ed in edNCLS.find_overlap(contigID, mrnaFeature.start-radiusSize, mrnaFeature.end+radiusSize):
                    if ed > closest[2]:
                        closest = [startPos, endPos-1, ed] # -1 to offset NCLS range behaviour
                
                # Categorise deletion location
                if closest[0] <= mrnaFeature.start and closest[1] >= mrnaFeature.end:
                    location = "within"
                else:
                    location = "adjacent"
                
                # Write to file
                fileOut.write(f"{mrnaFeature.ID}\t{contigID}\t{mrnaFeature.strand}\t{mrnaFeature.start}-{mrnaFeature.end}\t" +
                              f"{closest[0]}-{closest[1]}\t{location}\t{closest[2]:.4f}\n")
