import os, gzip

from .parsing import read_gz_file

def recode_vcf(vcfFile, outputFileName):
    '''
    Recode a VCF file to a format suitable for sPLS-DA analysis. Genotypes are encoded as
    an integer (0, 1, or 2) based on the number of minor alleles present. The output file
    is a TSV file with the following format:
    [chrom, pos, sample1EncodedGT, sample2EncodedGT, ...]
    
    Parameters:
        vcfFile -- a string indicating the location of the input VCF file; can be gzipped
        outputFileName -- a string indicating the location of the output file; will be gzipped
    '''
    with read_gz_file(vcfFile) as fileIn, gzip.open(outputFileName, "wt") as fileOut:
        for line in fileIn:
            sl = line.strip().split("\t")
            
            # Handle #CHROM line
            if line.startswith("#CHROM"):
                fileOut.write("\t".join(["chrom", "pos"] + sl[9:]) + "\n")
                continue
            
            # Skip comment lines
            if line.startswith("#"):
                continue
            
            # Skip multiallelic lines
            if "," in sl[4]:
                continue
            
            # Identify genotype position
            gtIndex = sl[8].split(":").index("GT")
            
            # Count the number of alleles to determine major/minor allele
            ref = 0
            alt = 0
            for sampleData in sl[9:]:
                gt = sampleData.split(":")[gtIndex]
                for allele in gt.split("/"):
                    if allele == "0":
                        ref += 1
                    elif allele == "1":
                        alt += 1
            minor = "0" if alt > ref else "1"
            
            # Encode genotypes
            encodedLine = [sl[0], sl[1]]
            for sampleData in sl[9:]:
                gt = sampleData.split(":")[gtIndex]
                encodedGT = 0
                if "." in gt:
                    encodedGT = "."
                else:
                    encodedGT = gt.count(minor)
                encodedLine.append(str(encodedGT))
            
            # Write to output file
            fileOut.write("\t".join(encodedLine) + "\n")

## TBD: Implement function to call a utilities R script to run PLS-DA in windows
## TBD: Implement a function to run sPLS-DA on selected features from windows (capable of 1 or more inputs)