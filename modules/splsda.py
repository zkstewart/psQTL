import shutil, subprocess, re
from collections import Counter

from .parsing import read_gz_file, WriteGzFile
from .ncls import RangeNCLS
from .ed import gt_median_adjustment
from .parsing import vcf_header_to_metadata_validation

def validate_r_exists():
    if not shutil.which("R"):
        raise FileNotFoundError("R not found in PATH")
    if not shutil.which("Rscript"):
        raise FileNotFoundError("Rscript not found in PATH")

def validate_r_package(packageName):
    '''
    Checks if the specified R package is installed via command-line
    input to Rscript.
    
    Parameters:
        packageName -- a string indicating the name of the R package to check
    Returns:
        isInstalled -- a boolean indicating whether the package is installed
    '''
    # Format command
    cmd = ["echo", f'\'find.package("{packageName}")\'', "|", "Rscript", "-"]
    
    # Check if package is installed through Rscript interface
    run_Rscript = subprocess.Popen(" ".join(cmd), shell=True,
                                   stdout = subprocess.PIPE,
                                   stderr = subprocess.PIPE)
    rout, rerr = run_Rscript.communicate()
    
    # Check for errors
    errorMsg = rerr.decode("utf-8")
    if "there is no package" in errorMsg or "error" in errorMsg:
        return False
    else:
        # If the package is installed, rout will contain the path to the package
        return True if rout.decode("utf-8").strip() else False

def validate_r_packages_installation():
    '''
    Sequentially checks if the required R packages for sPLS-DA are installed.
    Raises an exception if any package is not installed.
    '''
    REQUIRED_PACKAGES = ["argparser", "mixOmics"]
    for package in REQUIRED_PACKAGES:
        if not validate_r_package(package):
            raise FileNotFoundError(f"The R package '{package}' is not installed. "
                                    "Please install it before running sPLS-DA.")

def recode_variant(gtIndex, sampleFields):
    '''
    Encodes the genotype of a sample as the number of minor alleles present.
    
    Parameters:
        gtIndex -- an integer indicating the index of the GT field in the VCF format
        sampleFields -- a list of strings, each coming from the VCF sample data fields
    Returns:
        encodedGTs -- a list of strings, each representing the encoded genotype for a sample
    '''
    genotypes = [
        sampleData.split(":")[gtIndex].replace("|", "/").split("/")
        for sampleData in sampleFields
    ]
    
    # Count the number of alleles to determine major allele
    alleles = Counter([ allele for gt in genotypes for allele in gt if allele != "." ])
    mostCommonAlleles = alleles.most_common()
    majorAlleles = [
        allele
        for allele, count in mostCommonAlleles
        if count == mostCommonAlleles[0][1]
    ]
    
    # Pick a single major allele if there are multiple
    if len(majorAlleles) > 1:
        majorAlleles.sort(key = lambda x: int(x)) # prefer ref allele if possible
        majorAlleles = [majorAlleles[0]]
    
    # Encode genotype as the number of minor alleles
    encodedGTs = []
    for gt in genotypes:
        if "." in gt:
            encodedGT = "."
        else:
            encodedGT = sum([ 1 for allele in gt if not allele in majorAlleles ])
        encodedGTs.append(str(encodedGT))
    return encodedGTs

def recode_cnv(gtIndex, sampleFields):
    '''
    Encodes the CNV of a sample according to whether it has an allele copy number
    that is below or equal to the median (0) or above the median (1).
    
    Parameters:
        gtIndex -- an integer indicating the index of the GT field in the VCF format
        sampleFields -- a list of strings, each coming from the VCF sample data fields
    Returns:
        encodedGTs -- a list of strings, each representing the encoded genotype for a sample
    '''
    genotypes = [
        sampleData.split(":")[gtIndex].replace("|", "/").split("/")
        for sampleData in sampleFields
    ]
    
    # Median adjust genotypes
    alleles = [ list(map(int, gt)) if not "." in gt else "." for gt in genotypes ]
    alleles = gt_median_adjustment([alleles])[0] # unpack the result
    
    # Encode genotype as presence or absence of copy number above median
    encodedGTs = []
    for allele in alleles:
        if "." in allele:
            encodedGT = "."
        else:
            encodedGT = "2" if 2 in allele else "1" if 1 in allele else "0"
        encodedGTs.append(encodedGT)
    return encodedGTs

def recode_vcf(vcfFile, outputFileName, metadataDict, isCNV=False, quiet=False):
    '''
    Recode a VCF file to a format suitable for sPLS-DA analysis.
    
    For variant calls:
    Genotypes are encoded as an integer (0, 1, or 2) based on the number of minor
    alleles present.
    
    For CNVs:
    Allele copy numbers are encoded as 0 or 1, according to whether their copy
    number is below or equal to the median (0) or above the median (1).
    
    The output file is a TSV file with the following format:
    [chrom, start, end, sample1EncodedGT, sample2EncodedGT, ...]
    
    Start and end give the borders of the variant (1-based and inclusive).
    
    Parameters:
        vcfFile -- a string indicating the location of the input VCF file; can be gzipped
        outputFileName -- a string indicating the location of the output file; will be gzipped
        isCNV -- (OPTIONAL) a boolean indicating whether the genotypes are for CNVs
                 (True) or SNPs/indels (False); default is False
    '''
    windowSize = 1 # default value; will be reset if parsing a CNV file
    windowSizeRegex = re.compile(r"windowSize=(\d+)")
    
    # Parse through the file and write the recoded output
    with read_gz_file(vcfFile) as fileIn, WriteGzFile(outputFileName) as fileOut:
        for line in fileIn:
            sl = line.strip().split("\t")
            
            # Handle #CHROM line
            if line.startswith("#CHROM"):
                # Validate VCF header with respect to metadataDict and extract sample names
                samples = sl[9:] # This gives us the ordered sample IDs
                b1Samples, b2Samples = vcf_header_to_metadata_validation(samples, metadataDict, strict=False, quiet=quiet)
                foundSampleNames = b1Samples + b2Samples
                
                # Identify the indices of the samples in the VCF header
                sampleIndices = [ sl[9:].index(name) for name in foundSampleNames ]
                
                # Write the header line to the output file
                fileOut.write("\t".join(["chrom", "start", "end"] + [ sl[9:][i] for i in sampleIndices ]) + "\n")
            
            # Parse the psQTL_prep comment line
            if line.startswith("##psQTL_prep"):
                regexHit = windowSizeRegex.search(line)
                if regexHit is None:
                    raise ValueError(f"Expected VCF-like header line '{line.strip()}' to have a 'windowSize=(\\d+)' regex match")
                windowSize = int(regexHit.groups()[0])
                continue
            
            # Skip any other comment lines
            if line.startswith("#"):
                continue
            
            # Skip multiallelic lines
            if "," in sl[4]:
                continue
            
            # Identify genotype position
            gtIndex = sl[8].split(":").index("GT")
            
            # Recode the genotype according to variant type
            if isCNV:
                gtFields = recode_cnv(gtIndex, [ sl[9:][i] for i in sampleIndices ])
            else:
                gtFields = recode_variant(gtIndex, [ sl[9:][i] for i in sampleIndices ])
            
            # Format the output line
            encodedLine = [sl[0], sl[1], str(int(sl[1]) + windowSize - 1), *gtFields]
            
            # Write to output file
            fileOut.write("\t".join(encodedLine) + "\n")

def run_windowed_splsda(metadataFile, encodedVcfFile, outputVariants, outputBER, outputRdata,
                        scriptLocation, threads=1, windowSize=1000000, windowSizeIsSNPs=False,
                        berCutoff=0.4, maf=0.05, nrepeat=10, maxiters=1000):
    '''
    Calls the windowed_plsda.R script to run sPLS-DA on the provided encoded VCF file.
    
    Parameters:
        metadataFile -- a string indicating the location of the metadata file; should be a TSV file
                        with columns: [sampleID, groupNum] with NO header
        encodedVcfFile -- a string indicating the location of the encoded VCF file 
                          (i.e., a file generated by recode_vcf)
        outputVariants -- a string indicating the location of the output variants file
        outputBER -- a string indicating the location of the output BER file
        outputRdata -- a string indicating the location of the output RData file
        scriptLocation -- a string indicating the location of the windowed_plsda.R script
        threads -- (OPTIONAL) an integer indicating the number of threads to use (default is 1)
        windowSize -- (OPTIONAL) an integer indicating the size of the windows to run local
                       PLS-DA within (default is 1000000)
        windowSizeIsSNPs -- (OPTIONAL) a boolean indicating whether the window size is in
                            number of SNPs (True) or in base pairs (False; the default)
        berCutoff -- (OPTIONAL) a float indicating the BER cutoff to filter on (default is 0.4)
        maf -- (OPTIONAL) a float indicating the minor allele frequency threshold to
               filter on (default is 0.05)
        nrepeat -- (OPTIONAL) an integer indicating the number of repeats for stability analysis
        maxiters -- (OPTIONAL) an integer indicating the maximum number of iterations when tuning sPLS-DA
    '''
    # Format command
    cmd = ["Rscript", scriptLocation, metadataFile, encodedVcfFile,
           outputVariants, outputBER, outputRdata,
           "--threads", str(threads), "--windowSize", str(windowSize),
           "--berCutoff", str(berCutoff), "--MAF", str(maf),
            "--nrepeat", str(nrepeat), "--maxiters", str(maxiters)]
    if windowSizeIsSNPs:
        cmd.append("--windowSizeIsSNPs")
    
    # Run bcftools index
    run_Rscript = subprocess.Popen(" ".join(cmd), shell=True,
                                   stdout = subprocess.DEVNULL,
                                   stderr = subprocess.PIPE)
    rout, rerr = run_Rscript.communicate()
    
    # Check for errors
    if run_Rscript.returncode == 0:
        return None
    else:
        errorMsg = rerr.decode("utf-8").rstrip("\r\n ")
        if "Y variable must have exactly two unique values" in errorMsg:
            raise ValueError(f"Metadata seems to not have two groups; see original error message:\n{errorMsg}")
        elif "Y variable is expected to only have" in errorMsg:
            raise ValueError(f"Metadata seems to have unrecognised group labels; see original error message:\n{errorMsg}")
        elif "has no variants" in errorMsg:
            raise ValueError(f"There seem to be no variants available for sPLS-DA; see original error message:\n{errorMsg}")
        elif "Error with ordering of encoded and selected features" in errorMsg:
            raise Exception(f"Unexpected internal logic error in sPLS-DA code; see original error message:\n{errorMsg}")
        else:
            raise Exception(("run_windowed_splsda encountered an unrecoverable problem when processing " + 
                             f"'{encodedVcfFile}'; have a look at the stderr to make sense of this:\n'{errorMsg}'"))

def run_integrative_splsda(callRdataFile, depthRdataFile, outputSelected, outputPredictions,
                           scriptLocation, threads=1, nrepeat=10, maxiters=1000):
    '''
    Calls the integrative_splsda.R script to run sPLS-DA on the outputs of call and depth sPLS-DA.
    
    Parameters:
        callRdataFile -- a string indicating the location of the call sPLS-DA RData file
        depthRdataFile -- a string indicating the location of the depth sPLS-DA RData file
        outputSelected -- a string indicating the location of the output selected features file
        outputPredictions -- a string indicating the location of the output predictions file
        scriptLocation -- a string indicating the location of the integrative_plsda.R script
        threads -- (OPTIONAL) an integer indicating the number of threads to use (default is 1)
        nrepeat -- (OPTIONAL) an integer indicating the number of repeats for
                   stability analysis (default is 10)
        maxiters -- (OPTIONAL) an integer indicating the maximum number of iterations when
                    tuning sPLS-DA (default is 1000)
    '''
    # Format command
    cmd = ["Rscript", scriptLocation, callRdataFile, depthRdataFile, outputSelected, outputPredictions,
           "--threads", str(threads), "--nrepeat", str(nrepeat), "--maxiters", str(maxiters)]
    
    # Run bcftools index
    run_Rscript = subprocess.Popen(" ".join(cmd), shell=True,
                                   stdout = subprocess.DEVNULL,
                                   stderr = subprocess.PIPE)
    rout, rerr = run_Rscript.communicate()
    
    # Check for errors
    if run_Rscript.returncode == 0:
        return None
    else:
        errorMsg = rerr.decode("utf-8").rstrip("\r\n ")
        if "Integrative sPLS-DA is not possible or necessary" in errorMsg:
            print(f"# Note: {errorMsg}; this does not constitute an error, and program operations will continue sans integrative predictions")
        elif "X sample names do not match across datasets" in errorMsg:
            raise Exception(f"'call' and 'depth' result sample names do not match (are you using combining results somehow?); see original error message:\n{errorMsg}")
        elif "Y variables do not match across datasets" in errorMsg:
            raise Exception(f"'call' and 'depth' result metadata do not match (are you using combining results somehow?); see original error message:\n{errorMsg}")
        else:
            raise Exception(("run_integrative_splsda encountered an unhandled situation when processing " + 
                             f"'{callRdataFile}' and '{depthRdataFile}'; have a look at the stderr to " + 
                             f"make sense of this:\n'{errorMsg}'"))

def parse_selected_to_windowed_ncls(selectedFileName, exclusionsNCLS=None):
    '''
    Parameters:
        selectedFileName -- a file name indicating the location of the selected variants file
        exclusionsNCLS -- (OPTIONAL); a RangeNCLS object indexing ranges within which data points
                          should be eliminated
    Returns:
        rangeNCLS -- a RangeNCLS object containing statistical values indexed
                     by chromosome, start, and end position
    '''
    EXPECTED_HEADER = ["chrom", "start", "end", "stability", "abs_loading", "direction"]
    
    # Parse the selected file
    rangeNCLS = RangeNCLS()
    with open(selectedFileName, "r") as fileIn:
        # Read and validate the header
        header = fileIn.readline().strip().split("\t")
        if header != EXPECTED_HEADER:
            raise ValueError(f"Invalid header in file '{selectedFileName}', should be: {EXPECTED_HEADER}")
        
        # Store each line in the windowedNCLS object
        for line in fileIn:
            # Parse relevant details
            chrom, start, end, stability, abs_loading, direction = line.strip().split("\t")
            try:
                start = int(float(start))
                end = int(float(end))
            except:
                raise ValueError(f"Start '{start}' and/or end '{end}' is not an integer in file '{selectedFileName}'")
            try:
                stability = float(stability)
            except:
                raise ValueError(f"Stability '{stability}' is not a float in file '{selectedFileName}'")
            try:
                abs_loading = float(abs_loading)
            except:
                raise ValueError(f"abs_loading '{abs_loading}' is not a float in file '{selectedFileName}'")
            
            # Skip if we are excluding this position
            if exclusionsNCLS != None:
                if exclusionsNCLS.is_overlapping(chrom, start, end):
                    continue
            
            # Compute the stability*abs_loading value
            statProduct = stability * abs_loading
            
            # Store the values in the RangeNCLS object
            rangeNCLS.add(chrom, start, end, statProduct)
    
    rangeNCLS.build()
    return rangeNCLS

def parse_integrated_to_windowed_ncls(selectedFileName):
    '''
    Parameters:
        selectedFileName -- a file name indicating the location of the selected variants file
    Returns:
        nclsList -- a list of two WindowedNCLS objects:
            callWindowedNCLS -- a WindowedNCLS object containing statistical values indexed by chromosome
                                and position specifically for call variants that were selected
            depthWindowedNCLS -- a WindowedNCLS object containing statistical values indexed by chromosome
                                and position specifically for depth variants that were selected
    '''
    EXPECTED_HEADER = ["chrom", "start", "end", "type", "stability", "abs_loading", "direction"]
    
    # Parse the selected file
    callRangeNCLS = RangeNCLS()
    depthRangeNCLS = RangeNCLS()
    with open(selectedFileName, "r") as fileIn:
        # Read and validate the header
        header = fileIn.readline().strip().split("\t")
        if header != EXPECTED_HEADER:
            raise ValueError(f"Invalid header in file '{selectedFileName}', should be: {EXPECTED_HEADER}")
        
        # Store each line in the windowedNCLS object
        for line in fileIn:
            # Parse relevant details
            chrom, start, end, featuretype, stability, abs_loading, direction = line.strip().split("\t")
            try:
                start = int(float(start))
                end = int(float(end))
            except:
                raise ValueError(f"Start '{start}' and/or end '{end}' is not an integer in file '{selectedFileName}'")
            try:
                stability = float(stability)
            except:
                raise ValueError(f"Stability '{stability}' is not a float in file '{selectedFileName}'")
            try:
                abs_loading = float(abs_loading)
            except:
                raise ValueError(f"abs_loading '{abs_loading}' is not a float in file '{selectedFileName}'")
            
            # Compute the stability*abs_loading value
            statProduct = stability * abs_loading
            
            # Store the values in the RangeNCLS object
            if featuretype == "call":
                callRangeNCLS.add(chrom, start, end, statProduct)
            else:
                depthRangeNCLS.add(chrom, start, end, statProduct)
    
    callRangeNCLS.build()
    depthRangeNCLS.build()
    return callRangeNCLS, depthRangeNCLS

def parse_ber_to_windowed_ncls(berFileName, balancedAccuracy=True, exclusionsNCLS=None):
    '''
    Parameters:
        berFileName -- a file name indicating the location of the BER windows file
        balancedAccuracy -- (OPTIONAL); a boolean indicating whether to convert the
                            BER to balanced accuracy (default is True)
        exclusionsNCLS -- (OPTIONAL); a RangeNCLS object indexing ranges within which data points
                          should be eliminated
    Returns:
        rangeNCLS -- a RangeNCLS object containing statistical values indexed
                     by chromosome, start, and end position
    '''
    EXPECTED_HEADER = ["chrom", "start", "end", "BER"]
    
    # Parse the selected file
    rangeNCLS = RangeNCLS()
    with open(berFileName, "r") as fileIn:
        # Read and validate the header
        header = fileIn.readline().strip().split("\t")
        if header != EXPECTED_HEADER:
            raise ValueError(f"Invalid header in file '{berFileName}', should be: {EXPECTED_HEADER}")
        
        # Store each line in the RangeNCLS object
        for line in fileIn:
            # Parse relevant details
            chrom, start, end, ber = line.strip().split("\t")
            try:
                start = int(float(start))
                end = int(float(end))
            except:
                raise ValueError(f"Start '{start}' and/or end '{end}' is not an integer in file '{berFileName}'")
            try:
                ber = float(ber)
            except:
                raise ValueError(f"BER '{ber}' is not a float in file '{berFileName}'")
            
            # Skip if we are excluding this position
            if exclusionsNCLS != None:
                if exclusionsNCLS.is_overlapping(chrom, start, end):
                    continue
            
            # Convert BER to balanced accuracy if needed
            if balancedAccuracy:
                ber = 1 - (ber*2)
            
            # Store the values in the RangeNCLS object
            rangeNCLS.add(chrom, start, end, ber)
    
    rangeNCLS.build()
    return rangeNCLS
