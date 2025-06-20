#! python3

import os, sys, unittest, time, math, shutil, subprocess
import numpy as np
from collections import Counter

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.parsing import parse_metadata, vcf_header_to_metadata_validation, parse_vcf_genotypes, \
    parse_vcf_stats, parse_samtools_depth_tsv, parse_binned_tsv, read_gz_file
from modules.ncls import WindowedNCLS
from modules.ed import parse_vcf_for_ed, gt_median_adjustment, filter_impossible_genotypes, \
    calculate_allele_frequency_ed, calculate_genotype_frequency_ed, calculate_inheritance_ed
from modules.depth import get_median_value, convert_depth_to_alleles, call_cnvs_from_depth
from modules.samtools_handling import depth_to_histoDict
from modules.splsda import recode_variant, recode_cnv, recode_vcf
from modules.gff3 import GFF3Graph

# Specify data locations
dataDir = os.path.join(os.getcwd(), "data")
metadataFile = os.path.join(dataDir, "metadata.tsv")
baseDir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

MAXIMAL_SEGREGATION = 1.4142135623730951 # sqrt(2)

# Define utility functions
def run_subprocess(command):
    '''
    Parameters:
        command -- a list of strings representing the command to run
    '''
    process = subprocess.Popen(" ".join(command), shell = True,
                               stdout = subprocess.PIPE,
                               stderr = subprocess.PIPE)
    stdout, stderr = process.communicate()
    return process.returncode, stdout.decode("utf-8").rstrip("\r\n "), stderr.decode("utf-8").rstrip("\r\n ")

# Define unit tests
class TestParsing(unittest.TestCase):
    def test_parse_metadata(self):
        "Test parsing metadata file as many subsequent tests depend on it"
        # Arrange & Act
        metadataDict = parse_metadata(metadataFile)
        
        # Assert
        self.assertIsInstance(metadataDict, dict, "Metadata should be a dictionary")
        self.assertIn("group1", metadataDict, "Metadata should contain 'group1'")
        self.assertIn("group2", metadataDict, "Metadata should contain 'group2'")
        self.assertEqual(len(metadataDict), 2, "Should have 2 groups in metadata")
        self.assertEqual(len(metadataDict["group1"]), 11, "Group1 should have 11 entries")
        self.assertEqual(len(metadataDict["group2"]), 41, "Group2 should have 41 entries")
    
    def test_vcf_header_to_metadata_validation(self):
        "Test that the VCF header matches the metadata"
        # Arrange
        metadataDict = parse_metadata(metadataFile)
        vcfFile = os.path.join(dataDir, "deletions.1.vcf")
        with open(vcfFile, "r") as fileIn:
            for line in fileIn:
                l = line.strip('\r\n\t "')
                sl = l.replace('"', '').split("\t")
                if line.startswith("#CHROM"):
                    samples = sl[9:]
        
        # Act & Assert
        ## This will raise an error if the samples in the VCF header do not match the metadata
        vcf_header_to_metadata_validation(samples, metadataDict, strict=True, quiet=True)
    
    def test_parse_vcf_genotypes(self):
        "Test that VCF parsing works correctly"
        # Arrange
        metadataDict = parse_metadata(metadataFile)
        vcfFile = os.path.join(dataDir, "deletions.1.vcf")
        with open(vcfFile, "r") as fileIn:
            for line in fileIn:
                l = line.strip('\r\n\t "')
                sl = l.replace('"', '').split("\t")
                if line.startswith("#CHROM"):
                    samples = sl[9:] # This gives us the ordered sample IDs
                if not line.startswith("#"):
                    break
        
        # Act
        formatField = sl[8]
        sampleFields = sl[9:]
        snpDict = parse_vcf_genotypes(formatField, sampleFields, samples)
        
        # Assert
        self.assertIsInstance(snpDict, dict, "SNP dictionary should be a dictionary")
        self.assertEqual(len(snpDict), 52, "SNP dictionary should have 52 samples")
        self.assertEqual(snpDict["S48"], [0, 1], "Sample S48 should have genotype [0, 1]")
        self.assertEqual(snpDict["S56"], [0, 1], "Sample S56 should have genotype [0, 1]")
        self.assertEqual(snpDict["S59"], [0, 1], "Sample S59 should have genotype [0, 1]")
        self.assertEqual(snpDict["S44"], [0, 1], "Sample S44 should have genotype [0, 1]")
        self.assertEqual(snpDict["S77"], [0, 1], "Sample S77 should have genotype [0, 1]")
        self.assertEqual(snpDict["S64"], [0, 1], "Sample S64 should have genotype [0, 1]")
        self.assertEqual(snpDict["S54"], [0, 1], "Sample S54 should have genotype [0, 1]")
        self.assertEqual(snpDict["S76"], [0, 1], "Sample S76 should have genotype [0, 1]")
        self.assertEqual(snpDict["S50"], [0, 1], "Sample S50 should have genotype [0, 1]")
        self.assertEqual(snpDict["S69"], [1, 1], "Sample S69 should have genotype [1, 1]")
        self.assertEqual(snpDict["S95"], [0, 1], "Sample S95 should have genotype [0, 1]")
    
    def test_parse_vcf_stats(self):
        "Test that statistics generated from VCF parsing are correct"
        # Arrange
        vcfFile = os.path.join(dataDir, "deletions.1.vcf")
        
        # Act
        positions, variants, samples, contigs = parse_vcf_stats(vcfFile, refAllele="1")
        
        # Assert
        self.assertEqual(positions, 10, "Should have 10 positions in the VCF")
        self.assertEqual(variants, 9, "Should have 9 CNVs in the VCF")
        self.assertEqual(len(samples), 52, "Should have 52 samples in the VCF")
        self.assertEqual(len(contigs), 1, "Should have 1 contig in the VCF")
    
    def test_parse_samtools_depth_tsv(self):
        "Test that parsing a samtools depth TSV file works correctly"
        # Arrange
        depthTsvFile = os.path.join(dataDir, "depth.1.tsv")
        
        # Act
        depthData = list(parse_samtools_depth_tsv(depthTsvFile))
        
        # Assert
        self.assertEqual(len(depthData), 10, "Should have 10 positions in the depth dictionary")
        self.assertEqual(len(depthData[0]), 3, "Depth file should have 3 columns")
        self.assertEqual(depthData[0][0], "C.glau_01", "Depth file should contain 'C.glau_01'")
    
    def test_parse_binned_tsv(self):
        "Test that parsing a binned TSV file works correctly"
        # Arrange
        binnedTsvFile = os.path.join(dataDir, "depth.binned.1.tsv")
        
        # Act
        histoDict = parse_binned_tsv(binnedTsvFile)
        
        # Assert
        self.assertEqual(len(histoDict), 1, "Should have 1 contig in the binned depth file")
        self.assertIn("C.glau_01", histoDict, "Binned depth file should contain 'C.glau_01'")
        self.assertEqual(len(histoDict["C.glau_01"]), 18, "Should have 18 bins for 'C.glau_01'")

class TestDepth(unittest.TestCase):
    def test_get_median_value_1(self):
        "Test that the median value calculation works correctly with zeros"
        # Arrange
        depthList = np.array([0, 1, 2, 3, 4, 5])
        
        # Act
        medianValue = get_median_value(depthList)
        
        # Assert
        self.assertEqual(medianValue, 2.5, "Median of [0, 1, 2, 3, 4, 5] should be 2.5")
    
    def test_get_median_value_2(self):
        "Test that the median value calculation works correctly without zeros"
        # Arrange
        depthList = np.array([1, 2, 3, 4, 5])
        
        # Act
        medianValue = get_median_value(depthList)
        
        # Assert
        self.assertEqual(medianValue, 3.0, "Median of [1, 2, 3, 4, 5] should be 3.0")
    
    def test_get_median_value_3(self):
        "Test that the median value calculation works correctly when all values are 0"
        # Arrange
        depthList = np.array([0, 0, 0, 0, 0])
        
        # Act
        medianValue = get_median_value(depthList)
        
        # Assert
        self.assertEqual(medianValue, 1, "Median should fall back to 1 if all values are 0")
    
    def test_get_median_value_4(self):
        "Test that the median value calculation works correctly when falling back to a non-zero value"
        # Arrange
        depthList = np.array([0, 0, 0, 5, 0, 0])
        
        # Act
        medianValue = get_median_value(depthList)
        
        # Assert
        self.assertEqual(medianValue, 5.0, "Median should fall back to only non-zero value if all others are 0")
    
    def test_convert_depth_to_alleles(self):
        "Test that median-normalisation works correctly for deletions"
        # Arrange
        binnedTsvFile = os.path.join(dataDir, "depth.binned.1.tsv")
        histoDict = parse_binned_tsv(binnedTsvFile)
        truth = np.array([0., 0., 0., 0., 1., 1., 2., 2., 2., 2., 2., 2., 3., 3., 4., 4., 4., 4.])
        
        # Act
        alleles = convert_depth_to_alleles(histoDict["C.glau_01"])
        
        # Assert
        self.assertEqual(len(alleles), 18, "Should have 18 alleles")
        self.assertTrue(all(alleles == truth), "Predicted alleles should match the expected truth values")
    
    def test_samtools_binning_1(self):
        "Test that samtools binning works correctly with binSize 1"
        # Arrange
        depthFile = os.path.join(dataDir, "depth.1.tsv")
        lengthsDict = {"C.glau_01": 10}  # Example contig length
        binSize = 1
        truth = np.array([1,  2,  3,  4,  5,  6,  7,  8,  9, 10])
        
        # Act
        histoDict = depth_to_histoDict(depthFile, lengthsDict, binSize)
        
        # Assert
        self.assertIn("C.glau_01", histoDict, f"'{depthFile}' should contain 'C.glau_01'")
        self.assertEqual(len(histoDict["C.glau_01"]), 10, f"Should have 10 bins for 'C.glau_01' in {depthFile}")
        self.assertTrue(all(histoDict['C.glau_01'] == truth), f"Expected '{truth}' but got '{histoDict['C.glau_01']}'")
    
    def test_samtools_binning_2(self):
        "Test that samtools binning works correctly with divisible binSize"
        # Arrange
        depthFile = os.path.join(dataDir, "depth.1.tsv")
        lengthsDict = {"C.glau_01": 10}  # Example contig length
        binSize = 2
        truth = np.array([3,  7, 11, 15, 19])
        
        # Act
        histoDict = depth_to_histoDict(depthFile, lengthsDict, binSize)
        
        # Assert
        self.assertIn("C.glau_01", histoDict, f"'{depthFile}' should contain 'C.glau_01'")
        self.assertEqual(len(histoDict["C.glau_01"]), 5, f"Should have 5 bins for 'C.glau_01' in {depthFile}")
        self.assertTrue(all(histoDict['C.glau_01'] == truth), f"Expected '{truth}' but got '{histoDict['C.glau_01']}'")
    
    def test_samtools_binning_3(self):
        "Test that samtools binning works correctly with unevenly divisible binSize"
        # Arrange
        depthFile = os.path.join(dataDir, "depth.1.tsv")
        lengthsDict = {"C.glau_01": 10}  # Example contig length
        binSize = 3
        truth = np.array([6, 15, 24, 10])
        
        # Act
        histoDict = depth_to_histoDict(depthFile, lengthsDict, binSize)
        
        # Assert
        self.assertIn("C.glau_01", histoDict, f"'{depthFile}' should contain 'C.glau_01'")
        self.assertEqual(len(histoDict["C.glau_01"]), 4, f"Should have 4 bins for 'C.glau_01' in {depthFile}")
        self.assertTrue(all(histoDict['C.glau_01'] == truth), f"Expected '{truth}' but got '{histoDict['C.glau_01']}'")
    
    def test_samtools_binning_4(self):
        "Test that samtools binning works correctly with excess binSize"
        # Arrange
        depthFile = os.path.join(dataDir, "depth.1.tsv")
        lengthsDict = {"C.glau_01": 10}  # Example contig length
        binSize = 1000
        truth = np.array([55])
        
        # Act
        histoDict = depth_to_histoDict(depthFile, lengthsDict, binSize)
        
        # Assert
        self.assertIn("C.glau_01", histoDict, f"'{depthFile}' should contain 'C.glau_01'")
        self.assertEqual(len(histoDict["C.glau_01"]), 1, f"Should have 1 bin for 'C.glau_01' in {depthFile}")
        self.assertTrue(all(histoDict['C.glau_01'] == truth), f"Expected '{truth}' but got '{histoDict['C.glau_01']}'")
    
    def test_call_cnvs_from_depth(self):
        "Test that call_cnvs_from_depth() works correctly with a binned depth file"
        # Arrange
        depthFile = os.path.join(dataDir, "depth.binned.1.tsv")
        samplePairs = [["test", depthFile]]
        
        windowSize = 1000
        ploidy = 2
        
        workDir = os.path.join(dataDir, "tmp")
        outputFileName = os.path.join(workDir, "test_call_cnvs_from_depth.tsv.gz")
        
        num00 = 4
        num01 = 2
        num11 = 6
        num12 = 2
        num22 = 4
        
        # Arrange: cleanup any previous work directory
        if os.path.exists(workDir):
            shutil.rmtree(workDir)
        if not os.path.exists(workDir):
            os.makedirs(workDir)
        
        # Act
        call_cnvs_from_depth(samplePairs, outputFileName, windowSize, ploidy=ploidy)
        
        # Assert
        deletionContents = []
        with read_gz_file(outputFileName) as fileIn:
            for line in fileIn:
                deletionContents.append(line.strip())
        
        genotypes = []
        for line in deletionContents:
            if line.startswith("#"):
                continue
            sl = line.split("\t")
            genotypes.append(sl[-1])
        numGenotypes = Counter(genotypes)
        
        self.assertEqual(numGenotypes["0/0"], num00, f"Expected {num00} '0/0' genotypes but got {numGenotypes['0/0']}")
        self.assertEqual(numGenotypes["0/1"], num01, f"Expected {num01} '0/1' genotypes but got {numGenotypes['0/1']}")
        self.assertEqual(numGenotypes["1/1"], num11, f"Expected {num11} '1/1' genotypes but got {numGenotypes['1/1']}")
        self.assertEqual(numGenotypes["1/2"], num12, f"Expected {num12} '1/2' genotypes but got {numGenotypes['1/2']}")
        self.assertEqual(numGenotypes["2/2"], num22, f"Expected {num22} '2/2' genotypes but got {numGenotypes['2/2']}")
    
    def test_call_cnvs_from_depth_triploid(self):
        "Test that call_cnvs_from_depth() works correctly with a binned depth file"
        # Arrange
        depthFile = os.path.join(dataDir, "depth.binned.1.tsv")
        samplePairs = [["test", depthFile]]
        
        windowSize = 1000
        ploidy = 3
        
        workDir = os.path.join(dataDir, "tmp")
        outputFileName = os.path.join(workDir, "test_call_cnvs_from_depth.tsv.gz")
        
        num000 = 2
        num001 = 2
        num011 = 4
        num111 = 2
        num112 = 4
        num222 = 2
        
        # Arrange: cleanup any previous work directory
        if os.path.exists(workDir):
            shutil.rmtree(workDir)
        if not os.path.exists(workDir):
            os.makedirs(workDir)
        
        # Act
        call_cnvs_from_depth(samplePairs, outputFileName, windowSize, ploidy=ploidy)
        
        # Assert
        deletionContents = []
        with read_gz_file(outputFileName) as fileIn:
            for line in fileIn:
                deletionContents.append(line.strip())
        
        genotypes = []
        for line in deletionContents:
            if line.startswith("#"):
                continue
            sl = line.split("\t")
            genotypes.append(sl[-1])
        numGenotypes = Counter(genotypes)
        
        self.assertEqual(numGenotypes["0/0/0"], num000, f"Expected {num000} '0/0/0' genotypes but got {numGenotypes['0/0/0']}")
        self.assertEqual(numGenotypes["0/0/1"], num001, f"Expected {num001} '0/0/1' genotypes but got {numGenotypes['0/0/1']}")
        self.assertEqual(numGenotypes["0/1/1"], num011, f"Expected {num011} '0/1/1' genotypes but got {numGenotypes['0/1/1']}")
        self.assertEqual(numGenotypes["1/1/1"], num111, f"Expected {num111} '1/1/1' genotypes but got {numGenotypes['1/1/1']}")
        self.assertEqual(numGenotypes["1/1/2"], num112, f"Expected {num112} '1/1/2' genotypes but got {numGenotypes['1/1/2']}")
        self.assertEqual(numGenotypes["2/2/2"], num222, f"Expected {num222} '2/2/2' genotypes but got {numGenotypes['2/2/2']}")

class TestED(unittest.TestCase):
    def test_parse_vcf_for_ed(self):
        "Test parsing a VCF with deletions where isCNV is True and False"
        # Arrange
        metadataDict = parse_metadata(metadataFile)
        vcfFile = os.path.join(dataDir, "deletions.1.vcf")
        notCNV_truth = [0.478198599250326,0.5942205544782738,0.6083495165377463,
                        0.7877780740982955,0.6366455643787936,0.3733361028931871,
                        0.14110778338534205,0.396669657738795]
        isCNV_truth = [0.478198599250326,0.5942205544782738,0.6914281385881762,
                       0.888979035327655,0.6626410942564203,0.3733361028931871,
                       0.14110778338534205,0.396669657738795]
        # notCNV_truth = [0.7716249967879256,0.7586519135256539,0.6083495165377463, # these values
        #                 0.8019437941637759,0.6764821026903794,0.6951352124261518, # are for if
        #                 0.2324038163770089,0.79333931547759]
        # isCNV_truth = [0.7716249967879256,0.8644095090819802,0.7223550682708708, # we are using
        #                0.8988827108558033,0.7694428635984445,0.28221556677068416, # genotype frequency
        #                0.0,0.79333931547759]
        
        # Act & Assert
        notCNV_results = []
        for resultsDict in parse_vcf_for_ed(vcfFile, metadataDict, isCNV=False, ignoreIdentical=True, quiet=True):
            notCNV_results.append(resultsDict["alleleED"])
        
        isCNV_results = []
        for resultsDict in parse_vcf_for_ed(vcfFile, metadataDict, isCNV=True, ignoreIdentical=True, quiet=True):
            isCNV_results.append(resultsDict["alleleED"])
        
        # Assert
        for v1, v2 in zip(notCNV_results, notCNV_truth):
            self.assertAlmostEqual(v1, v2, places=5, msg=f"Expected {v2} but got {v1}")
        for v1, v2 in zip(isCNV_results, isCNV_truth):
            self.assertAlmostEqual(v1, v2, places=5, msg=f"Expected {v2} but got {v1}")
    
    def test_parse_vcf_for_ed_with_subset_metadata(self):
        "Test parsing a VCF with deletions before and after subsetting metadata"
        # Arrange #1
        metadataDict = parse_metadata(metadataFile)
        vcfFile = os.path.join(dataDir, "deletions.1.vcf")
        
       # Act & Assert
        notCNV_results_1 = []
        for resultsDict in parse_vcf_for_ed(vcfFile, metadataDict, isCNV=False, ignoreIdentical=True, quiet=True):
            notCNV_results_1.append(resultsDict["alleleED"])
        
        isCNV_results_1 = []
        for resultsDict in parse_vcf_for_ed(vcfFile, metadataDict, isCNV=True, ignoreIdentical=True, quiet=True):
            isCNV_results_1.append(resultsDict["alleleED"])
        
        # Arrange #2
        metadataDict["group1"] = ['S44', 'S48', 'S50', 'S54', 'S56', 'S59', 'S64', 'S69', 'S76', 'S77'] # drop S95
        
        # Act #2
        notCNV_results_2 = []
        for resultsDict in parse_vcf_for_ed(vcfFile, metadataDict, isCNV=False, ignoreIdentical=True, quiet=True):
            notCNV_results_2.append(resultsDict["alleleED"])
        
        isCNV_results_2 = []
        for resultsDict in parse_vcf_for_ed(vcfFile, metadataDict, isCNV=True, ignoreIdentical=True, quiet=True):
            isCNV_results_2.append(resultsDict["alleleED"])
        
        # Assert
        for v1, v2 in zip(notCNV_results_1[:-1], notCNV_results_2[:-1]): # last value will not change
            self.assertNotEqual(v1, v2, f"Subsetting metadata should change results when isCNV=False")
        for v1, v2 in zip(isCNV_results_1[:-2], isCNV_results_2[:-2]): # last two values will not change
            self.assertNotEqual(v1, v2, f"Subsetting metadata should change results when isCNV=True")
    
    def test_parse_vcf_for_ed_tetraploid(self):
        "Test parsing a VCF (with tetraploid variants) with deletions where isCNV is True and False"
        # Arrange
        metadataDict = parse_metadata(metadataFile)
        vcfFile = os.path.join(dataDir, "deletions.2.vcf")
        notCNV_truth = [0.478198599250326,0.5942205544782738,0.6083495165377463,
                        0.7877780740982955,0.6366455643787936,0.3733361028931871,
                        0.14110778338534205,0.396669657738795]
        isCNV_truth = [0.478198599250326,0.5942205544782738,0.6914281385881762,
                       0.888979035327655,0.6626410942564203,0.3733361028931871,
                       0.14110778338534205,0.396669657738795]
        # notCNV_truth = [0.7716249967879256,0.7586519135256539,0.6083495165377463,
        #                 0.8019437941637759,0.6764821026903794,0.6951352124261518,
        #                 0.2324038163770089,0.79333931547759]
        # isCNV_truth = [0.7716249967879256,0.8644095090819802,0.7223550682708708,
        #                0.8988827108558033,0.7694428635984445,0.28221556677068416,
        #                0.0,0.79333931547759]
        
        # Act & Assert
        notCNV_results = []
        for resultsDict in parse_vcf_for_ed(vcfFile, metadataDict, isCNV=False, ignoreIdentical=True, quiet=True):
            notCNV_results.append(resultsDict["alleleED"])
        
        isCNV_results = []
        for resultsDict in parse_vcf_for_ed(vcfFile, metadataDict, isCNV=True, ignoreIdentical=True, quiet=True):
            isCNV_results.append(resultsDict["alleleED"])
        
        # Assert
        for v1, v2 in zip(notCNV_results, notCNV_truth):
            self.assertAlmostEqual(v1, v2, places=5, msg=f"Expected {v2} but got {v1}")
        for v1, v2 in zip(isCNV_results, isCNV_truth):
            self.assertAlmostEqual(v1, v2, places=5, msg=f"Expected {v2} but got {v1}")
    
    def test_calculate_allele_frequency_ed_1(self):
        "Test that allele frequency ED is different for isCNV True and False under certain conditions"
        # Arrange
        g1Gt = [[2, 2],[0, 1],[0, 0],[1, 2],[0, 1],[0, 0],[0, 1],[0, 0],[0, 0],[1, 2],[0, 0]]
        g2Gt = [[0, 0],[2, 2],[0, 0],[0, 0],[0, 0],[0, 0],[0, 0],[0, 0],[0, 0],[0, 0],[0, 1],
                [1, 2],[0, 0],[2, 2],[1, 2],[0, 0],[1, 2],[2, 2],[0, 0],[1, 2],[0, 0],[1, 2],
                [0, 0],[0, 0],[0, 0],[0, 0],[0, 1],[1, 2],[1, 2],[2, 2],[1, 2],[1, 2],[2, 2],
                [0, 0],[0, 0],[0, 0],[0, 0],[0, 0],[0, 0],[2, 2],[1, 2]]
        g1Gt_cnv, g2Gt_cnv = gt_median_adjustment([g1Gt, g2Gt])
        
        # Act
        g1Alleles_1, g2Alleles_1, aED_1 = calculate_allele_frequency_ed(g1Gt_cnv, g2Gt_cnv)
        g1Alleles_2, g2Alleles_2, aED_2 = calculate_allele_frequency_ed(g1Gt, g2Gt)
        
        # Assert
        self.assertNotEqual(aED_1, aED_2, "ED should not be the same for isCNV True and False")
    
    def test_calculate_allele_frequency_ed_2(self):
        "Test that segregation ED is the same for isCNV True and False under certain conditions"
        # Arrange
        g1Gt = [[0, 1],[0, 1],[0, 1],[0, 1],[0, 1],[0, 1],[0, 1],[0, 1],[0, 1],[1, 1],[0, 1]]
        g2Gt = [[0, 0],[0, 1],[0, 1],[0, 1],[0, 0],[0, 0],[0, 1],[0, 1],[0, 1],[0, 0],[0, 1],
                [0, 1],[0, 0],[0, 1],[0, 1],[0, 0],[0, 0],[0, 0],[0, 1],[0, 1],[0, 1],[0, 0],
                [0, 0],[0, 0],[0, 0],[0, 0],[0, 1],[0, 0],[0, 1],[0, 0],[0, 0],[0, 0],[0, 0],
                [0, 1],[0, 0],[0, 1],[0, 0],[0, 0],[0, 0],[0, 0],[0, 0]]
        g1Gt_cnv, g2Gt_cnv = gt_median_adjustment([g1Gt, g2Gt])
        
        # Act
        g1Alleles_1, g2Alleles_1, aED_1 = calculate_allele_frequency_ed(g1Gt_cnv, g2Gt_cnv)
        g1Alleles_2, g2Alleles_2, aED_2 = calculate_allele_frequency_ed(g1Gt, g2Gt)
        
        # Assert
        self.assertEqual(aED_1, aED_2, "ED should be the same for isCNV True and False")
    
    def test_genotype_versus_allele_frequency(self):
        "Test that genotype frequency differentiates heterozygotes from a mixture of different homozygotes"
        # Arrange
        parentsGT = [ [0, 1], [0, 1] ]
        g1Gt = [[0, 1],[0, 1],[0, 1],[0, 1],[0, 1],[0, 1]]
        g2Gt = [[0, 0],[1, 1],[0, 0],[1, 1],[0, 0],[1, 1]]
        gEDTruth = 1.224744871391589 # not maximal since it is a mixture of homozygotes
        aEDTruth = 0
        
        # Act
        g1Alleles_1, g2Alleles_1, aED = calculate_allele_frequency_ed(g1Gt, g2Gt)
        g1Alleles_2, g2Alleles_2, gED = calculate_genotype_frequency_ed(g1Gt, g2Gt)
        
        # Assert
        self.assertEqual(gED, gEDTruth, f"Expected genotypes ED to {gEDTruth} but got {gED}")
        self.assertEqual(aED, aEDTruth, f"Expected alleles ED to be {aEDTruth} but got {aED}")
    
    def test_inheritance_versus_allele_frequency_equality_1(self):
        "Both ED methods should be equal for zero segregation"
        # Arrange
        parentsGT = [ [0, 1], [0, 1] ]
        g1Gt = [[0, 1], [0, 1], [0, 1], [0, 1]]
        g2Gt = [[0, 1], [0, 1], [0, 1], [0, 1]]
        truth = 0
        
        # Act
        g1Alleles_1, g2Alleles_1, aED = calculate_allele_frequency_ed(g1Gt, g2Gt)
        g1Gt, g2Gt = filter_impossible_genotypes(g1Gt, g2Gt, parentsGT)
        g1Filtered, g2Filtered, iED = calculate_inheritance_ed(g1Gt, g2Gt, parentsGT)
        
        # Assert
        self.assertEqual(aED, truth, "Allele ED should be 0 for identical genotypes")
        self.assertEqual(iED, truth, "Inheritance ED should be 0 for identical genotypes")
        self.assertEqual(aED, iED, "Both ED methods should be equal for zero segregation")
    
    def test_inheritance_versus_allele_frequency_equality_2(self):
        "Both ED methods should be equal for maximal segregation"
        # Arrange
        parentsGT = [ [0, 1], [0, 1] ]
        g1Gt = [[0, 0], [0, 0], [0, 0], [0, 0]]
        g2Gt = [[1, 1], [1, 1], [1, 1], [1, 1]]
        truth = MAXIMAL_SEGREGATION
        
        # Act
        g1Alleles_1, g2Alleles_1, aED = calculate_allele_frequency_ed(g1Gt, g2Gt)
        g1Gt, g2Gt = filter_impossible_genotypes(g1Gt, g2Gt, parentsGT)
        g1Filtered, g2Filtered, iED = calculate_inheritance_ed(g1Gt, g2Gt, parentsGT)
        
        # Assert
        self.assertAlmostEqual(aED, truth, places=5, msg="Alleles ED should be ~1.41421 for identical genotypes")
        self.assertAlmostEqual(iED, truth, places=5, msg="Inheritance ED should be ~1.41421 for identical genotypes")
        self.assertEqual(aED, iED, "Both ED methods should be equal for maximal segregation")
    
    def test_calculate_inheritance_ed_dip_dip_parents(self):
        # Arrange
        parentsGT = [ [0, 1], [0, 2] ] # represents "A/T" and "A/G" parents
        g1Gt = [[0, 2], [1, 2], [0, 2], [0, 0]] # represents "A/G", "T/G", "A/G", "A/A" genotypes in group 1
        g2Gt = [[1, 2], [0, 1], [0, 0], [0, 0]] # represents "T/G", "A/T", "A/A", "A/A" genotypes in group 2
        
        # truth = 0.7905694150420949 # calculated by hand on the plane from WA to Brisbane!
        truth = 0.5590169943749475 # after adjustment of dividing sum by 2 prior to sqrt() operation
        
        # Act
        g1Gt, g2Gt = filter_impossible_genotypes(g1Gt, g2Gt, parentsGT)
        g1Filtered, g2Filtered, iED = calculate_inheritance_ed(g1Gt, g2Gt, parentsGT)
                
        # Assert
        self.assertAlmostEqual(iED, truth, places=5, 
                             msg=f"Expected ED to be approximately {truth} but got {iED}")
        self.assertEqual(g1Filtered, 8, "Expected 8 alleles in group 1")
        self.assertEqual(g2Filtered, 8, "Expected 8 alleles in group 2")
    
    def test_calculate_inheritance_ed_tet_tet_parents_1(self):
        "Non-segregation should still give ED==0 with tetraploid parents and samples"
        # Arrange
        parentsGT = [ [0, 0, 1, 1], [0, 0, 2, 2] ]
        g1Gt = [[0, 0, 1, 2], [0, 0, 1, 2], [0, 0, 1, 2], [0, 0, 1, 2]]
        g2Gt = [[0, 0, 1, 2], [0, 0, 1, 2], [0, 0, 1, 2], [0, 0, 1, 2]]
        truth = 0
        
        # Act
        g1Gt, g2Gt = filter_impossible_genotypes(g1Gt, g2Gt, parentsGT)
        g1Filtered, g2Filtered, iED = calculate_inheritance_ed(g1Gt, g2Gt, parentsGT)
                
        # Assert
        self.assertEqual(iED, truth, f"Expected ED to be zero but got {iED}")
    
    def test_calculate_inheritance_ed_tet_tet_parents_2(self):
        "Maximal-segregation should still give ED==1.41421 with tetraploid parents and samples"
        # Arrange
        parentsGT = [ [0, 1, 2, 3], [0, 1, 4, 5] ]
        g1Gt = [[0, 1, 2, 4], [0, 1, 2, 4], [0, 1, 2, 4], [0, 1, 2, 4]]
        g2Gt = [[0, 1, 3, 5], [0, 1, 3, 5], [0, 1, 3, 5], [0, 1, 3, 5]]
        truth = MAXIMAL_SEGREGATION
        
        # Act
        g1Gt, g2Gt = filter_impossible_genotypes(g1Gt, g2Gt, parentsGT)
        g1Filtered, g2Filtered, iED = calculate_inheritance_ed(g1Gt, g2Gt, parentsGT)
        
        # Assert
        self.assertEqual(iED, truth, f"Expected ED to be zero but got {iED}")
    
    def test_calculate_inheritance_ed_tet_tet_parents_3(self):
        "Half-shuffled samples from the above test should give ED==0 with tetraploid parents and samples"
        # Arrange
        parentsGT = [ [0, 1, 2, 3], [0, 1, 4, 5] ]
        g1Gt = [[0, 1, 2, 4], [0, 1, 2, 4], [0, 1, 3, 5], [0, 1, 3, 5]]
        g2Gt = [[0, 1, 3, 5], [0, 1, 3, 5], [0, 1, 2, 4], [0, 1, 2, 4]]
        truth = 0
        
        # Act
        g1Gt, g2Gt = filter_impossible_genotypes(g1Gt, g2Gt, parentsGT)
        g1Filtered, g2Filtered, iED = calculate_inheritance_ed(g1Gt, g2Gt, parentsGT)
        
        # Assert
        self.assertEqual(iED, truth, f"Expected ED to be zero but got {iED}")
    
    def test_calculate_inheritance_ed_with_impossible_progeny_1(self):
        "Test for impossible progeny (alleles do not exist in parents)"
        # Arrange
        parentsGT = [ [0, 0, 1, 1], [0, 0, 2, 2] ]
        g1Gt = [[0, 1, 2, 4], [0, 1, 2, 4], [0, 1, 3, 5], [0, 1, 3, 5]]
        g2Gt = [[0, 1, 3, 5], [0, 1, 3, 5], [0, 1, 2, 4], [0, 1, 2, 4]]
        edTruth = 0
        postFilterTruth = 0
        
        # Act
        g1Gt, g2Gt = filter_impossible_genotypes(g1Gt, g2Gt, parentsGT)
        g1Filtered, g2Filtered, iED = calculate_inheritance_ed(g1Gt, g2Gt, parentsGT)
        
        # Assert
        self.assertEqual(iED, edTruth, f"Expected ED to be {edTruth} but got {iED}")
        self.assertEqual(g1Filtered, postFilterTruth, f"Expected {postFilterTruth} alleles in group 1 after filtering but got {g1Filtered}")
        self.assertEqual(g2Filtered, postFilterTruth, f"Expected {postFilterTruth} alleles in group 2 after filtering but got {g2Filtered}")
    
    def test_calculate_inheritance_ed_with_impossible_progeny_2(self):
        "Test for impossible progeny (allele combination cannot be inherited from parents)"
        # Arrange
        parentsGT = [ [1, 1], [0, 1] ]
        g1Gt = [[0, 0],[0, 0],[0, 0],[0, 0]]
        g2Gt = [[1, 1],[0, 1],[1, 1],[0, 1]]
        edTruth = 0
        preFilterTruth = 8
        postFilterTruth = 0 # only group 1 should be filtered out
        
        # Act
        g1Alleles, g2Alleles, aED = calculate_allele_frequency_ed(g1Gt, g2Gt) # just to get the allele counts
        g1Gt, g2Gt = filter_impossible_genotypes(g1Gt, g2Gt, parentsGT)
        g1Filtered, g2Filtered, iED = calculate_inheritance_ed(g1Gt, g2Gt, parentsGT)
        
        # Assert
        self.assertEqual(iED, edTruth, f"Expected ED to be {edTruth} but got {iED}")
        self.assertEqual(g1Alleles, preFilterTruth, f"Expected {preFilterTruth} alleles in group 1 before filtering but got {g1Alleles}")
        self.assertEqual(g2Alleles, preFilterTruth, f"Expected {preFilterTruth} alleles in group 2 before filtering but got {g2Alleles}")
        self.assertEqual(g1Filtered, postFilterTruth, f"Expected {postFilterTruth} alleles in group 1 after filtering but got {g1Filtered}")
        self.assertEqual(g2Filtered, preFilterTruth, f"Expected {preFilterTruth} alleles in group 2 after filtering but got {g2Filtered}")
    
    def test_calculate_inheritance_ed_with_impossible_progeny_3(self):
        "Test for impossible progeny (alleles could only be from a clone)"
        # Arrange
        parentsGT = [ [2, 1], [0, 0] ]
        g1Gt = [[2, 1]]
        g2Gt = [[0, 0]]
        edTruth = 0
        preFilterTruth = 2
        postFilterTruth = 0
        
        # Act
        g1Alleles, g2Alleles, aED = calculate_allele_frequency_ed(g1Gt, g2Gt) # just to get the allele counts
        g1Gt, g2Gt = filter_impossible_genotypes(g1Gt, g2Gt, parentsGT)
        g1Filtered, g2Filtered, iED = calculate_inheritance_ed(g1Gt, g2Gt, parentsGT)
        
        # Assert
        self.assertEqual(iED, edTruth, f"Expected ED to be {edTruth} but got {iED}")
        self.assertEqual(g1Alleles, preFilterTruth, f"Expected {preFilterTruth} alleles in group 1 before filtering but got {g1Alleles}")
        self.assertEqual(g2Alleles, preFilterTruth, f"Expected {preFilterTruth} alleles in group 2 before filtering but got {g2Alleles}")
        self.assertEqual(g1Filtered, postFilterTruth, f"Expected {postFilterTruth} alleles in group 1 after filtering but got {g1Filtered}")
        self.assertEqual(g2Filtered, postFilterTruth, f"Expected {postFilterTruth} alleles in group 2 after filtering but got {g2Filtered}")
    
    def test_calculate_ed_for_results_stability_1(self):
        "Tests done with the program in a state where the results are trustworthy; future changes should not affect the results"
        # Arrange
        parentsGT = [ [0, 1], [1, 1] ]
        g1Gt = [[0, 1], [0, 1], [0, 1], [0, 1]]
        g2Gt = [[1, 1], [0, 1], [0, 1], [1, 1]]
        inheritanceTruth = 0.42898458920779164
        alleleTruth = 0.3535533905932738
        genotypeTruth = 0.7071067811865476
        
        # Act
        g1Alleles, g2Alleles, aED = calculate_allele_frequency_ed(g1Gt, g2Gt) # just to get the allele counts
        _, _, gED = calculate_genotype_frequency_ed(g1Gt, g2Gt)
        g1Gt, g2Gt = filter_impossible_genotypes(g1Gt, g2Gt, parentsGT)
        g1Filtered, g2Filtered, iED = calculate_inheritance_ed(g1Gt, g2Gt, parentsGT)
        
        # Assert
        self.assertAlmostEqual(aED, alleleTruth, places=5, msg=f"Expected allele frequency ED to be {alleleTruth} but got {aED}")
        self.assertAlmostEqual(gED, genotypeTruth, places=5, msg=f"Expected genotype frequency ED to be {genotypeTruth} but got {gED}")
        self.assertAlmostEqual(iED, inheritanceTruth, places=5, msg=f"Expected inheritance ED to be {inheritanceTruth} but got {iED}")
    
    def test_calculate_ed_for_results_stability_2(self):
        "Tests done with the program in a state where the results are trustworthy; future changes should not affect the results"
        # Arrange
        parentsGT = [ [0, 0], [1, 1] ]
        g1Gt = [[0, 1], [0, 1], [0, 1], [0, 1]]
        g2Gt = [[1, 1], [0, 1], [0, 1], [1, 1]]
        inheritanceTruth = 0.0 # 0.25 without parental filtering
        alleleTruth = 0.3535533905932738
        genotypeTruth = 0.7071067811865476
        
        # Act
        g1Alleles, g2Alleles, aED = calculate_allele_frequency_ed(g1Gt, g2Gt) # just to get the allele counts
        _, _, gED = calculate_genotype_frequency_ed(g1Gt, g2Gt)
        g1Gt, g2Gt = filter_impossible_genotypes(g1Gt, g2Gt, parentsGT)
        g1Filtered, g2Filtered, iED = calculate_inheritance_ed(g1Gt, g2Gt, parentsGT)
        
        # Assert
        self.assertAlmostEqual(aED, alleleTruth, places=5, msg=f"Expected allele frequency ED to be {alleleTruth} but got {aED}")
        self.assertAlmostEqual(gED, genotypeTruth, places=5, msg=f"Expected genotype frequency ED to be {genotypeTruth} but got {gED}")
        self.assertAlmostEqual(iED, inheritanceTruth, places=5, msg=f"Expected inheritance ED to be {inheritanceTruth} but got {iED}")
    
    def test_calculate_ed_for_results_stability_3(self):
        "Tests done with the program in a state where the results are trustworthy; future changes should not affect the results"
        # Arrange
        parentsGT = [ [0, 0], [0, 1] ]
        g1Gt = [[0, 0], [0, 1], [0, 0], [0, 1]]
        g2Gt = [[1, 1], [0, 1], [0, 1], [0, 0]]
        inheritanceTruth = 0.14299486306926387
        alleleTruth = 0.3535533905932738
        genotypeTruth = 0.3535533905932738
        
        # Act
        g1Alleles, g2Alleles, aED = calculate_allele_frequency_ed(g1Gt, g2Gt) # just to get the allele counts
        _, _, gED = calculate_genotype_frequency_ed(g1Gt, g2Gt)
        g1Gt, g2Gt = filter_impossible_genotypes(g1Gt, g2Gt, parentsGT)
        g1Filtered, g2Filtered, iED = calculate_inheritance_ed(g1Gt, g2Gt, parentsGT)
        
        # Assert
        self.assertAlmostEqual(aED, alleleTruth, places=5, msg=f"Expected allele frequency ED to be {alleleTruth} but got {aED}")
        self.assertAlmostEqual(gED, genotypeTruth, places=5, msg=f"Expected genotype frequency ED to be {genotypeTruth} but got {gED}")
        self.assertAlmostEqual(iED, inheritanceTruth, places=5, msg=f"Expected inheritance ED to be {inheritanceTruth} but got {iED}")
    
    def test_calculate_inheritance_ed_with_empty_group(self):
        "Test for impossible progeny (alleles could only be from a clone)"
        # Arrange
        parentsGT = [ [1, 2], [0, 0] ]
        g1Gt = []
        g2Gt = [[0, 1]]
        edTruth = 0
        g1Truth = 0
        g2Truth = 2
        
        # Act
        g1Gt, g2Gt = filter_impossible_genotypes(g1Gt, g2Gt, parentsGT)
        g1Filtered, g2Filtered, iED = calculate_inheritance_ed(g1Gt, g2Gt, parentsGT)
        
        # Assert
        self.assertEqual(iED, edTruth, f"Expected ED to be zero but got {iED}")
        self.assertEqual(g1Filtered, g1Truth, f"Expected {g1Truth} alleles in group 1 but got {g1Filtered}")
        self.assertEqual(g2Filtered, g2Truth, f"Expected {g2Truth} alleles in group 2 but got {g2Filtered}")
    
    def test_calculate_gt_median_adjustment_1(self):
        "Test that median adjustment for CNVs works correctly (groups segregate evenly)"
        # Arrange
        g1Gt = [[0, 0], [0, 0], [0, 0], [0, 0]]
        g2Gt = [[100, 100], [100, 100], [100, 100], [100, 100]]
        g1Truth = 0 # below median == 0
        g2Truth = 2 # above median == 2
        
        # Act
        g1AdjGt, g2AdjGt = gt_median_adjustment([g1Gt, g2Gt])
        
        # Assert
        self.assertTrue(all([ allele == g1Truth for gt in g1AdjGt for allele in gt ]), f"Expected adjusted group 1 to be {g1Truth}")
        self.assertTrue(all([ allele == g2Truth for gt in g2AdjGt for allele in gt ]), f"Expected adjusted group 2 to be {g2Truth}")
    
    def test_calculate_gt_median_adjustment_2(self):
        "Test that median adjustment for CNVs works correctly (all alleles are zero)"
        # Arrange
        g1Gt = [[0, 0], [0, 0], [0, 0], [0, 0]]
        g2Gt = [[0, 0], [0, 0], [0, 0], [0, 0]]
        g1Truth = 1 # equal to median == 1
        g2Truth = 1 # equal to median == 1
        
        # Act
        g1AdjGt, g2AdjGt = gt_median_adjustment([g1Gt, g2Gt])
        
        # Assert
        self.assertTrue(all([ allele == g1Truth for gt in g1AdjGt for allele in gt ]), f"Expected adjusted group 1 to be {g1Truth}")
        self.assertTrue(all([ allele == g2Truth for gt in g2AdjGt for allele in gt ]), f"Expected adjusted group 2 to be {g2Truth}")
    
    def test_calculate_gt_median_adjustment_3(self):
        "Test that median adjustment for CNVs works correctly (all alleles are one)"
        # Arrange
        g1Gt = [[1, 1], [1, 1], [1, 1], [1, 1]]
        g2Gt = [[1, 1], [1, 1], [1, 1], [1, 1]]
        g1Truth = 1 # equal to median == 1
        g2Truth = 1 # equal to median == 1
        
        # Act
        g1AdjGt, g2AdjGt = gt_median_adjustment([g1Gt, g2Gt])
        
        # Assert
        self.assertTrue(all([ allele == g1Truth for gt in g1AdjGt for allele in gt ]), f"Expected adjusted group 1 to be {g1Truth}")
        self.assertTrue(all([ allele == g2Truth for gt in g2AdjGt for allele in gt ]), f"Expected adjusted group 2 to be {g2Truth}")
    
    def test_calculate_gt_median_adjustment_4(self):
        "Test that median adjustment for CNVs works correctly (one genotype is different)"
        # Arrange
        g1Gt = [[0, 0], [0, 0], [0, 0], [0, 0]]
        g2Gt = [[0, 0], [0, 0], [0, 0], [0, 1]]
        g1Truth = 1 # equal to median == 1
        g2Truth = 1 # equal to median == 1
        
        # Act
        g1AdjGt, g2AdjGt = gt_median_adjustment([g1Gt, g2Gt])
        
        # Assert
        self.assertTrue(all([ allele == g1Truth for gt in g1AdjGt for allele in gt ]), f"Expected adjusted group 1 to be {g1Truth}")
        self.assertFalse(all([ allele == g2Truth for gt in g2AdjGt for allele in gt ]), f"Expected adjusted group 2 to NOT be all {g2Truth}")
    
    def test_calculate_gt_median_adjustment_with_dots(self):
        "Test that median adjustment handles dot values safely"
        # Arrange
        g1Gt = [[0, 0], ["."], ".", [1, 1]]
        g1Truth = [[0, 0], ['.'], ['.'], [2, 2]] # below median == 0, skip over dots, median is 0.5, so above median == 2
        
        # Act
        g1AdjGt = gt_median_adjustment([g1Gt])[0]
        
        # Assert
        self.assertEqual(g1AdjGt, g1Truth, f"Expected adjusted group 1 to be {g1Truth}")

class TestSPLSDA(unittest.TestCase):
    def test_recode_variant_with_pipes(self):
        "Test that recoding a variant works correctly when pipes separate alleles"
        # Arrange
        gtIndex = 1
        sampleFields = ["example:0|0", "example2:0/0", "example3:1|1:otherstuff"]
        truth = ["0", "0", "2"]
        
        # Act
        recoded_variant = recode_variant(gtIndex, sampleFields)
        
        # Assert
        self.assertEqual(recoded_variant, truth, f"Expected {truth} but got {recoded_variant}")
    
    def test_recode_variant_with_dots(self):
        "Test that recoding a variant works correctly when dots are present"
        # Arrange
        gtIndex = 0
        sampleFields = ["0/0", "0/.", "./."]
        truth = ["0", ".", "."]
        
        # Act
        recoded_variant = recode_variant(gtIndex, sampleFields)
        
        # Assert
        self.assertEqual(recoded_variant, truth, f"Expected {truth} but got {recoded_variant}")
    
    def test_recode_variant_1(self):
        "Test that recoding a variant works correctly with heterozygotes"
        # Arrange
        gtIndex = 0
        sampleFields = ["0/0", "0/1", "0/1"]
        truth = ["0", "1", "1"]
        
        # Act
        recoded_variant = recode_variant(gtIndex, sampleFields)
        
        # Assert
        self.assertEqual(recoded_variant, truth, f"Expected {truth} but got {recoded_variant}")
    
    def test_recode_variant_2(self):
        "Test that recoding a variant works correctly with a homozygote"
        # Arrange
        gtIndex = 0
        sampleFields = ["0/0", "0/1", "1/1"]
        truth = ["0", "1", "2"]
        
        # Act
        recoded_variant = recode_variant(gtIndex, sampleFields)
        
        # Assert
        self.assertEqual(recoded_variant, truth, f"Expected {truth} but got {recoded_variant}")
    
    def test_recode_cnv_with_pipes(self):
        "Test that recoding a CNV works correctly when pipes separate alleles"
        # Arrange
        gtIndex = 1
        sampleFields = ["example:0|0", "example2:0/0", "example3:1|1:otherstuff"]
        truth = ["1", "1", "2"] # median is 0, so equals median == 1, above median == 2
        
        # Act
        recoded_cnv = recode_cnv(gtIndex, sampleFields)
        
        # Assert
        self.assertEqual(recoded_cnv, truth, f"Expected {truth} but got {recoded_cnv}")
    
    def test_recode_cnv_with_dots(self):
        "Test that recoding a CNV works correctly when dots are present"
        # Arrange
        gtIndex = 0
        sampleFields = ["0/0", "0/.", "./."]
        truth = ["1", ".", "."] # median is 0, so equals median == 1
        
        # Act
        recoded_cnv = recode_cnv(gtIndex, sampleFields)
        
        # Assert
        self.assertEqual(recoded_cnv, truth, f"Expected {truth} but got {recoded_cnv}")
    
    def test_recode_cnv_1(self):
        "Test that median adjustment for CNVs works correctly (samples segregate strongly)"
        # Arrange
        gtIndex = 0
        sampleFields = ["0/0", "1/1", "100/100", "100/100"]
        truth = ["0", "0", "2", "2"]
        
        # Act
        recoded_cnv = recode_cnv(gtIndex, sampleFields)
        
        # Assert
        self.assertEqual(recoded_cnv, truth, f"Expected {truth} but got {recoded_cnv}")
    
    def test_recode_cnv_2(self):
        "Test that median adjustment for CNVs works correctly (CNV presence is evenly distributed)"
        # Arrange
        gtIndex = 0
        sampleFields = ["0/0", "0/1", "1/1", "1/2", "2/2"]
        truth = ["0", "1", "1", "2", "2"]
        
        # Act
        recoded_cnv = recode_cnv(gtIndex, sampleFields)
        
        # Assert
        self.assertEqual(recoded_cnv, truth, f"Expected {truth} but got {recoded_cnv}")
    
    def test_recode_vcf(self):
        "Test that recoding a standard, diploid deletion VCF works correctly"
        # Arrange
        metadataDict = parse_metadata(metadataFile)
        
        workDir = os.path.join(dataDir, "tmp")
        outputFile = os.path.join(workDir, "recode.vcf.gz")
        vcfFile = os.path.join(dataDir, "deletions.1.vcf")
        
        numLines = 11
        numHeaderColumns = 54
        
        # Arrange: cleanup any previous work directory
        if os.path.exists(workDir):
            shutil.rmtree(workDir)
        if not os.path.exists(workDir):
            os.makedirs(workDir)
        
        # Act
        recode_vcf(vcfFile, outputFile, metadataDict, isCNV=True, quiet=True)
        recodeContents = []
        with read_gz_file(outputFile) as fileIn:
            for line in fileIn:
                recodeContents.append(line.strip())
        
        # Assert
        self.assertTrue(len(recodeContents) == numLines, 
                        f"Expected recode file to be {numLines} lines long but got {len(recodeContents)} lines")
        lengthOfHeader = len(recodeContents[0].split("\t"))
        self.assertTrue(lengthOfHeader == numHeaderColumns, 
                        f"Expected recode file to have {numHeaderColumns} columns but got {lengthOfHeader} columns")
    
    def test_recode_vcf_with_subset_metadata(self):
        "Test that recoding a standard, diploid deletion VCF works correctly when using a subset of metadata"
        # Arrange
        metadataDict = parse_metadata(metadataFile)
        metadataDict["group1"] = ['S44', 'S48', 'S50', 'S54', 'S56', 'S59', 'S64', 'S69', 'S76', 'S77'] # drop S95
        
        workDir = os.path.join(dataDir, "tmp")
        outputFile = os.path.join(workDir, "recode.vcf.gz")
        vcfFile = os.path.join(dataDir, "deletions.1.vcf")
        
        numLines = 11
        numHeaderColumns = 53 # 1 less than the full metadata
        
        # Arrange: cleanup any previous work directory
        if os.path.exists(workDir):
            shutil.rmtree(workDir)
        if not os.path.exists(workDir):
            os.makedirs(workDir)
        
        # Act
        recode_vcf(vcfFile, outputFile, metadataDict, isCNV=True, quiet=True)
        recodeContents = []
        with read_gz_file(outputFile) as fileIn:
            for line in fileIn:
                recodeContents.append(line.strip())
        
        # Assert
        self.assertTrue(len(recodeContents) == numLines, 
                        f"Expected recode file to be {numLines} lines long but got {len(recodeContents)} lines")
        lengthOfHeader = len(recodeContents[0].split("\t"))
        self.assertTrue(lengthOfHeader == numHeaderColumns, 
                        f"Expected recode file to have {numHeaderColumns} columns but got {lengthOfHeader} columns")
    
    def test_recode_vcf_tetraploid(self):
        "Test that recoding a tetraploid deletion VCF works correctly; note that tetraploid VCFs won't be created by psQTL but may be input by user"
        # Arrange
        metadataDict = parse_metadata(metadataFile)
        
        workDir = os.path.join(dataDir, "tmp")
        outputFile = os.path.join(workDir, "recode.vcf.gz")
        vcfFile = os.path.join(dataDir, "deletions.2.vcf")
        
        numLines = 11
        numHeaderColumns = 54
        
        # Arrange: cleanup any previous work directory
        if os.path.exists(workDir):
            shutil.rmtree(workDir)
        if not os.path.exists(workDir):
            os.makedirs(workDir)
        
        # Act
        recode_vcf(vcfFile, outputFile, metadataDict, isCNV=True, quiet=True)
        recodeContents = []
        with read_gz_file(outputFile) as fileIn:
            for line in fileIn:
                recodeContents.append(line.strip())
        
        # Assert
        self.assertTrue(len(recodeContents) == numLines, 
                        f"Expected recode file to be {numLines} lines long but got {len(recodeContents)} lines")
        lengthOfHeader = len(recodeContents[0].split("\t"))
        self.assertTrue(lengthOfHeader == numHeaderColumns, 
                        f"Expected recode file to have {numHeaderColumns} columns but got {lengthOfHeader} columns")

class TestMain(unittest.TestCase):
    def test_empty_vcf(self):
        "Test the psQTL analysis pipeline with an empty VCF input"
        # Arrange: set variables
        workDir = os.path.join(dataDir, "tmp")
        fulltestMetadata = os.path.join(dataDir, "fulltest.metadata.1.tsv")
        vcfFile = os.path.join(dataDir, "empty.vcf")
        genomeFile = os.path.join(dataDir, "genome.fasta")
        expectedError = "VCF file is empty or has no variants"
        
        # Arrange: cleanup any previous work directory
        if os.path.exists(workDir):
            shutil.rmtree(workDir)
        if not os.path.exists(workDir):
            os.makedirs(workDir)
        
        # Act&Assert: run psQTL_prep.py initialise
        cmd = [
            "python", os.path.join(baseDir, "psQTL_prep.py"), "initialise",
            "-d", workDir,
            "--meta", fulltestMetadata,
            "--fvcf", vcfFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(expectedError in stderr, (f"Expected '{expectedError}' message in " +
                                                  f"stderr output but got: {stderr}"))
    
    def test_diploid_variants_no_parents_full_metadata(self):
        "Test a full psQTL analysis pipeline with a test set of variants and metadata"
        # Arrange: set variables
        workDir = os.path.join(dataDir, "tmp")
        fulltestMetadata = os.path.join(dataDir, "fulltest.metadata.1.tsv")
        vcfFile = os.path.join(dataDir, "fulltest.variants.1.vcf")
        genomeFile = os.path.join(dataDir, "genome.fasta")
        gff3File = os.path.join(dataDir, "fulltest.1.gff3")
        
        edTruth = "chr1\t10\tsnp\t42\t42\t98\t98\t1.0678755470980512" # 1.0967370483709717 if genotype inheritance is used
        berTruth = "chr1\t0\t0.0714285714285715\n"
        selectedTruth = "chr1\t10\t1\t1\tleft\n"
        recodeTruth = ['chrom\tpos\tbulk1_1\tbulk1_10\tbulk1_11\tbulk1_12\tbulk1_13\tbulk1_14\tbulk1_15\tbulk1_16\tbulk1_17\tbulk1_18\tbulk1_19\tbulk1_2\tbulk1_20\tbulk1_21\tbulk1_3\tbulk1_4\tbulk1_5\tbulk1_6\tbulk1_7\tbulk1_8\tbulk1_9\tbulk2_1\tbulk2_10\tbulk2_11\tbulk2_12\tbulk2_13\tbulk2_14\tbulk2_15\tbulk2_16\tbulk2_17\tbulk2_18\tbulk2_19\tbulk2_2\tbulk2_20\tbulk2_21\tbulk2_22\tbulk2_23\tbulk2_24\tbulk2_25\tbulk2_26\tbulk2_27\tbulk2_28\tbulk2_29\tbulk2_3\tbulk2_30\tbulk2_31\tbulk2_32\tbulk2_33\tbulk2_34\tbulk2_35\tbulk2_36\tbulk2_37\tbulk2_38\tbulk2_39\tbulk2_4\tbulk2_40\tbulk2_41\tbulk2_42\tbulk2_43\tbulk2_44\tbulk2_45\tbulk2_46\tbulk2_47\tbulk2_48\tbulk2_49\tbulk2_5\tbulk2_6\tbulk2_7\tbulk2_8\tbulk2_9',
                       'chr1\t10\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t1\t1\t2\t2\t2\t1\t2\t2\t1\t0\t2\t1\t1\t2\t1\t1\t0\t1\t0\t1\t0\t0\t0\t0\t0']
        
        # Arrange: cleanup any previous work directory
        if os.path.exists(workDir):
            shutil.rmtree(workDir)
        if not os.path.exists(workDir):
            os.makedirs(workDir)
        
        # Act&Assert: run psQTL_prep.py initialise
        cmd = [
            "python", os.path.join(baseDir, "psQTL_prep.py"), "initialise",
            "-d", workDir,
            "--meta", fulltestMetadata,
            "--fvcf", vcfFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Act&Assert: run psQTL_proc.py ed
        cmd = [
            "python", os.path.join(baseDir, "psQTL_proc.py"), "ed",
            "-d", workDir,
            "-i", "call"
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        edFile = os.path.join(workDir, "psQTL_call.alleles_ed.tsv.gz")
        edContents = []
        with read_gz_file(edFile) as fileIn:
            for line in fileIn:
                edContents.append(line.strip())
        self.assertTrue(edContents[1] == edTruth, f"Expected ED file to contain '{edTruth}' but got: {edContents[1]}")
        
        # Act&Assert: run psQTL_proc.py splsda
        cmd = [
            "python", os.path.join(baseDir, "psQTL_proc.py"), "splsda",
            "-d", workDir,
            "-i", "call"
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Assert: check that the output files are correctly generated
        berFile = os.path.join(workDir, "splsda", "psQTL_call.BER.tsv")
        with open(berFile) as fileIn:
            berContents = fileIn.readlines()
        self.assertTrue(berContents[1] == berTruth, f"Expected BER file to contain '{berTruth}' but got: {berContents[1]}")
        
        selectedFile = os.path.join(workDir, "splsda", "psQTL_call.selected.tsv")
        with open(selectedFile) as fileIn:
            selectedContents = fileIn.readlines()
        self.assertTrue(selectedContents[1] == selectedTruth, f"Expected selected file to contain '{selectedTruth}' but got: {selectedContents[1]}")
        
        recodeFile = os.path.join(workDir, "splsda", "psQTL_call.recode.tsv.gz")
        recodeContents = []
        with read_gz_file(recodeFile) as fileIn:
            for line in fileIn:
                recodeContents.append(line.strip())
        self.assertTrue(recodeContents == recodeTruth, f"Expected recode file to be '{recodeTruth}' but got: {recodeContents}")
    
    def test_diploid_variants_parents_full_metadata(self):
        "Test a full psQTL analysis pipeline with a test set of variants (parents calculation) and metadata"
        # Arrange: set variables
        workDir = os.path.join(dataDir, "tmp")
        fulltestMetadata = os.path.join(dataDir, "fulltest.metadata.1.tsv")
        vcfFile = os.path.join(dataDir, "fulltest.variants.1.vcf")
        genomeFile = os.path.join(dataDir, "genome.fasta")
        
        edTruth = "chr1\t10\tsnp\t0\t42\t20\t98\t0"
        berTruth = "chr1\t0\t0.0714285714285715\n"
        selectedTruth = "chr1\t10\t1\t1\tleft\n"
        recodeTruth = ['chrom\tpos\tbulk1_1\tbulk1_10\tbulk1_11\tbulk1_12\tbulk1_13\tbulk1_14\tbulk1_15\tbulk1_16\tbulk1_17\tbulk1_18\tbulk1_19\tbulk1_2\tbulk1_20\tbulk1_21\tbulk1_3\tbulk1_4\tbulk1_5\tbulk1_6\tbulk1_7\tbulk1_8\tbulk1_9\tbulk2_1\tbulk2_10\tbulk2_11\tbulk2_12\tbulk2_13\tbulk2_14\tbulk2_15\tbulk2_16\tbulk2_17\tbulk2_18\tbulk2_19\tbulk2_2\tbulk2_20\tbulk2_21\tbulk2_22\tbulk2_23\tbulk2_24\tbulk2_25\tbulk2_26\tbulk2_27\tbulk2_28\tbulk2_29\tbulk2_3\tbulk2_30\tbulk2_31\tbulk2_32\tbulk2_33\tbulk2_34\tbulk2_35\tbulk2_36\tbulk2_37\tbulk2_38\tbulk2_39\tbulk2_4\tbulk2_40\tbulk2_41\tbulk2_42\tbulk2_43\tbulk2_44\tbulk2_45\tbulk2_46\tbulk2_47\tbulk2_48\tbulk2_49\tbulk2_5\tbulk2_6\tbulk2_7\tbulk2_8\tbulk2_9',
                       'chr1\t10\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t1\t1\t2\t2\t2\t1\t2\t2\t1\t0\t2\t1\t1\t2\t1\t1\t0\t1\t0\t1\t0\t0\t0\t0\t0']
        
        # Arrange: cleanup any previous work directory
        if os.path.exists(workDir):
            shutil.rmtree(workDir)
        if not os.path.exists(workDir):
            os.makedirs(workDir)
        
        # Act&Assert: run psQTL_prep.py initialise
        cmd = [
            "python", os.path.join(baseDir, "psQTL_prep.py"), "initialise",
            "-d", workDir,
            "--meta", fulltestMetadata,
            "--fvcf", vcfFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        #self.assertTrue(returncode == 0, f"Expected returncode 0 but got: {returncode}")
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Act&Assert: run psQTL_proc.py ed
        cmd = [
            "python", os.path.join(baseDir, "psQTL_proc.py"), "ed",
            "-d", workDir,
            "-i", "call", "--parents", "bulk1_1", "bulk2_29"
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        #self.assertTrue(returncode == 0, f"Expected returncode 0 but got: {returncode}")
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        edFile = os.path.join(workDir, "psQTL_call.inheritance_ed.tsv.gz")
        edContents = []
        with read_gz_file(edFile) as fileIn:
            for line in fileIn:
                edContents.append(line.strip())
        self.assertTrue(edContents[1] == edTruth, f"Expected ED file to contain '{edTruth}' but got: {edContents[1]}")
        
        # Act&Assert: run psQTL_proc.py splsda
        cmd = [
            "python", os.path.join(baseDir, "psQTL_proc.py"), "splsda",
            "-d", workDir,
            "-i", "call"
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        #self.assertTrue(returncode == 0, f"Expected returncode 0 but got: {returncode}")
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Assert: check that the output files are correctly generated
        berFile = os.path.join(workDir, "splsda", "psQTL_call.BER.tsv")
        with open(berFile) as fileIn:
            berContents = fileIn.readlines()
        self.assertTrue(berContents[1] == berTruth, f"Expected BER file to contain '{berTruth}' but got: {berContents[1]}")
        
        selectedFile = os.path.join(workDir, "splsda", "psQTL_call.selected.tsv")
        with open(selectedFile) as fileIn:
            selectedContents = fileIn.readlines()
        self.assertTrue(selectedContents[1] == selectedTruth, f"Expected selected file to contain '{selectedTruth}' but got: {selectedContents[1]}")
        
        recodeFile = os.path.join(workDir, "splsda", "psQTL_call.recode.tsv.gz")
        recodeContents = []
        with read_gz_file(recodeFile) as fileIn:
            for line in fileIn:
                recodeContents.append(line.strip())
        self.assertTrue(recodeContents == recodeTruth, f"Expected recode file to be '{recodeTruth}' but got: {recodeContents}")
    
    def test_diploid_variants_no_parents_subset_metadata(self):
        "Test a full psQTL analysis pipeline with a test set of variants and subsetted metadata"
        # Arrange: set variables
        workDir = os.path.join(dataDir, "tmp")
        fulltestMetadata = os.path.join(dataDir, "subset.metadata.1.tsv")
        vcfFile = os.path.join(dataDir, "fulltest.variants.1.vcf")
        genomeFile = os.path.join(dataDir, "genome.fasta")
        
        edTruth = "chr1\t10\tsnp\t42\t42\t58\t58\t1.4142135623730951"
        berTruth = "chr1\t0\t0\n"
        selectedTruth = "chr1\t10\t1\t1\tleft\n"
        recodeTruth = ['chrom\tpos\tbulk1_1\tbulk1_10\tbulk1_11\tbulk1_12\tbulk1_13\tbulk1_14\tbulk1_15\tbulk1_16\tbulk1_17\tbulk1_18\tbulk1_19\tbulk1_2\tbulk1_20\tbulk1_21\tbulk1_3\tbulk1_4\tbulk1_5\tbulk1_6\tbulk1_7\tbulk1_8\tbulk1_9\tbulk2_1\tbulk2_10\tbulk2_11\tbulk2_12\tbulk2_13\tbulk2_14\tbulk2_15\tbulk2_16\tbulk2_17\tbulk2_18\tbulk2_19\tbulk2_2\tbulk2_20\tbulk2_21\tbulk2_22\tbulk2_23\tbulk2_24\tbulk2_25\tbulk2_26\tbulk2_27\tbulk2_28\tbulk2_29\tbulk2_3\tbulk2_4\tbulk2_5\tbulk2_6\tbulk2_7\tbulk2_8\tbulk2_9',
                       'chr1\t10\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0']
        
        # Arrange: cleanup any previous work directory
        if os.path.exists(workDir):
            shutil.rmtree(workDir)
        if not os.path.exists(workDir):
            os.makedirs(workDir)
        
        # Act&Assert: run psQTL_prep.py initialise
        cmd = [
            "python", os.path.join(baseDir, "psQTL_prep.py"), "initialise",
            "-d", workDir,
            "--meta", fulltestMetadata,
            "--fvcf", vcfFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        #self.assertTrue(returncode == 0, f"Expected returncode 0 but got: {returncode}")
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Act&Assert: run psQTL_proc.py ed
        cmd = [
            "python", os.path.join(baseDir, "psQTL_proc.py"), "ed",
            "-d", workDir,
            "-i", "call"
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        #self.assertTrue(returncode == 0, f"Expected returncode 0 but got: {returncode}")
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        edFile = os.path.join(workDir, "psQTL_call.alleles_ed.tsv.gz")
        edContents = []
        with read_gz_file(edFile) as fileIn:
            for line in fileIn:
                edContents.append(line.strip())
        self.assertTrue(edContents[1] == edTruth, f"Expected ED file to contain '{edTruth}' but got: {edContents[1]}")
        
        # Act&Assert: run psQTL_proc.py splsda
        cmd = [
            "python", os.path.join(baseDir, "psQTL_proc.py"), "splsda",
            "-d", workDir,
            "-i", "call"
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        #self.assertTrue(returncode == 0, f"Expected returncode 0 but got: {returncode}")
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Assert: check that the output files are correctly generated
        berFile = os.path.join(workDir, "splsda", "psQTL_call.BER.tsv")
        with open(berFile) as fileIn:
            berContents = fileIn.readlines()
        self.assertTrue(berContents[1] == berTruth, f"Expected BER file to contain '{berTruth}' but got: {berContents[1]}")
        
        selectedFile = os.path.join(workDir, "splsda", "psQTL_call.selected.tsv")
        with open(selectedFile) as fileIn:
            selectedContents = fileIn.readlines()
        self.assertTrue(selectedContents[1] == selectedTruth, f"Expected selected file to contain '{selectedTruth}' but got: {selectedContents[1]}")
        
        recodeFile = os.path.join(workDir, "splsda", "psQTL_call.recode.tsv.gz")
        recodeContents = []
        with read_gz_file(recodeFile) as fileIn:
            for line in fileIn:
                recodeContents.append(line.strip())
        self.assertTrue(recodeContents == recodeTruth, f"Expected recode file to be '{recodeTruth}' but got: {recodeContents}")
    
    def test_diploid_variants_parents_subset_metadata(self):
        "Test a full psQTL analysis pipeline with a test set of variants (parents calculation) and metadata"
        # Arrange: set variables
        workDir = os.path.join(dataDir, "tmp")
        fulltestMetadata = os.path.join(dataDir, "subset.metadata.1.tsv")
        vcfFile = os.path.join(dataDir, "fulltest.variants.1.vcf")
        genomeFile = os.path.join(dataDir, "genome.fasta")
        
        edTruth = "chr1\t10\tsnp\t0\t42\t0\t58\t0"
        berTruth = "chr1\t0\t0\n"
        selectedTruth = "chr1\t10\t1\t1\tleft\n"
        recodeTruth = ['chrom\tpos\tbulk1_1\tbulk1_10\tbulk1_11\tbulk1_12\tbulk1_13\tbulk1_14\tbulk1_15\tbulk1_16\tbulk1_17\tbulk1_18\tbulk1_19\tbulk1_2\tbulk1_20\tbulk1_21\tbulk1_3\tbulk1_4\tbulk1_5\tbulk1_6\tbulk1_7\tbulk1_8\tbulk1_9\tbulk2_1\tbulk2_10\tbulk2_11\tbulk2_12\tbulk2_13\tbulk2_14\tbulk2_15\tbulk2_16\tbulk2_17\tbulk2_18\tbulk2_19\tbulk2_2\tbulk2_20\tbulk2_21\tbulk2_22\tbulk2_23\tbulk2_24\tbulk2_25\tbulk2_26\tbulk2_27\tbulk2_28\tbulk2_29\tbulk2_3\tbulk2_4\tbulk2_5\tbulk2_6\tbulk2_7\tbulk2_8\tbulk2_9',
                       'chr1\t10\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0']
        
        # Arrange: cleanup any previous work directory
        if os.path.exists(workDir):
            shutil.rmtree(workDir)
        if not os.path.exists(workDir):
            os.makedirs(workDir)
        
        # Act&Assert: run psQTL_prep.py initialise
        cmd = [
            "python", os.path.join(baseDir, "psQTL_prep.py"), "initialise",
            "-d", workDir,
            "--meta", fulltestMetadata,
            "--fvcf", vcfFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        #self.assertTrue(returncode == 0, f"Expected returncode 0 but got: {returncode}")
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Act&Assert: run psQTL_proc.py ed
        cmd = [
            "python", os.path.join(baseDir, "psQTL_proc.py"), "ed",
            "-d", workDir,
            "-i", "call", "--parents", "bulk1_1", "bulk2_29"
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        #self.assertTrue(returncode == 0, f"Expected returncode 0 but got: {returncode}")
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        edFile = os.path.join(workDir, "psQTL_call.inheritance_ed.tsv.gz")
        edContents = []
        with read_gz_file(edFile) as fileIn:
            for line in fileIn:
                edContents.append(line.strip())
        self.assertTrue(edContents[1] == edTruth, f"Expected ED file to contain '{edTruth}' but got: {edContents[1]}")
        
        # Act&Assert: run psQTL_proc.py splsda
        cmd = [
            "python", os.path.join(baseDir, "psQTL_proc.py"), "splsda",
            "-d", workDir,
            "-i", "call"
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        #self.assertTrue(returncode == 0, f"Expected returncode 0 but got: {returncode}")
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Assert: check that the output files are correctly generated
        berFile = os.path.join(workDir, "splsda", "psQTL_call.BER.tsv")
        with open(berFile) as fileIn:
            berContents = fileIn.readlines()
        self.assertTrue(berContents[1] == berTruth, f"Expected BER file to contain '{berTruth}' but got: {berContents[1]}")
        
        selectedFile = os.path.join(workDir, "splsda", "psQTL_call.selected.tsv")
        with open(selectedFile) as fileIn:
            selectedContents = fileIn.readlines()
        self.assertTrue(selectedContents[1] == selectedTruth, f"Expected selected file to contain '{selectedTruth}' but got: {selectedContents[1]}")
        
        recodeFile = os.path.join(workDir, "splsda", "psQTL_call.recode.tsv.gz")
        recodeContents = []
        with read_gz_file(recodeFile) as fileIn:
            for line in fileIn:
                recodeContents.append(line.strip())
        self.assertTrue(recodeContents == recodeTruth, f"Expected recode file to be '{recodeTruth}' but got: {recodeContents}")
    
    def test_tetraploid_variants_no_parents_full_metadata(self):
        "Test a full psQTL analysis pipeline with a test set of tetraploid variants and metadata"
        # Arrange: set variables
        workDir = os.path.join(dataDir, "tmp")
        fulltestMetadata = os.path.join(dataDir, "fulltest.metadata.1.tsv")
        vcfFile = os.path.join(dataDir, "fulltest.variants.2.vcf")
        genomeFile = os.path.join(dataDir, "genome.fasta")
        
        edTruth = "chr1\t10\tsnp\t84\t84\t196\t196\t1.0678755470980512" # twice as many alleles, same ED as diploid version
        berTruth = "chr1\t0\t0.0714285714285715\n"
        selectedTruth = "chr1\t10\t1\t1\tleft\n"
        recodeTruth = ['chrom\tpos\tbulk1_1\tbulk1_10\tbulk1_11\tbulk1_12\tbulk1_13\tbulk1_14\tbulk1_15\tbulk1_16\tbulk1_17\tbulk1_18\tbulk1_19\tbulk1_2\tbulk1_20\tbulk1_21\tbulk1_3\tbulk1_4\tbulk1_5\tbulk1_6\tbulk1_7\tbulk1_8\tbulk1_9\tbulk2_1\tbulk2_10\tbulk2_11\tbulk2_12\tbulk2_13\tbulk2_14\tbulk2_15\tbulk2_16\tbulk2_17\tbulk2_18\tbulk2_19\tbulk2_2\tbulk2_20\tbulk2_21\tbulk2_22\tbulk2_23\tbulk2_24\tbulk2_25\tbulk2_26\tbulk2_27\tbulk2_28\tbulk2_29\tbulk2_3\tbulk2_30\tbulk2_31\tbulk2_32\tbulk2_33\tbulk2_34\tbulk2_35\tbulk2_36\tbulk2_37\tbulk2_38\tbulk2_39\tbulk2_4\tbulk2_40\tbulk2_41\tbulk2_42\tbulk2_43\tbulk2_44\tbulk2_45\tbulk2_46\tbulk2_47\tbulk2_48\tbulk2_49\tbulk2_5\tbulk2_6\tbulk2_7\tbulk2_8\tbulk2_9',
                       'chr1\t10\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t2\t2\t4\t4\t4\t2\t4\t4\t2\t0\t4\t2\t2\t4\t2\t2\t0\t2\t0\t2\t0\t0\t0\t0\t0']
        
        # Arrange: cleanup any previous work directory
        if os.path.exists(workDir):
            shutil.rmtree(workDir)
        if not os.path.exists(workDir):
            os.makedirs(workDir)
        
        # Act&Assert: run psQTL_prep.py initialise
        cmd = [
            "python", os.path.join(baseDir, "psQTL_prep.py"), "initialise",
            "-d", workDir,
            "--meta", fulltestMetadata,
            "--fvcf", vcfFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(returncode == 0, f"Expected returncode 0 but got: {returncode}")
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Act&Assert: run psQTL_proc.py ed
        cmd = [
            "python", os.path.join(baseDir, "psQTL_proc.py"), "ed",
            "-d", workDir,
            "-i", "call"
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(returncode == 0, f"Expected returncode 0 but got: {returncode}")
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        edFile = os.path.join(workDir, "psQTL_call.alleles_ed.tsv.gz")
        edContents = []
        with read_gz_file(edFile) as fileIn:
            for line in fileIn:
                edContents.append(line.strip())
        self.assertTrue(edContents[1] == edTruth, f"Expected ED file to contain '{edTruth}' but got: {edContents[1]}")
        
        # Act&Assert: run psQTL_proc.py splsda
        cmd = [
            "python", os.path.join(baseDir, "psQTL_proc.py"), "splsda",
            "-d", workDir,
            "-i", "call"
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(returncode == 0, f"Expected returncode 0 but got: {returncode}")
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Assert: check that the output files are correctly generated
        berFile = os.path.join(workDir, "splsda", "psQTL_call.BER.tsv")
        with open(berFile) as fileIn:
            berContents = fileIn.readlines()
        self.assertTrue(berContents[1] == berTruth, f"Expected BER file to contain '{berTruth}' but got: {berContents[1]}")
        
        selectedFile = os.path.join(workDir, "splsda", "psQTL_call.selected.tsv")
        with open(selectedFile) as fileIn:
            selectedContents = fileIn.readlines()
        self.assertTrue(selectedContents[1] == selectedTruth, f"Expected selected file to contain '{selectedTruth}' but got: {selectedContents[1]}")
        
        recodeFile = os.path.join(workDir, "splsda", "psQTL_call.recode.tsv.gz")
        recodeContents = []
        with read_gz_file(recodeFile) as fileIn:
            for line in fileIn:
                recodeContents.append(line.strip())
        self.assertTrue(recodeContents == recodeTruth, f"Expected recode file to be '{recodeTruth}' but got: {recodeContents}")
    
    def test_tetraploid_variants_parents_full_metadata(self):
        "Test a full psQTL analysis pipeline with a test set of variants (parents calculation) and metadata"
        # Arrange: set variables
        workDir = os.path.join(dataDir, "tmp")
        fulltestMetadata = os.path.join(dataDir, "fulltest.metadata.1.tsv")
        vcfFile = os.path.join(dataDir, "fulltest.variants.2.vcf")
        genomeFile = os.path.join(dataDir, "genome.fasta")
        
        edTruth = "chr1\t10\tsnp\t0\t84\t40\t196\t0"
        berTruth = "chr1\t0\t0.0714285714285715\n"
        selectedTruth = "chr1\t10\t1\t1\tleft\n"
        recodeTruth = ['chrom\tpos\tbulk1_1\tbulk1_10\tbulk1_11\tbulk1_12\tbulk1_13\tbulk1_14\tbulk1_15\tbulk1_16\tbulk1_17\tbulk1_18\tbulk1_19\tbulk1_2\tbulk1_20\tbulk1_21\tbulk1_3\tbulk1_4\tbulk1_5\tbulk1_6\tbulk1_7\tbulk1_8\tbulk1_9\tbulk2_1\tbulk2_10\tbulk2_11\tbulk2_12\tbulk2_13\tbulk2_14\tbulk2_15\tbulk2_16\tbulk2_17\tbulk2_18\tbulk2_19\tbulk2_2\tbulk2_20\tbulk2_21\tbulk2_22\tbulk2_23\tbulk2_24\tbulk2_25\tbulk2_26\tbulk2_27\tbulk2_28\tbulk2_29\tbulk2_3\tbulk2_30\tbulk2_31\tbulk2_32\tbulk2_33\tbulk2_34\tbulk2_35\tbulk2_36\tbulk2_37\tbulk2_38\tbulk2_39\tbulk2_4\tbulk2_40\tbulk2_41\tbulk2_42\tbulk2_43\tbulk2_44\tbulk2_45\tbulk2_46\tbulk2_47\tbulk2_48\tbulk2_49\tbulk2_5\tbulk2_6\tbulk2_7\tbulk2_8\tbulk2_9',
                       'chr1\t10\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t4\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t2\t2\t4\t4\t4\t2\t4\t4\t2\t0\t4\t2\t2\t4\t2\t2\t0\t2\t0\t2\t0\t0\t0\t0\t0']
        
        # Arrange: cleanup any previous work directory
        if os.path.exists(workDir):
            shutil.rmtree(workDir)
        if not os.path.exists(workDir):
            os.makedirs(workDir)
        
        # Act&Assert: run psQTL_prep.py initialise
        cmd = [
            "python", os.path.join(baseDir, "psQTL_prep.py"), "initialise",
            "-d", workDir,
            "--meta", fulltestMetadata,
            "--fvcf", vcfFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        #self.assertTrue(returncode == 0, f"Expected returncode 0 but got: {returncode}")
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Act&Assert: run psQTL_proc.py ed
        cmd = [
            "python", os.path.join(baseDir, "psQTL_proc.py"), "ed",
            "-d", workDir,
            "-i", "call", "--parents", "bulk1_1", "bulk2_29"
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        #self.assertTrue(returncode == 0, f"Expected returncode 0 but got: {returncode}")
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        edFile = os.path.join(workDir, "psQTL_call.inheritance_ed.tsv.gz")
        edContents = []
        with read_gz_file(edFile) as fileIn:
            for line in fileIn:
                edContents.append(line.strip())
        self.assertTrue(edContents[1] == edTruth, f"Expected ED file to contain '{edTruth}' but got: {edContents[1]}")
        
        # Act&Assert: run psQTL_proc.py splsda
        cmd = [
            "python", os.path.join(baseDir, "psQTL_proc.py"), "splsda",
            "-d", workDir,
            "-i", "call"
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        #self.assertTrue(returncode == 0, f"Expected returncode 0 but got: {returncode}")
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Assert: check that the output files are correctly generated
        berFile = os.path.join(workDir, "splsda", "psQTL_call.BER.tsv")
        with open(berFile) as fileIn:
            berContents = fileIn.readlines()
        self.assertTrue(berContents[1] == berTruth, f"Expected BER file to contain '{berTruth}' but got: {berContents[1]}")
        
        selectedFile = os.path.join(workDir, "splsda", "psQTL_call.selected.tsv")
        with open(selectedFile) as fileIn:
            selectedContents = fileIn.readlines()
        self.assertTrue(selectedContents[1] == selectedTruth, f"Expected selected file to contain '{selectedTruth}' but got: {selectedContents[1]}")
        
        recodeFile = os.path.join(workDir, "splsda", "psQTL_call.recode.tsv.gz")
        recodeContents = []
        with read_gz_file(recodeFile) as fileIn:
            for line in fileIn:
                recodeContents.append(line.strip())
        self.assertTrue(recodeContents == recodeTruth, f"Expected recode file to be '{recodeTruth}' but got: {recodeContents}")
    
    def test_reporting_normal(self):
        "Run a full psQTL analysis pipeline with a test set of variants and metadata, then check the reports generated"
        # Arrange: set variables
        workDir = os.path.join(dataDir, "tmp")
        fulltestMetadata = os.path.join(dataDir, "fulltest.metadata.1.tsv")
        vcfFile = os.path.join(dataDir, "fulltest.variants.1.vcf")
        genomeFile = os.path.join(dataDir, "genome.fasta")
        gff3File = os.path.join(dataDir, "fulltest.1.gff3")
        
        allelesReportFile = os.path.join(workDir, "alleles_report.tsv")
        genotypesReportFile = os.path.join(workDir, "genotypes_report.tsv")
        splsdaReportFile = os.path.join(workDir, "splsda_report.tsv")
        alleleMarkersReportFile = os.path.join(workDir, "alleles_report.markers.tsv")
        genotypesMarkersReportFile = os.path.join(workDir, "genotypes_report.markers.tsv")
        splsdaMarkersReportFile = os.path.join(workDir, "splsda_report.markers.tsv")
        reportRadius = 10
        
        allelesReportTruth = ['gene_id\tcontig\tstrand\tgene_start\tgene_end\tN1\tmax_ED_within_radius\tnum_exon\texon_mean\tnum_intron\tintron_mean\tnum_adjacent\tadjacent_mean',
                              'gene1.1\tchr1\t+\t1\t3\t1\t1.3004\t0\t-1.0000\t0\t-1.0000\t1\t1.3004',
                              'gene2.1\tchr1\t+\t3\t7\t1\t1.3004\t0\t-1.0000\t0\t-1.0000\t1\t1.3004',
                              'gene3.1\tchr1\t+\t8\t11\t1\t1.3004\t1\t1.3004\t0\t-1.0000\t0\t-1.0000']
        genotypesReportTruth = ['gene_id\tcontig\tstrand\tgene_start\tgene_end\tN1\tmax_ED_within_radius\tnum_exon\texon_mean\tnum_intron\tintron_mean\tnum_adjacent\tadjacent_mean',
                                'gene1.1\tchr1\t+\t1\t3\t1\t1.4468\t0\t-1.0000\t0\t-1.0000\t1\t1.4468',
                                'gene2.1\tchr1\t+\t3\t7\t1\t1.4468\t0\t-1.0000\t0\t-1.0000\t1\t1.4468',
                                'gene3.1\tchr1\t+\t8\t11\t1\t1.4468\t1\t1.4468\t0\t-1.0000\t0\t-1.0000']
        splsdaReportTruth = ['gene_id\tcontig\tstrand\tgene_start\tgene_end\tnum_adjacent_call_selected\tnum_overlapping_call_selected\tnum_adjacent_depth_selected\tnum_overlapping_depth_selected\tnum_adjacent_integrated_selected\tnum_overlapping_integrated_selected',
                             'gene1.1\tchr1\t+\t1\t3\t1\t0\t0\t0\t0\t0',
                             'gene2.1\tchr1\t+\t3\t7\t1\t0\t0\t0\t0\t0',
                             'gene3.1\tchr1\t+\t8\t11\t0\t1\t0\t0\t0\t0']
        allelesMarkersReportTruth = ['contig\tposition\tstatistic\tnearest_gene',
                                     'chr1\t10\t1.3004\tgene3.1']
        genotypesMarkersReportTruth = ['contig\tposition\tstatistic\tnearest_gene',
                                       'chr1\t10\t1.4468\tgene3.1']
        splsdaMarkersReportTruth = ['contig\tposition\tstatistic\tnearest_gene',
                                    'chr1\t10\tcall\tgene3.1']
        
        # Arrange: cleanup any previous work directory
        if os.path.exists(workDir):
            shutil.rmtree(workDir)
        if not os.path.exists(workDir):
            os.makedirs(workDir)
        
        # Arrange: run full pipeline to generate the necessary files
        cmd = [
            "python", os.path.join(baseDir, "psQTL_prep.py"), "initialise",
            "-d", workDir,
            "--meta", fulltestMetadata,
            "--fvcf", vcfFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        cmd = [
            "python", os.path.join(baseDir, "psQTL_proc.py"), "ed",
            "-d", workDir,
            "-i", "call"
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        cmd = [
            "python", os.path.join(baseDir, "psQTL_proc.py"), "splsda",
            "-d", workDir,
            "-i", "call"
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Act&Assert: run psQTL_post.py report (alleles ED)
        cmd = [
            "python", os.path.join(baseDir, "psQTL_post.py"), "report",
            "-d", workDir, "-f", genomeFile, "-a", gff3File,
            "-m", "ed-call", "-t", "genes", "--ed", "alleles",
            "--radius", str(reportRadius),
            "-o", allelesReportFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Assert: check that the alleles report file is generated correctly
        reportContents = []
        with read_gz_file(allelesReportFile) as fileIn:
            for line in fileIn:
                reportContents.append(line.strip())
        self.assertTrue(reportContents == allelesReportTruth, f"Expected alleles report file to contain '{allelesReportTruth}' but got: {reportContents}")
        
        # Act&Assert: run psQTL_post.py report (genotypes ED)
        cmd = [
            "python", os.path.join(baseDir, "psQTL_post.py"), "report",
            "-d", workDir, "-f", genomeFile, "-a", gff3File,
            "-m", "ed-call", "-t", "genes", "--ed", "genotypes",
            "--radius", str(reportRadius),
            "-o", genotypesReportFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Assert: check that the alleles report file is generated correctly
        reportContents = []
        with read_gz_file(genotypesReportFile) as fileIn:
            for line in fileIn:
                reportContents.append(line.strip())
        self.assertTrue(reportContents == genotypesReportTruth, f"Expected genotypes report file to contain '{genotypesReportTruth}' but got: {reportContents}")
        
        # Act&Assert: run psQTL_post.py report (sPLS-DA)
        cmd = [
            "python", os.path.join(baseDir, "psQTL_post.py"), "report",
            "-d", workDir, "-f", genomeFile, "-a", gff3File,
            "-m", "splsda", "-t", "genes",
            "--radius", str(reportRadius),
            "-o", splsdaReportFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Assert: check that the alleles report file is generated correctly
        reportContents = []
        with read_gz_file(splsdaReportFile) as fileIn:
            for line in fileIn:
                reportContents.append(line.strip())
        self.assertTrue(reportContents == splsdaReportTruth, f"Expected sPLS-DA report file to contain '{splsdaReportTruth}' but got: {reportContents}")
        
        # Act&Assert: run psQTL_post.py report (allele markers)
        cmd = [
            "python", os.path.join(baseDir, "psQTL_post.py"), "report",
            "-d", workDir, "-f", genomeFile, "-a", gff3File,
            "-m", "ed-call", "-t", "markers",
            "--radius", str(reportRadius),
            "-o", alleleMarkersReportFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Assert: check that the alleles report file is generated correctly
        reportContents = []
        with read_gz_file(alleleMarkersReportFile) as fileIn:
            for line in fileIn:
                reportContents.append(line.strip())
        self.assertTrue(reportContents == allelesMarkersReportTruth, f"Expected allele markers report file to contain '{allelesMarkersReportTruth}' but got: {reportContents}")
        
        # Act&Assert: run psQTL_post.py report (genotype markers)
        cmd = [
            "python", os.path.join(baseDir, "psQTL_post.py"), "report",
            "-d", workDir, "-f", genomeFile, "-a", gff3File,
            "-m", "ed-call", "-t", "markers", "--ed", "genotypes",
            "--radius", str(reportRadius),
            "-o", genotypesMarkersReportFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Assert: check that the alleles report file is generated correctly
        reportContents = []
        with read_gz_file(genotypesMarkersReportFile) as fileIn:
            for line in fileIn:
                reportContents.append(line.strip())
        self.assertTrue(reportContents == genotypesMarkersReportTruth, f"Expected genotypes markers report file to contain '{genotypesMarkersReportTruth}' but got: {reportContents}")
        
        # Act&Assert: run psQTL_post.py report (sPLS-DA markers)
        cmd = [
            "python", os.path.join(baseDir, "psQTL_post.py"), "report",
            "-d", workDir, "-f", genomeFile, "-a", gff3File,
            "-m", "splsda", "-t", "markers",
            "--radius", str(reportRadius),
            "-o", splsdaMarkersReportFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Assert: check that the  report file is generated correctly
        reportContents = []
        with read_gz_file(splsdaMarkersReportFile) as fileIn:
            for line in fileIn:
                reportContents.append(line.strip())
        self.assertTrue(reportContents == splsdaMarkersReportTruth, f"Expected sPLS-DA markers report file to contain '{splsdaMarkersReportTruth}' but got: {reportContents}")
    
    def test_reporting_edge_cases(self):
        "Attempt to run psQTL_post.py report with potential edge cases"
        # Arrange: set variables
        workDir = os.path.join(dataDir, "tmp")
        fulltestMetadata = os.path.join(dataDir, "fulltest.metadata.1.tsv")
        vcfFile = os.path.join(dataDir, "fulltest.variants.1.vcf")
        genomeFile = os.path.join(dataDir, "genome.fasta")
        reportFile = os.path.join(workDir, "report.tsv")
        noneWithinGff3 = os.path.join(dataDir, "fulltest.2.gff3")
        smallRadius = 1
        reportRadius = 10
        
        smallRadiusMarkerTruth = ['contig\tposition\tstatistic\tnearest_gene', 'chr1\t10\t1.3004\t.']
        smallRadiusSplsdaTruth = ['contig\tposition\tstatistic\tnearest_gene', 'chr1\t10\tcall\t.']
        
        # Arrange: cleanup any previous work directory
        if os.path.exists(workDir):
            shutil.rmtree(workDir)
        if not os.path.exists(workDir):
            os.makedirs(workDir)
        
        # Arrange: run full pipeline to generate the necessary files
        cmd = [
            "python", os.path.join(baseDir, "psQTL_prep.py"), "initialise",
            "-d", workDir,
            "--meta", fulltestMetadata,
            "--fvcf", vcfFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        cmd = [
            "python", os.path.join(baseDir, "psQTL_proc.py"), "ed",
            "-d", workDir,
            "-i", "call"
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        cmd = [
            "python", os.path.join(baseDir, "psQTL_proc.py"), "splsda",
            "-d", workDir,
            "-i", "call"
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Act&Assert: run psQTL_post.py report (alleles ED) with no overlaps
        cmd = [
            "python", os.path.join(baseDir, "psQTL_post.py"), "report",
            "-d", workDir, "-f", genomeFile, "-a", noneWithinGff3,
            "-m", "ed-call", "-t", "genes", "--ed", "alleles",
            "--radius", str(reportRadius),
            "-o", reportFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Assert: check that the report file is generated correctly (only adjacent features found)
        reportContents = []
        with read_gz_file(reportFile) as fileIn:
            for line in fileIn:
                reportContents.append(line.strip())
        self.assertTrue(reportContents[1].endswith("1\t1.3004"), f"Expected no overlaps report file to end with '1\t1.3004' but got: {reportContents[1]}")
        self.assertTrue(reportContents[2].endswith("1\t1.3004"), f"Expected no overlaps report file to end with '1\t1.3004' but got: {reportContents[2]}")
        
        # Act&Assert: run psQTL_post.py report (alleles ED) with a small radius
        os.remove(reportFile) # remove the previous report file
        cmd = [
            "python", os.path.join(baseDir, "psQTL_post.py"), "report",
            "-d", workDir, "-f", genomeFile, "-a", noneWithinGff3,
            "-m", "ed-call", "-t", "genes", "--ed", "alleles",
            "--radius", str(smallRadius),
            "-o", reportFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Assert: check that the report file is generated correctly (radius too small)
        reportContents = []
        with read_gz_file(reportFile) as fileIn:
            for line in fileIn:
                reportContents.append(line.strip())
        self.assertTrue(len(reportContents) == 1, f"Expected report file with small radius to just be a header line but got: {reportContents}")
        
        # Act&Assert: run psQTL_post.py report (sPLS-DA) with no overlaps
        os.remove(reportFile)
        cmd = [
            "python", os.path.join(baseDir, "psQTL_post.py"), "report",
            "-d", workDir, "-f", genomeFile, "-a", noneWithinGff3,
            "-m", "splsda", "-t", "genes",
            "--radius", str(reportRadius),
            "-o", reportFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Assert: check that the report file is generated correctly (only adjacent features found)
        reportContents = []
        with read_gz_file(reportFile) as fileIn:
            for line in fileIn:
                reportContents.append(line.strip())
        self.assertTrue("\t1\t0\t0\t0\t0\t0" in reportContents[1], f"Expected no overlaps report file to contain '\t1\t0\t0\t0\t0\t0' but got: {reportContents[1]}")
        self.assertTrue("\t1\t0\t0\t0\t0\t0" in reportContents[2], f"Expected no overlaps report file to contain '\t1\t0\t0\t0\t0\t0' but got: {reportContents[2]}")
        
        # Act&Assert: run psQTL_post.py report (sPLS-DA) with a small radius
        os.remove(reportFile)
        cmd = [
            "python", os.path.join(baseDir, "psQTL_post.py"), "report",
            "-d", workDir, "-f", genomeFile, "-a", noneWithinGff3,
            "-m", "splsda", "-t", "genes",
            "--radius", str(smallRadius),
            "-o", reportFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Assert: check that the report file is generated correctly (radius too small)
        reportContents = []
        with read_gz_file(reportFile) as fileIn:
            for line in fileIn:
                reportContents.append(line.strip())
        self.assertTrue(len(reportContents) == 1, f"Expected report file with small radius to just be a header line but got: {reportContents}")
        
        # Act&Assert: run psQTL_post.py report (allele markers) with a small radius
        os.remove(reportFile)
        cmd = [
            "python", os.path.join(baseDir, "psQTL_post.py"), "report",
            "-d", workDir, "-f", genomeFile, "-a", noneWithinGff3,
            "-m", "ed-call", "-t", "markers",
            "--radius", str(smallRadius),
            "-o", reportFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Assert: check that the report file is generated correctly (radius too small)
        reportContents = []
        with read_gz_file(reportFile) as fileIn:
            for line in fileIn:
                reportContents.append(line.strip())
        self.assertTrue(reportContents == smallRadiusMarkerTruth, f"Expected markers report file to contain '{smallRadiusMarkerTruth}' but got: {reportContents}")
        
        # Act&Assert: run psQTL_post.py report (sPLS-DA markers) with a small radius
        os.remove(reportFile)
        cmd = [
            "python", os.path.join(baseDir, "psQTL_post.py"), "report",
            "-d", workDir, "-f", genomeFile, "-a", noneWithinGff3,
            "-m", "splsda", "-t", "markers",
            "--radius", str(smallRadius),
            "-o", reportFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue(stderr == "", f"Expected no stderr output but got: {stderr}")
        
        # Assert: check that the report file is generated correctly (radius too small)
        reportContents = []
        with read_gz_file(reportFile) as fileIn:
            for line in fileIn:
                reportContents.append(line.strip())
        self.assertTrue(reportContents == smallRadiusSplsdaTruth, f"Expected markers report file to contain '{smallRadiusSplsdaTruth}' but got: {reportContents}")
        
class TestGFF3(unittest.TestCase):
    def test_normal_gff3_handling_1(self):
        "Test parsing GFF3 which is fully specified but has some minor formatting issues"
        # Arrange
        gff3File = os.path.join(dataDir, "normal.1.gff3")
        numGenes = 2
        numMrnas = 2
        numExons = 10
        numCDS = 10
        
        # Act
        gff3Obj = GFF3Graph(gff3File)
        
        # Assert
        self.assertTrue(len(gff3Obj.ftypes["gene"]) == numGenes,
                        f"Expected GFF3Graph object to have {numGenes} genes but got {len(gff3Obj.ftypes['gene'])}")
        self.assertTrue(len(gff3Obj.ftypes["mRNA"]) == numMrnas,
                        f"Expected GFF3Graph object to have {numMrnas} mRNAs but got {len(gff3Obj.ftypes['mRNA'])}")
        self.assertTrue(len(gff3Obj.ftypes["exon"]) == numExons,
                        f"Expected GFF3Graph object to have {numExons} exons but got {len(gff3Obj.ftypes['exon'])}")
        self.assertTrue(len(gff3Obj.ftypes["CDS"]) == numCDS,
                        f"Expected GFF3Graph object to have {numCDS} CDS but got {len(gff3Obj.ftypes['CDS'])}")
    
    def test_gff3_longest_isoform(self):
        "Test parsing GFF3 which is fully specified but has some minor formatting issues"
        # Arrange
        gff3File = os.path.join(dataDir, "normal.2.gff3")
        longestIsoformTruth = "Soltu.DM.S001600.1" # 939 bp versus 671 for the other mRNA
        
        # Act
        gff3Obj = GFF3Graph(gff3File)
        longestMrnaFeature = GFF3Graph.longest_isoform(gff3Obj["Soltu.DM.S001600"])
        
        # Assert
        self.assertTrue(longestMrnaFeature.ID == longestIsoformTruth,
                        f"Expected longest mRNA feature to be {longestIsoformTruth} genes but got {longestMrnaFeature.ID}")
    
    def test_faulty_gff3_handling_1(self):
        "Test parsing GFF3 which lacks gene lines"
        # Arrange
        gff3File = os.path.join(dataDir, "faulty.1.gff3")
        numGenes = 2
        numMrnas = 2
        numExons = 10
        numCDS = 10
        
        # Act
        gff3Obj = GFF3Graph(gff3File)
        
        # Assert
        self.assertTrue(len(gff3Obj.ftypes["gene"]) == numGenes,
                        f"Expected GFF3Graph object to have {numGenes} genes but got {len(gff3Obj.ftypes['gene'])}")
        self.assertTrue(len(gff3Obj.ftypes["mRNA"]) == numMrnas,
                        f"Expected GFF3Graph object to have {numMrnas} mRNAs but got {len(gff3Obj.ftypes['mRNA'])}")
        self.assertTrue(len(gff3Obj.ftypes["exon"]) == numExons,
                        f"Expected GFF3Graph object to have {numExons} exons but got {len(gff3Obj.ftypes['exon'])}")
        self.assertTrue(len(gff3Obj.ftypes["CDS"]) == numCDS,
                        f"Expected GFF3Graph object to have {numCDS} CDS but got {len(gff3Obj.ftypes['CDS'])}")
    
    def test_faulty_gff3_handling_2(self):
        "Test parsing GFF3 which has gene lines but no mRNA lines"
        # Arrange
        gff3File = os.path.join(dataDir, "faulty.2.gff3")
        numGenes = 2
        numExons = 10
        numCDS = 10
        
        # Act
        gff3Obj = GFF3Graph(gff3File)
        
        # Assert
        self.assertTrue(len(gff3Obj.ftypes["gene"]) == numGenes,
                        f"Expected GFF3Graph object to have {numGenes} genes but got {len(gff3Obj.ftypes['gene'])}")
        self.assertFalse(hasattr(gff3Obj, "mRNA"),
                        f"Expected GFF3Graph object to have no mRNAs")
        self.assertTrue(len(gff3Obj.ftypes["exon"]) == numExons,
                        f"Expected GFF3Graph object to have {numExons} exons but got {len(gff3Obj.ftypes['exon'])}")
        self.assertTrue(len(gff3Obj.ftypes["CDS"]) == numCDS,
                        f"Expected GFF3Graph object to have {numCDS} CDS but got {len(gff3Obj.ftypes['CDS'])}")
    
    def test_faulty_gff3_handling_3(self):
        "Test parsing GFF3 which has a missing link between gene and exon (mRNA missing)"
        # Arrange
        gff3File = os.path.join(dataDir, "faulty.3.gff3")
        numGenes = 2
        numMrnas = 2
        numExons = 10
        numCDS = 10
        
        mrnaID1 = "Soltu.DM.S001600.1"
        mrnaTruth1 = [6135, 8473]
        
        mrnaID2 = "Soltu.DM.S001610.1"
        mrnaTruth2 = [9694, 11432]
        
        # Act
        gff3Obj = GFF3Graph(gff3File)
        
        # Assert
        self.assertTrue(len(gff3Obj.ftypes["gene"]) == numGenes,
                        f"Expected GFF3Graph object to have {numGenes} genes but got {len(gff3Obj.ftypes['gene'])}")
        self.assertTrue(len(gff3Obj.ftypes["mRNA"]) == numMrnas,
                        f"Expected GFF3Graph object to have {numMrnas} mRNAs but got {len(gff3Obj.ftypes['mRNA'])}")
        self.assertTrue(len(gff3Obj.ftypes["exon"]) == numExons,
                        f"Expected GFF3Graph object to have {numExons} exons but got {len(gff3Obj.ftypes['exon'])}")
        self.assertTrue(len(gff3Obj.ftypes["CDS"]) == numCDS,
                        f"Expected GFF3Graph object to have {numCDS} CDS but got {len(gff3Obj.ftypes['CDS'])}")
        self.assertTrue(all([ len(gff3Obj[x].parents) == 0 for x in gff3Obj.ftypes["mRNA"] ]),
                        "Expected mRNAs to be unlinked from their parents, but some have parents")
        self.assertTrue(gff3Obj[mrnaID1].start == mrnaTruth1[0] and gff3Obj[mrnaID1].end == mrnaTruth1[1],
                        f"Expected mRNA {mrnaID1} to have start {mrnaTruth1[0]} and end {mrnaTruth1[1]}, but got start {gff3Obj[mrnaID1].start} and end {gff3Obj[mrnaID1].end}")
        self.assertTrue(gff3Obj[mrnaID2].start == mrnaTruth2[0] and gff3Obj[mrnaID2].end == mrnaTruth2[1],
                        f"Expected mRNA {mrnaID2} to have start {mrnaTruth2[0]} and end {mrnaTruth2[1]}, but got start {gff3Obj[mrnaID2].start} and end {gff3Obj[mrnaID2].end}")
    
    def test_faulty_gff3_handling_4(self):
        "Test parsing GFF3 which has exon/CDS lines (with no gene or mRNA lines)"
        # Arrange
        gff3File = os.path.join(dataDir, "faulty.4.gff3")
        numMrnas = 2
        numExons = 10
        numCDS = 10
        
        # Act
        gff3Obj = GFF3Graph(gff3File)
        
        # Assert
        self.assertTrue(len(gff3Obj.ftypes["mRNA"]) == numMrnas,
                        f"Expected GFF3Graph object to have {numMrnas} mRNAs but got {len(gff3Obj.ftypes['mRNA'])}")
        self.assertTrue(len(gff3Obj.ftypes["exon"]) == numExons,
                        f"Expected GFF3Graph object to have {numExons} exons but got {len(gff3Obj.ftypes['exon'])}")
        self.assertTrue(len(gff3Obj.ftypes["CDS"]) == numCDS,
                        f"Expected GFF3Graph object to have {numCDS} CDS but got {len(gff3Obj.ftypes['CDS'])}")
    
    def test_faulty_gff3_handling_5(self):
        "Test parsing GFF3 which has exon/CDS lines preceeding mRNA lines (with no gene lines)"
        # Arrange
        gff3File = os.path.join(dataDir, "faulty.5.gff3")
        numGenes = 2 # genes should be inferred from mRNAs without duplication
        numMrnas = 2
        numExons = 10
        numCDS = 10
        
        # Act
        gff3Obj = GFF3Graph(gff3File)
        
        # Assert
        self.assertTrue(len(gff3Obj.ftypes["gene"]) == numGenes,
                        f"Expected GFF3Graph object to have {numGenes} genes but got {len(gff3Obj.ftypes['gene'])}")
        self.assertTrue(len(gff3Obj.ftypes["mRNA"]) == numMrnas,
                        f"Expected GFF3Graph object to have {numMrnas} mRNAs but got {len(gff3Obj.ftypes['mRNA'])}")
        self.assertTrue(len(gff3Obj.ftypes["exon"]) == numExons,
                        f"Expected GFF3Graph object to have {numExons} exons but got {len(gff3Obj.ftypes['exon'])}")
        self.assertTrue(len(gff3Obj.ftypes["CDS"]) == numCDS,
                        f"Expected GFF3Graph object to have {numCDS} CDS but got {len(gff3Obj.ftypes['CDS'])}")
    
    def test_faulty_gff3_handling_6(self):
        "Test parsing GFF3 which has exon/CDS lines preceeding mRNA lines (with no gene lines)"
        # Arrange
        gff3File = os.path.join(dataDir, "faulty.6.gff3")
        numGenes = 2 # neither genes nor mRNAs should have duplicates being inferred
        numMrnas = 2
        numExons = 10
        numCDS = 10
        
        # Act
        gff3Obj = GFF3Graph(gff3File)
        
        # Assert
        self.assertTrue(len(gff3Obj.ftypes["gene"]) == numGenes,
                        f"Expected GFF3Graph object to have {numGenes} genes but got {len(gff3Obj.ftypes['gene'])}")
        self.assertTrue(len(gff3Obj.ftypes["mRNA"]) == numMrnas,
                        f"Expected GFF3Graph object to have {numMrnas} mRNAs but got {len(gff3Obj.ftypes['mRNA'])}")
        self.assertTrue(len(gff3Obj.ftypes["exon"]) == numExons,
                        f"Expected GFF3Graph object to have {numExons} exons but got {len(gff3Obj.ftypes['exon'])}")
        self.assertTrue(len(gff3Obj.ftypes["CDS"]) == numCDS,
                        f"Expected GFF3Graph object to have {numCDS} CDS but got {len(gff3Obj.ftypes['CDS'])}")

    def test_faulty_gff3_handling_7(self):
        "Test parsing GFF3 which has exon/CDS lines preceeding gene lines (with no mRNA lines)"
        # Arrange
        gff3File = os.path.join(dataDir, "faulty.7.gff3")
        numGenes = 2 # neither genes nor mRNAs should have duplicates being inferred
        numMrnas = 2
        numExons = 10
        numCDS = 10
        
        # Act
        gff3Obj = GFF3Graph(gff3File)
        
        # Assert
        self.assertTrue(len(gff3Obj.ftypes["gene"]) == numGenes,
                        f"Expected GFF3Graph object to have {numGenes} genes but got {len(gff3Obj.ftypes['gene'])}")
        self.assertTrue(len(gff3Obj.ftypes["mRNA"]) == numMrnas,
                        f"Expected GFF3Graph object to have {numMrnas} mRNAs but got {len(gff3Obj.ftypes['mRNA'])}")
        self.assertTrue(len(gff3Obj.ftypes["exon"]) == numExons,
                        f"Expected GFF3Graph object to have {numExons} exons but got {len(gff3Obj.ftypes['exon'])}")
        self.assertTrue(len(gff3Obj.ftypes["CDS"]) == numCDS,
                        f"Expected GFF3Graph object to have {numCDS} CDS but got {len(gff3Obj.ftypes['CDS'])}")
        self.assertTrue(all([ len(gff3Obj[x].parents) == 0 for x in gff3Obj.ftypes["mRNA"] ]),
                        "Expected mRNAs to be unlinked from their gene parents, but some have parents")

if __name__ == '__main__':
    unittest.main()
