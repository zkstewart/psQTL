#! python3

import os, sys, unittest, time, math
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.parsing import parse_metadata, vcf_header_to_metadata_validation, parse_vcf_genotypes, \
    parse_vcf_stats, parse_samtools_depth_tsv, parse_binned_tsv
from modules.ncls import WindowedNCLS
from modules.ed import parse_vcf_for_ed, calculate_segregant_ed, calculate_inheritance_ed
from modules.depth import get_median_value, predict_deletions
from modules.samtools_handling import depth_to_histoDict

# Specify data locations
dataDir = os.path.join(os.getcwd(), "data")
metadataFile = os.path.join(dataDir, "metadata.tsv")

MAXIMAL_SEGREGATION = 1.4142135623730951 # sqrt(2)

# Define unit tests
class TestParsing(unittest.TestCase):
    def test_parse_metadata(self):
        "Test parsing metadata file as many subsequent tests depend on it"
        # Arrange & Act
        metadataDict = parse_metadata(metadataFile)
        
        # Assert
        self.assertIsInstance(metadataDict, dict, "Metadata should be a dictionary")
        self.assertIn("bulk1", metadataDict, "Metadata should contain 'bulk1'")
        self.assertIn("bulk2", metadataDict, "Metadata should contain 'bulk2'")
        self.assertEqual(len(metadataDict), 2, "Should have 2 bulks in metadata")
        self.assertEqual(len(metadataDict["bulk1"]), 11, "Bulk1 should have 11 entries")
        self.assertEqual(len(metadataDict["bulk2"]), 41, "Bulk2 should have 41 entries")
    
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
    
    def test_predict_deletions(self):
        "Test that median-normalisation works correctly for deletions"
        # Arrange
        binnedTsvFile = os.path.join(dataDir, "depth.binned.1.tsv")
        histoDict = parse_binned_tsv(binnedTsvFile)
        truth = np.array([0., 0., 0., 0., 1., 1., 2., 2., 2., 2., 2., 2., 3., 3., 4., 4., 4., 4.])
        
        # Act
        alleles = predict_deletions(histoDict["C.glau_01"])
        
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

class TestED(unittest.TestCase):
    def test_parse_vcf_for_ed(self):
        "Test parsing a VCF with deletions where isCNV is True and False"
        # Arrange
        metadataDict = parse_metadata(metadataFile)
        vcfFile = os.path.join(dataDir, "deletions.1.vcf")
        notCNV_truth = [0.478198599250326,0.5942205544782738,0.6083495165377463,
                        0.7877780740982955,0.6366455643787936,0.3733361028931871,
                        0.14110778338534205,0.396669657738795]
        isCNV_truth = [0.478198599250326,0.40607684329781774,0.6914281385881762,
                       0.888979035327655,0.7651177588005215,0.14110778338534202,
                       0.0,0.396669657738795]
        
        # Act & Assert
        notCNV_results = []
        for contig, pos, variant, numAllelesB1, numAllelesB2, euclideanDist \
        in parse_vcf_for_ed(vcfFile, metadataDict, isCNV=False, ignoreIdentical=True, quiet=True):
            notCNV_results.append(euclideanDist)
        
        isCNV_results = []
        for contig, pos, variant, numAllelesB1, numAllelesB2, euclideanDist \
        in parse_vcf_for_ed(vcfFile, metadataDict, isCNV=True, ignoreIdentical=True, quiet=True):
            isCNV_results.append(euclideanDist)
        
        # Assert
        for v1, v2 in zip(notCNV_results, notCNV_truth):
            self.assertAlmostEqual(v1, v2, places=5, msg=f"Expected {v2} but got {v1}")
        for v1, v2 in zip(isCNV_results, isCNV_truth):
            self.assertAlmostEqual(v1, v2, places=5, msg=f"Expected {v2} but got {v1}")
    
    def test_calculate_segregant_ed_1(self):
        "Test that segregation ED is different for isCNV True and False under certain conditions"
        # Arrange
        b1Gt = [[2, 2],[0, 1],[0, 0],[1, 2],[0, 1],[0, 0],[0, 1],[0, 0],[0, 0],[1, 2],[0, 0]]
        b2Gt = [[0, 0],[2, 2],[0, 0],[0, 0],[0, 0],[0, 0],[0, 0],[0, 0],[0, 0],[0, 0],[0, 1],
                [1, 2],[0, 0],[2, 2],[1, 2],[0, 0],[1, 2],[2, 2],[0, 0],[1, 2],[0, 0],[1, 2],
                [0, 0],[0, 0],[0, 0],[0, 0],[0, 1],[1, 2],[1, 2],[2, 2],[1, 2],[1, 2],[2, 2],
                [0, 0],[0, 0],[0, 0],[0, 0],[0, 0],[0, 0],[2, 2],[1, 2]]
        
        # Act
        b1Alleles_1, b2Alleles_1, ed_1 = calculate_segregant_ed(b1Gt, b2Gt, isCNV=True)
        b1Alleles_2, b2Alleles_2, ed_2 = calculate_segregant_ed(b1Gt, b2Gt, isCNV=False)
        
        # Assert
        self.assertNotEqual(ed_1, ed_2, "ED should not be the same for isCNV True and False")
    
    def test_calculate_segregant_ed_2(self):
        "Test that segregation ED is the same for isCNV True and False under certain conditions"
        # Arrange
        b1Gt = [[0, 1],[0, 1],[0, 1],[0, 1],[0, 1],[0, 1],[0, 1],[0, 1],[0, 1],[1, 1],[0, 1]]
        b2Gt = [[0, 0],[0, 1],[0, 1],[0, 1],[0, 0],[0, 0],[0, 1],[0, 1],[0, 1],[0, 0],[0, 1],
                [0, 1],[0, 0],[0, 1],[0, 1],[0, 0],[0, 0],[0, 0],[0, 1],[0, 1],[0, 1],[0, 0],
                [0, 0],[0, 0],[0, 0],[0, 0],[0, 1],[0, 0],[0, 1],[0, 0],[0, 0],[0, 0],[0, 0],
                [0, 1],[0, 0],[0, 1],[0, 0],[0, 0],[0, 0],[0, 0],[0, 0]]
        
        # Act
        b1Alleles_1, b2Alleles_1, ed_1 = calculate_segregant_ed(b1Gt, b2Gt, isCNV=True)
        b1Alleles_2, b2Alleles_2, ed_2 = calculate_segregant_ed(b1Gt, b2Gt, isCNV=False)
        
        # Assert
        self.assertEqual(ed_1, ed_2, "ED should be the same for isCNV True and False")
    
    def test_calculate_inheritance_and_segregant_equality_1(self):
        "Both ED methods should be equal for zero segregation"
        # Arrange
        parentsGT = [ [0, 1], [0, 1] ]
        b1Gt = [[0, 1], [0, 1], [0, 1], [0, 1]]
        b2Gt = [[0, 1], [0, 1], [0, 1], [0, 1]]
        truth = 0
        
        # Act
        b1Alleles_1, b2Alleles_1, ed_1 = calculate_segregant_ed(b1Gt, b2Gt, isCNV=False)
        b1Alleles_2, b2Alleles_2, ed_2 = calculate_inheritance_ed(b1Gt, b2Gt, parentsGT)
        
        # Assert
        self.assertEqual(ed_1, truth, "Segregant ED should be 0 for identical genotypes")
        self.assertEqual(ed_2, truth, "Inheritance ED should be 0 for identical genotypes")
        self.assertEqual(ed_1, ed_2, "Both ED methods should be equal for zero segregation")
    
    def test_calculate_inheritance_and_segregant_equality_2(self):
        "Both ED methods should be equal for maximal segregation"
        # Arrange
        parentsGT = [ [0, 1], [0, 1] ]
        b1Gt = [[0, 0], [0, 0], [0, 0], [0, 0]]
        b2Gt = [[1, 1], [1, 1], [1, 1], [1, 1]]
        truth = MAXIMAL_SEGREGATION
        
        # Act
        b1Alleles_1, b2Alleles_1, ed_1 = calculate_segregant_ed(b1Gt, b2Gt, isCNV=False)
        b1Alleles_2, b2Alleles_2, ed_2 = calculate_inheritance_ed(b1Gt, b2Gt, parentsGT)
        
        # Assert
        self.assertAlmostEqual(ed_1, truth, places=5, msg="Segregant ED should be ~1.41421 for identical genotypes")
        self.assertAlmostEqual(ed_2, truth, places=5, msg="Inheritance ED should be ~1.41421 for identical genotypes")
        self.assertEqual(ed_1, ed_2, "Both ED methods should be equal for maximal segregation")
    
    def test_calculate_inheritance_ed_dip_dip_parents(self):
        # Arrange
        parentsGT = [ [0, 1], [0, 2] ] # represents "A/T" and "A/G" parents
        b1Gt = [[0, 2], [1, 2], [0, 2], [0, 0]] # represents "A/G", "T/G", "A/G", "A/A" genotypes in bulk 1
        b2Gt = [[1, 2], [0, 1], [0, 0], [0, 0]] # represents "T/G", "A/T", "A/A", "A/A" genotypes in bulk 2
        # truth = 0.7905694150420949 # calculated by hand on the plane from WA to Brisbane!
        truth = 0.5590169943749475 # after adjustment of dividing sum by 2 prior to sqrt() operation
        
        # Act
        b1Alleles, b2Alleles, ed = calculate_inheritance_ed(b1Gt, b2Gt, parentsGT)
        
        # Assert
        self.assertAlmostEqual(ed, truth, places=5, 
                             msg=f"Expected ED to be approximately {truth} but got {ed}")
        self.assertEqual(b1Alleles, 8, "Expected 8 alleles in bulk 1")
        self.assertEqual(b2Alleles, 8, "Expected 8 alleles in bulk 2")
    
    def test_calculate_inheritance_ed_tet_tet_parents_1(self):
        "Non-segregation should still give ED==0 with tetraploid parents and samples"
        # Arrange
        parentsGT = [ [0, 0, 1, 1], [0, 0, 2, 2] ]
        b1Gt = [[0, 0, 1, 2], [0, 0, 1, 2], [0, 0, 1, 2], [0, 0, 1, 2]]
        b2Gt = [[0, 0, 1, 2], [0, 0, 1, 2], [0, 0, 1, 2], [0, 0, 1, 2]]
        truth = 0
        
        # Act
        b1Alleles, b2Alleles, ed = calculate_inheritance_ed(b1Gt, b2Gt, parentsGT)
        
        # Assert
        self.assertEqual(ed, truth, f"Expected ED to be zero but got {ed}")
    
    def test_calculate_inheritance_ed_tet_tet_parents_2(self):
        "Maximal-segregation should still give ED==1.41421 with tetraploid parents and samples"
        # Arrange
        parentsGT = [ [0, 1, 2, 3], [0, 1, 4, 5] ]
        b1Gt = [[0, 1, 2, 4], [0, 1, 2, 4], [0, 1, 2, 4], [0, 1, 2, 4]]
        b2Gt = [[0, 1, 3, 5], [0, 1, 3, 5], [0, 1, 3, 5], [0, 1, 3, 5]]
        truth = MAXIMAL_SEGREGATION
        
        # Act
        b1Alleles, b2Alleles, ed = calculate_inheritance_ed(b1Gt, b2Gt, parentsGT)
        
        # Assert
        self.assertEqual(ed, truth, f"Expected ED to be zero but got {ed}")
    
    def test_calculate_inheritance_ed_tet_tet_parents_3(self):
        "Half-shuffled samples from the above test should give ED==0 with tetraploid parents and samples"
        # Arrange
        parentsGT = [ [0, 1, 2, 3], [0, 1, 4, 5] ]
        b1Gt = [[0, 1, 2, 4], [0, 1, 2, 4], [0, 1, 3, 5], [0, 1, 3, 5]]
        b2Gt = [[0, 1, 3, 5], [0, 1, 3, 5], [0, 1, 2, 4], [0, 1, 2, 4]]
        truth = 0
        
        # Act
        b1Alleles, b2Alleles, ed = calculate_inheritance_ed(b1Gt, b2Gt, parentsGT)
        
        # Assert
        self.assertEqual(ed, truth, f"Expected ED to be zero but got {ed}")
    
    def test_calculate_inheritance_ed_with_impossible_progeny_1(self):
        "Test for impossible progeny (alleles do not exist in parents)"
        # Arrange
        parentsGT = [ [0, 0, 1, 1], [0, 0, 2, 2] ]
        b1Gt = [[0, 1, 2, 4], [0, 1, 2, 4], [0, 1, 3, 5], [0, 1, 3, 5]]
        b2Gt = [[0, 1, 3, 5], [0, 1, 3, 5], [0, 1, 2, 4], [0, 1, 2, 4]]
        truth = 0
        
        # Act
        b1Alleles, b2Alleles, ed = calculate_inheritance_ed(b1Gt, b2Gt, parentsGT)
        
        # Assert
        self.assertEqual(ed, truth, f"Expected ED to be zero but got {ed}")
        self.assertEqual(b1Alleles, 0, "Expected 0 alleles in bulk 1")
        self.assertEqual(b2Alleles, 0, "Expected 0 alleles in bulk 2")
    
    def test_calculate_inheritance_ed_with_impossible_progeny_2(self):
        "Test for impossible progeny (allele combination cannot be inherited from parents)"
        # Arrange
        parentsGT = [ [1, 1], [0, 1] ]
        b1Gt = [[0, 0],[0, 0],[0, 0],[0, 0]]
        b2Gt = [[1, 1],[0, 1],[1, 1],[0, 1]]
        truth = 0
        
        # Act
        b1Alleles, b2Alleles, ed = calculate_inheritance_ed(b1Gt, b2Gt, parentsGT)
        
        # Assert
        self.assertEqual(ed, truth, f"Expected ED to be zero but got {ed}")
        self.assertEqual(b1Alleles, 0, "Expected 0 alleles in bulk 1")
        self.assertEqual(b2Alleles, 8, "Expected 0 alleles in bulk 2")
    
    def test_calculate_inheritance_ed_with_impossible_progeny_3(self):
        "Test for impossible progeny (alleles could only be from a clone)"
        # Arrange
        parentsGT = [ [2, 1], [0, 0] ]
        b1Gt = [[2, 1]]
        b2Gt = [[0, 0]]
        truth = 0
        
        # Act
        b1Alleles, b2Alleles, ed = calculate_inheritance_ed(b1Gt, b2Gt, parentsGT)
        
        # Assert
        self.assertEqual(ed, truth, f"Expected ED to be zero but got {ed}")
        self.assertEqual(b1Alleles, 0, "Expected 0 alleles in bulk 1")
        self.assertEqual(b2Alleles, 0, "Expected 0 alleles in bulk 2")
    
    def test_calculate_ed_for_results_stability_1(self):
        "Tests done with the program in a state where the results are trustworthy; future changes should not affect the results"
        # Arrange
        parentsGT = [ [0, 1], [1, 1] ]
        b1Gt = [[0, 1], [0, 1], [0, 1], [0, 1]]
        b2Gt = [[1, 1], [0, 1], [0, 1], [1, 1]]
        truth = 0.42898458920779164
        
        # Act
        b1Alleles, b2Alleles, ed = calculate_inheritance_ed(b1Gt, b2Gt, parentsGT)
        
        # Assert
        self.assertAlmostEqual(ed, truth, places=5, msg=f"Expected ED to be {truth} but got {ed}")
    
    def test_calculate_ed_for_results_stability_2(self):
        "Tests done with the program in a state where the results are trustworthy; future changes should not affect the results"
        # Arrange
        parentsGT = [ [0, 0], [1, 1] ]
        b1Gt = [[0, 1], [0, 1], [0, 1], [0, 1]]
        b2Gt = [[1, 1], [0, 1], [0, 1], [1, 1]]
        truth = 0.25
        
        # Act
        b1Alleles, b2Alleles, ed = calculate_inheritance_ed(b1Gt, b2Gt, parentsGT)
        
        # Assert
        self.assertAlmostEqual(ed, truth, places=5, msg=f"Expected ED to be {truth} but got {ed}")
    
    def test_calculate_ed_for_results_stability_2(self):
        "Tests done with the program in a state where the results are trustworthy; future changes should not affect the results"
        # Arrange
        parentsGT = [ [0, 0], [0, 1] ]
        b1Gt = [[0, 0], [0, 1], [0, 0], [0, 1]]
        b2Gt = [[1, 1], [0, 1], [0, 1], [0, 0]]
        truth = 0.14299486306926387
        
        # Act
        b1Alleles, b2Alleles, ed = calculate_inheritance_ed(b1Gt, b2Gt, parentsGT)
        
        # Assert
        self.assertAlmostEqual(ed, truth, places=5, msg=f"Expected ED to be {truth} but got {ed}")
    
    def test_calculate_inheritance_ed_with_empty_bulk(self):
        "Test for impossible progeny (alleles could only be from a clone)"
        # Arrange
        parentsGT = [ [1, 2], [0, 0] ]
        b1Gt = []
        b2Gt = [[0, 1]]
        truth = 0
        
        # Act
        b1Alleles, b2Alleles, ed = calculate_inheritance_ed(b1Gt, b2Gt, parentsGT)
        
        # Assert
        self.assertEqual(ed, truth, f"Expected ED to be zero but got {ed}")
        self.assertEqual(b1Alleles, 0, "Expected 0 alleles in bulk 1")
        self.assertEqual(b2Alleles, 2, "Expected 0 alleles in bulk 2")
    
    def test_calculate_segregant_ed_with_empty_bulk(self):
        "Test for impossible progeny (alleles could only be from a clone)"
        # Arrange
        parentsGT = [ [1, 2], [0, 0] ]
        b1Gt = []
        b2Gt = [[0, 1]]
        truth = 0
        
        # Act
        b1Alleles, b2Alleles, ed = calculate_segregant_ed(b1Gt, b2Gt, parentsGT)
        
        # Assert
        self.assertEqual(ed, truth, f"Expected ED to be zero but got {ed}")
        self.assertEqual(b1Alleles, 0, "Expected 0 alleles in bulk 1")
        self.assertEqual(b2Alleles, 2, "Expected 0 alleles in bulk 2")

if __name__ == '__main__':
    unittest.main()
