import os, gzip, codecs
import numpy as np
from math import sqrt, ceil
from ncls import NCLS
from contextlib import contextmanager

def get_codec(fileName):
    try:
        f = codecs.open(fileName, encoding='utf-8', errors='strict')
        for line in f:
            break
        return "utf-8"
    except:
        try:
            f = codecs.open(fileName, encoding='utf-16', errors='strict')
            for line in f:
                break
            return "utf-16"
        except UnicodeDecodeError:
            print(f"Can't tell what codec '{fileName}' is!!")

@contextmanager
def read_gz_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename, "r", encoding=get_codec(filename)) as f:
            yield f

class DeltaNCLS:
    def __init__(self):
        self.ncls = {}
        self.values = {}
        self.windowSize = 0
        
        self.numPositions = 0
        self.longestContig = 0
    
    @property
    def contigs(self):
        return list(self.ncls.keys())
    
    def add(self, chrom, positions, deltaValues):
        '''
        Parameters:
            chrom -- a string indicating the chromosome name
            positions -- a numpy array of integers indicating the positions of the delta values
            deltaValues -- a numpy array of floats indicating the delta SNP-index values
        '''
        if chrom in self.ncls:
            raise ValueError(f"Chromosome '{chrom}' already exists in this DeltaNCLS object")
        
        ends = positions + 1 + self.windowSize
        self.values[chrom] = deltaValues
        
        self.ncls[chrom] = NCLS(positions, ends, np.arange(len(deltaValues)))
        self.numPositions += len(positions)
        self.longestContig = max(self.longestContig, ends[-1])
    
    def find_overlap(self, chrom, start, end):
        '''
        Parameters:
            chrom -- a string indicating the chromosome name
            start -- an integer indicating the start position of the query
            end -- an integer indicating the end position of the query
        Returns:
            results -- an iterator of tuples containing the start position, end position, and delta value
        '''
        if not chrom in self.ncls:
            raise ValueError(f"Chromosome '{chrom}' does not exist in this DeltaNCLS object")
        return (
            (windowStart, windowEnd, self.values[chrom][valueIndex])
            for windowStart, windowEnd, valueIndex in self.ncls[chrom].find_overlap(start, end)
        )
    
    def find_all(self, chrom):
        '''
        Parameters:
            chrom -- a string indicating the chromosome name
        Returns:
            results -- an iterator of tuples containing the start position, end position, and delta value
        '''
        if not chrom in self.ncls:
            raise ValueError(f"Chromosome '{chrom}' does not exist in this DeltaNCLS object")
        return (
            (windowStart, windowEnd, self.values[chrom][valueIndex])
            for windowStart, windowEnd, valueIndex in self.find_overlap(chrom, 0, self.longestContig+1)
        )
    
    def __contains__(self, value):
        return value in self.ncls
    
    def __repr__(self):
        return "<DeltaNCLS object;num_contigs={0};numPositions={1};windowSize={2}>".format(
            len(self.ncls),
            self.numPositions,
            self.windowSize
        )

def parse_qtlseq_as_dict(qtlseqFile):
    '''
    Parameters:
        qtlseqFile -- a string indicating the path to QTL-seq file
    Returns:
        deltaDict -- a dictionary with structure like:
                     {
                         "chr1": [[pos1, pos2, ...], [delta1, delta2, ...]],
                         "chr2": [[pos1, pos2, ...], [delta1, delta2, ...]],
                         ...
                     }
    '''
    HEADER = ["CHROM", "POSI", "variant", "bulk1_alleles", "bulk2_alleles", "euclideanDist"]
    
    # Parse the ED file
    deltaDict = {}
    starts, ends = [], []
    with read_gz_file(qtlseqFile) as fileIn:
        firstLine = True
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            # Optionally skip header
            "Not every QTL-seq index file has a header"
            if firstLine:
                firstLine = False
                if sl[0] == "CHROM":
                    continue
            
            # Parse relevant details and validate format
            chrom, posi, variant, bulk1_depth, bulk2_depth, p99, p95, \
                bulk1_SNPindex, bulk2_SNPindex, delta_SNPindex = sl
            try:
                posi = int(posi)
            except:
                raise ValueError(f"Position '{posi}' is not an integer; offending line is '{line}'")
            try:
                deltaValue = float(delta_SNPindex)
            except:
                raise ValueError(f"Delta SNP-index '{delta_SNPindex}' is not a float; offending line is '{line}'")
            
            # Store in dictionary
            deltaDict.setdefault(chrom, [[], []])
            deltaDict[chrom][0].append(posi)
            deltaDict[chrom][1].append(abs(deltaValue))
    return deltaDict

def convert_dict_to_deltancls(deltaDict):
    '''
    Parameters:
        deltaDict -- a dictionary with structure like:
                     {
                         "chr1": [[pos1, pos2, ...], [delta1, delta2, ...]],
                         "chr2": [[pos1, pos2, ...], [delta1, delta2, ...]],
                         ...
                     }
    Returns:
        deltaNCLS -- a DeltaNCLS object containing the delta SNP-index values by chromosome
                     and position
    '''
    deltaNCLS = DeltaNCLS()
    for chrom, value in deltaDict.items():
        positions = np.array(value[0])
        deltaValues = np.array(value[1])
        deltaNCLS.add(chrom, positions, deltaValues)
    return deltaNCLS
