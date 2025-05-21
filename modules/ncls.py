import numpy as np
from ncls import NCLS

class WindowedNCLS:
    def __init__(self, windowSize=0):
        self.ncls = {}
        self.values = {}
        self.windowSize = windowSize
        
        self.numPositions = 0
        self.longestContig = 0
    
    @property
    def contigs(self):
        return list(self.ncls.keys())
    
    def add(self, chrom, positions, edValues):
        '''
        Parameters:
            chrom -- a string indicating the chromosome name
            positions -- a numpy array of integers indicating the positions of the stat values
            edValues -- a numpy array of floats indicating the stat values
        '''
        if chrom in self.ncls:
            raise ValueError(f"Chromosome '{chrom}' already exists in this WindowedNCLS object")
        
        ends = positions + 1 + self.windowSize
        self.values[chrom] = edValues
        
        self.ncls[chrom] = NCLS(positions, ends, np.arange(len(edValues)))
        self.numPositions += len(positions)
        self.longestContig = max(self.longestContig, ends[-1])
    
    def find_overlap(self, chrom, start, end):
        '''
        Parameters:
            chrom -- a string indicating the chromosome name
            start -- an integer indicating the start position of the query; if negative, will
                     be set to 0 to prevent weird NCLS bug(?)
            end -- an integer indicating the end position of the query
        Returns:
            results -- an iterator of tuples containing the start position, end position,
                       and stat value
        '''
        if not chrom in self.ncls:
            raise KeyError(f"Chromosome '{chrom}' does not exist in this WindowedNCLS object")
        return (
            (windowStart, windowEnd, self.values[chrom][valueIndex])
            for windowStart, windowEnd, valueIndex in self.ncls[chrom].find_overlap(start if start >= 0 else 0, end)
        )
    
    def find_all(self, chrom):
        '''
        Parameters:
            chrom -- a string indicating the chromosome name
        Returns:
            results -- an iterator of tuples containing the start position, end position,
                       and stat value
        '''
        if not chrom in self.ncls:
            raise KeyError(f"Chromosome '{chrom}' does not exist in this WindowedNCLS object")
        return (
            (windowStart, windowEnd, self.values[chrom][valueIndex])
            for windowStart, windowEnd, valueIndex in self.find_overlap(chrom, 0, self.longestContig+1)
        )
    
    def __contains__(self, value):
        return value in self.ncls
    
    def __repr__(self):
        return "<WindowedNCLS object;num_contigs={0};numPositions={1};windowSize={2}>".format(
            len(self.ncls),
            self.numPositions,
            self.windowSize
        )
