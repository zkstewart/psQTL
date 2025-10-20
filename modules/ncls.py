import numpy as np
from ncls import NCLS

class WindowedNCLS:
    '''
    Note for how ranges are handled:
    - The start and end positions of the ranges are inclusive.
    - The window size must be a positive integer
    - Values are stored in "windows" beginning at the start position
      and extending to the end position, which is inclusive. If you
      query any position within that range, you will get the
      contents within that "window".
    '''
    def __init__(self, windowSize=None):
        self.ncls = {}
        self.values = {}
        self.windowSize = windowSize
        
        self.numPositions = 0
        self.longestContig = 0
        self.isWindowedNCLS = True # allows easy type checking
    
    @property
    def contigs(self):
        return list(self.ncls.keys())
    
    @property
    def windowSize(self):
        return self._windowSize
    
    @windowSize.setter
    def windowSize(self, value):
        if value is None:
            value = 1
        if not isinstance(value, int) or value < 1:
            raise ValueError("Window size must be a positive integer")
        self._windowSize = value
    
    def add(self, chrom, positions, edValues):
        '''
        Parameters:
            chrom -- a string indicating the chromosome name
            positions -- a numpy array of integers indicating the positions of the stat values
            edValues -- a numpy array of floats indicating the stat values
        '''
        if chrom in self.ncls:
            raise ValueError(f"Chromosome '{chrom}' already exists in this WindowedNCLS object")
        
        ends = positions + self.windowSize
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
            return iter(())
            #raise KeyError(f"Chromosome '{chrom}' does not exist in this WindowedNCLS object")
        return (
            (windowStart, windowEnd, self.values[chrom][valueIndex])
            for windowStart, windowEnd, valueIndex in self.ncls[chrom].find_overlap(start if start >= 0 else 0, end+1) # end+1 for inclusive end
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
            return iter(())
            #raise KeyError(f"Chromosome '{chrom}' does not exist in this WindowedNCLS object")
        return self.find_overlap(chrom, 0, self.longestContig+1)
    
    def __contains__(self, value):
        return value in self.ncls
    
    def __repr__(self):
        return "<WindowedNCLS object;num_contigs={0};numPositions={1};windowSize={2}>".format(
            len(self.ncls),
            self.numPositions,
            self.windowSize
        )

class RangeNCLS:
    '''
    Note for how ranges are handled:
    - Ranges begin at the start position and extend to the end position,
      which is inclusive. If you query any position encompassed by a range,
      you will receive the boolean True in response; otherwise, False.
    '''
    def __init__(self):
        self.ncls = {}
        self.values = {}
        self.ranges = {}
        
        self.numRanges = 0
        self.longestContig = 0
        self.isRangeNCLS = True # allows easy type checking
    
    @property
    def contigs(self):
        return list(self.ncls.keys())
    
    def add(self, chrom, start, end, value=None):
        '''
        Parameters:
            chrom -- a string indicating the chromosome name
            start -- an integer indicating the start of the range
            end -- an integer indicating the end of the range; paired with start
            value -- (OPTIONAL); any value to associate with the specified range
        '''
        # Index the range for NCLS range query
        self.ranges.setdefault(chrom, {"starts": [], "ends": []})
        self.ranges[chrom]["starts"].append(start)
        self.ranges[chrom]["ends"].append(end + 1) # +1 for inclusive indexing
        
        # Associate a value with this range
        if value == None:
            value = 0
        self.values.setdefault(chrom, [])
        self.values[chrom].append(value)
    
    def build(self):
        self.ncls = {}
        for chrom in self.ranges.keys():
            self.ncls[chrom] = NCLS(np.array(self.ranges[chrom]["starts"]),
                                    np.array(self.ranges[chrom]["ends"]),
                                    np.arange(len(self.values[chrom])))
            self.numRanges = len(self.ranges[chrom]["starts"])
            self.longestContig = max(self.longestContig, max(self.ranges[chrom]["ends"]))
    
    def find_overlap(self, chrom, start, end=None):
        '''
        Parameters:
            chrom -- a string indicating the chromosome name
            start -- an integer indicating the start position of the query; if negative, will
                     be set to 0 to prevent weird NCLS bug(?)
            end -- OPTIONAL; an integer indicating the end position of the query OR None
                   if a single position (start) is being queried
        Returns:
            isOverlapping -- a boolean indicating whether the 
        '''
        if not chrom in self.ncls:
            return iter(())
        
        if end == None:
            return (
                (windowStart, windowEnd, self.values[chrom][valueIndex])
                for windowStart, windowEnd, valueIndex in self.ncls[chrom].find_overlap(start if start >= 0 else 0, start+1) # start+1 for inclusive end
            )
        else:
            return (
                (windowStart, windowEnd, self.values[chrom][valueIndex])
                for windowStart, windowEnd, valueIndex in self.ncls[chrom].find_overlap(start if start >= 0 else 0, end+1) # end+1 for inclusive end
            )
    
    def is_overlapping(self, chrom, start, end=None):
        '''
        Parameters:
            chrom -- a string indicating the chromosome name
            start -- an integer indicating the start position of the query; if negative, will
                     be set to 0 to prevent weird NCLS bug(?)
            end -- OPTIONAL; an integer indicating the end position of the query OR None
                   if a single position (start) is being queried
        Returns:
            isOverlapping -- a boolean indicating whether the given position (or range) is overlapping
                             a range indexed by this NCLS object.
        '''
        if not chrom in self.ncls:
            return False
        
        if end == None:
            return True if list(self.ncls[chrom].find_overlap(start, start+1)) else False # start+1 for inclusive search
        else:
            return True if list(self.ncls[chrom].find_overlap(start, end+1)) else False # end+1 for inclusive search
    
    def __contains__(self, value):
        return value in self.ncls
    
    def __repr__(self):
        return "<RangeNCLS object;num_contigs={0};numRanges={1}>".format(
            len(self.ncls),
            self.numRanges
        )
