#! python3
# This is a succinct re-implementation of the ZS_GFF3IO Feature
# and GFF3 classes that are part of the https://github.com/zkstewart/Various_scripts
# repository. Only the attributes and methods necessary for psQTL
# to function are retained.

import re, sys, os
import pandas as pd
from collections import Counter
from ncls import NCLS

from .parsing import read_gz_file

class GFF3Feature:
    IMMUTABLE = ["ID", "ftype"] # these attributes should never change once set
    def __init__(self, ID, ftype, start=None, end=None, strand=None, contig=None, children=[], parents=set()):
        self.ID = ID
        self.ftype = ftype
        self.start = start
        self.end = end
        self.strand = strand
        self.contig = contig
        
        self._children = []
        self.children = children
        self.parents = parents if isinstance(parents, set) \
                       else set(parents) if isinstance(parents, list) \
                       else set([parents]) if isinstance(parents, str) \
                       else set()
        self.isGFF3Feature = True # flag for easier type checking
    
    @staticmethod
    def make_ftype_case_appropriate(ftype):
        if ftype.lower() == "gene":
            return "gene"
        elif ftype.lower() == "mrna":
            return "mRNA"
        elif ftype.lower() == "exon":
            return "exon"
        elif ftype.lower() == "cds":
            return "CDS"
        elif ftype.lower() == "lnc_rna":
            return "lnc_RNA"
        elif ftype.lower() == "product":
            return "Product"
        else:
            return ftype
    
    @property
    def start(self):
        if self._start != None:
            return self._start
        else:
            return min([ child.start for child in self.children ])
    
    @start.setter
    def start(self, value):
        if value != None:
            try:
                value = int(value)
            except ValueError:
                raise ValueError(f"Start value of '{value}' is not an integer and is invalid for GFF3 formatting")
            if value < 1:
                raise ValueError(f"Start value '{value}' cannot be zero or negative; GFF3 positions are 1-based")
            self._start = value
        else:
            self._start = None
    
    @property
    def end(self):
        if self._end != None:
            return self._end
        else:
            return max([ child.end for child in self.children ])
    
    @end.setter
    def end(self, value):
        if value != None:
            try:
                value = int(value)
            except ValueError:
                raise ValueError(f"End value of '{value}' is not an integer and is invalid for GFF3 formatting")
            if value < 1:
                raise ValueError(f"End value '{value}' cannot be zero or negative; GFF3 positions are 1-based")
            self._end = value
        else:
            self._end = None
    
    @property
    def strand(self):
        if self._strand != None and self._strand != ".":
            return self._strand
        else:
            childStrands = [ child.strand for child in self.children if child.strand != None and child.strand != "." ]
            if len(childStrands) == 0:
                return "+" # default strand if no children have a strand
            else:
                mostCommonStrand = Counter(childStrands).most_common(1)[0][0] # basic majority vote
                return mostCommonStrand
    
    @strand.setter
    def strand(self, value):
        ACCEPTED_STRANDS = ["+", "-", "."] # "." might represent unknown strand
        if value != None:
            if not value in ACCEPTED_STRANDS:
                raise ValueError(f"Strand value '{value}' is not recognised; should be one of {ACCEPTED_STRANDS}")
            self._strand = value
        else:
            self._strand = None
    
    @property
    def children(self):
        return self._children
    
    @children.setter
    def children(self, value):
        if isinstance(value, list):
            for child in value:
                self.add_child(child) # ensure each child is added properly
        else:
            self.add_child(value) # type is not validated, but we assume it's a GFF3Feature object
    
    def add_child(self, childFeature):
        '''
        Adds a child feature to this feature's children list.
        
        Parameters:
            childFeature -- a GFF3Feature object to add as a child of this feature.
        '''
        childFeature.parents.add(self.ID) # ensure the child knows its parent
        self.children.append(childFeature)
        self.__dict__.setdefault(childFeature.ftype, [])
        self.__dict__[childFeature.ftype].append(childFeature)
    
    def __repr__(self):
        reprPairs = []
        attrsToShow = ["ID", "ftype", "contig", "coords", "strand", "parents"]
        
        for attr in attrsToShow:
            if attr == "coords":
                reprPairs.append("coords=[{0}, {1}]".format(self.start, self.end))
            elif attr == "strand":
                reprPairs.append("strand={0}".format(self.strand if self.strand != None else "."))
            else:
                reprPairs.append("{0}={1}".format(attr, self.__dict__[attr]))
        
        return "<{0};{1}>".format(";".join(reprPairs),
                                  f"children=[{', '.join([child.ID for child in self.children])}]")

class GFF3Graph:
    PARENT_INFERENCE = {
        "CDS": "mRNA",
        "exon": "mRNA",
        "mRNA": "gene",
        "lnc_RNA": "gene",
        "Product": "gene", # Product is a special case, but we treat it as a gene parent
        "gene": None  # Gene is the top-level feature, no parent
    }
    
    def __init__(self, fileLocation):
        self.fileLocation = fileLocation
        self.ftypes = {}
        self.features = {}
        self.contigs = set()
        
        self.idRegex = re.compile(r"ID=(.+?)(;|$)")
        self.parentRegex = re.compile(r"Parent=(.+?)(;|$)")
        
        self.ncls = None
        self._nclsType = None
        self._nclsIndex = None
        
        self.parse_gff3(self.fileLocation)
        self.isGFF3Graph = True # flag for easier type checking
    
    @property
    def fileLocation(self):
        return self._fileLocation
    
    @fileLocation.setter
    def fileLocation(self, value):
        if not isinstance(value, str):
            raise ValueError("File location must be a string")
        if not os.path.isfile(value):
            raise FileNotFoundError(f"GFF3 file '{value}' is not a file")
        
        self._fileLocation = value
    
    @staticmethod
    def longest_isoform(geneFeature):
        '''
        We pick out the representative gene based on length. If length is identical,
        we'll end up picking the entry listed first in the gff3 file since our > condition
        won't be met. I doubt this will happen much or at all though.
        '''
        if not hasattr(geneFeature, "mRNA"):
            raise ValueError("Longest isoform finding can only occur on features that have mRNA children")
        
        longestMrna = [None, 0]
        for mrnaFeature in geneFeature.mRNA:
            mrnaLen = 0
            
            # Determine the features to use for length calculation
            if hasattr(mrnaFeature, "CDS"):
                features = mrnaFeature.CDS
            elif hasattr(mrnaFeature, "exon"):
                features = mrnaFeature.exon
            else:
                features = []
            
            # Sum the lengths of the CDS (or exon) features
            for subFeature in features:
                mrnaLen += (subFeature.end - subFeature.start + 1)
            
            # Update the longest mRNA if this one is longer
            if mrnaLen > longestMrna[1]:
                longestMrna = [mrnaFeature, mrnaLen]
        return longestMrna[0]
    
    def parse_gff3(self, gff3File):
        '''
        Parses a GFF3 file and populates the graph with features.
        
        Parameters:
            gff3 -- a GFF3 object to parse and populate this graph with.
        '''
        # Reset the graph
        self.fileLocation = gff3File
        self.ftypes = {}
        self.features = {}
        self.contigs = set()
        
        # Parse the GFF3 file into a graph structure
        lineCount = 0
        with read_gz_file(self.fileLocation) as fileIn:
            for line in fileIn:
                lineCount += 1
                sl = line.strip("\r\n\t;'\" ").split("\t")
                
                # Skip filler and comment lines
                if line.startswith("#") or len(sl) != 9:
                    continue
                
                # Extract information from this line
                contig, source, ftype, start, end, \
                    score, strand, frame, attributes = sl
                start = int(start)
                end = int(end)
                ftype = GFF3Feature.make_ftype_case_appropriate(ftype)
                
                # Establish or populate tracking containers
                self.ftypes.setdefault(ftype, [])
                self.contigs.add(contig)
                
                # Get the ID attribute
                featureID = [ x[0] for x in self.idRegex.findall(attributes) ]
                if len(featureID) == 1:
                    featureID = featureID[0]
                elif len(featureID) == 0:
                    featureID = f"{ftype}.{len(self.ftypes[ftype]) + 1}"
                else:
                    raise ValueError(f"GFF3 parsing failed since line #{lineCount} (\"{line}\") has multiple IDs")
                
                # Get the parent ID(s)
                parentIDs = [ x[0] for x in self.parentRegex.findall(attributes) ]
                
                # Create a feature object
                feature = GFF3Feature(ID=featureID, ftype=ftype,
                                      start=start, end=end, strand=strand,
                                      contig=contig, children=[], parents=parentIDs)
                
                # Add the feature to the graph
                self.add(feature)
    
    def add(self, feature):
        # Obtain a unique ID for the feature
        featureID = feature.ID
        ongoingCount = 1
        while featureID in self.features:
            featureID = f"{feature.ID}.{ongoingCount}"
            if not featureID in self.features:
                break
            ongoingCount += 1
        
        # Store feature within the graph
        self.ftypes.setdefault(feature.ftype, []) # necessary if first occurrence of a subfeature preceeds its parent type
        self.ftypes[feature.ftype].append(featureID)
        self.features[featureID] = feature
        
        # Update graph features with parent-child relationships
        for parentID in feature.parents:
            # Associate the feature with its existing parents
            if parentID in self.features:
                self.features[parentID].add_child(feature)
            # Create a placeholder for the parent if it doesn't exist
            else:
                if feature.ftype in GFF3Graph.PARENT_INFERENCE:
                    parentFeature = GFF3Feature(parentID, GFF3Graph.PARENT_INFERENCE[feature.ftype],
                                                contig=feature.contig,
                                                children=[feature])
                    self.add(parentFeature)
                else:
                    raise ValueError("Your GFF3 is not sorted in top-down hierarchical order which has caused an error; " +
                                     f"I encountered a {feature.ftype} with ID '{featureID}' that has a parent '{parentID}' which has " + 
                                     f"not yet appeared in your GFF3 file. I am unsure what parent type to infer for " +
                                     f"'{feature.ftype}' features, so I cannot continue parsing. Sort your GFF3 file in " +
                                     "conventional top-down hierarchical order before trying again.")
    
    def create_ncls_index(self, typeToIndex="mRNA"):
        '''
        Creates an indexed NCLS structure that can be used to find range overlaps
        for the feature types of interest.
        
        Associates the created index to the .ncls field of this object instance.
        A hidden ._nclsIndex dictionary links the ncls indices to feature objects.
        
        Parameters:
            typeToIndex -- a string (case-sensitive) indicating the entry type
                           to index OR an iterable of strings indicating multiple
                            types to index.
        '''
        if isinstance(typeToIndex, str):
            typeToIndex = [typeToIndex]
        
        for indexType in typeToIndex:
            assert indexType in self.ftypes, \
                "'{0}' not found as a feature type within the parsed GFF3 ('{1}')".format(indexType, self.fileLocation)
        
        nclsIndex = {}
        starts, ends, ids = [], [], []
        
        # Add features of the specified type to our NCLS structure
        ongoingCount = 0
        for indexType in typeToIndex:
            for featureID in self.ftypes[indexType]:
                feature = self.features[featureID]
                starts.append(feature.start)
                ends.append(feature.end + 1) # NCLS indexes 0-based like a range so +1 to make this more logically compliant with gff3 1-based system
                ids.append(ongoingCount)
                nclsIndex[ongoingCount] = feature
                ongoingCount += 1
        
        # Build the NCLS object
        starts = pd.Series(starts)
        ends = pd.Series(ends)
        ids = pd.Series(ids)
        ncls = NCLS(starts.values, ends.values, ids.values)
        
        # Associate it to this instance
        self.ncls = ncls
        self._nclsType = typeToIndex
        self._nclsIndex = nclsIndex
    
    def ncls_finder(self, start, stop, field, value):
        '''
        Queries the NCLS structure to find Features that exist within the given
        start->stop range. Specifying the field and value will narrow results
        to only those that have a Feature .field with an equal (==) value.
        
        Parameters:
            start -- an integer indicating the start position of the feature to check
                     for overlaps
            end -- an integer indicating the end positon of the feature to check for
                   overlaps; this should be 1-based in GFF3 style e.g., a first
                   position of a feature would be start=1, end=1.
            field -- a string (case-sensitive) indicating the field of the Feature
                     object that we want to check. For example, if you want to find
                     features that overlap the given start->stop range on the contig
                     "X", you'd provide "contig" as the field so this function knows
                     to check the Feature.contig field for the value of "X".
            value -- a string (case-sensitive) indicating the value of the Feature
                     field we want to find. As in the above example, if you want to
                     find the value "X" within the .contig field, you'd provide "X" as
                     the value here.
        Returns:
            features -- a list containing Features that overlap the specified range.
                        These Features are NOT deepcopied, so handle them carefully.
        '''
        assert self.ncls != None and self._nclsIndex != None, \
            "Run create_ncls_index before you call this method!"
        
        overlaps = self.ncls.find_overlap(start, stop+1) # Although our ncls is already 1-based, find_overlap acts as a range. We need to +1 to keep everything logically 1-based.
        
        features = []
        for result in overlaps: # result == [start, end, index]
            feature = self._nclsIndex[result[2]]
            if feature.__dict__[field] == value:
                features.append(feature)
        
        # Return list
        return features
    
    def __getitem__(self, key):
        return self.features[key]
    
    def __len__(self):
        return len(self.features)
    
    def __iter__(self):
        return iter(self.features.values())
    
    def __contains__(self, item):
        return item.ID in self.features
    
    def has_key(self, key):
        return key in self.features
    
    def keys(self):
        return self.features.keys()
    
    def values(self):
        return self.features.values()
    
    def items(self):
        return self.features.items()
    
    def __repr__(self):
        return "<GFF3Graph object;file='{0}';num_contigs={1};{2}>".format(
            self.fileLocation,
            len(self.contigs),
            ";".join(["num_{0}={1}".format(key, len(self.ftypes[key])) for key in self.ftypes.keys()])
        )
