#! python3
# This is a succinct re-implementation of the ZS_GFF3IO Feature
# and GFF3 classes that are part of the https://github.com/zkstewart/Various_scripts
# repository. Only the attributes and methods necessary for psQTL
# to function are retained.

import re, sys, os
import pandas as pd
from collections import OrderedDict
from ncls import NCLS

class Feature:
    def __init__(self):
        self.children = []
        self.types = {}
        self.isFeature = True
    
    def add_attributes(self, dict):
        '''
        Parameters:
            dict -- any dictionary with key: value pairs
        '''
        for key, value in dict.items():
            self.__dict__[key] = value
    
    def add_child(self, childFeature):
        '''
        Parameters:
            childFeature -- should be a Feature object, but any object
                            will be accepted.
        '''
        self.children.append(childFeature)
        try:
            self.__dict__.setdefault(childFeature.type, [])
            self.__dict__[childFeature.type].append(childFeature)
            
            self.types.setdefault(childFeature.type, [])
            self.types[childFeature.type].append(childFeature)
        except:
            pass
    
    def retrieve_child(self, childID):
        '''
        Parameters:
            childID -- a string indicating the ID field of the child to return
        Returns:
            child -- the child object/Feature that was found, OR a None value
                     when the child does not exist.
        '''
        childList = self.retrieve_all_children()
        for child in childList:
            try:
                if child.ID == childID:
                    return child
            except:
                pass
        return None
    
    def retrieve_all_children(self):
        childList = []
        for child in self.children:
            childList.append(child)
            childList += child.retrieve_all_children()
        return childList
    
    def __getitem__(self, key):
        return self.retrieve_child(key)
    
    def __repr__(self):
        reprPairs = []
        attrsToShow = ["ID", "type", "Parent", "coords"]
        
        for attr in attrsToShow:
            try:
                reprPairs.append("{0}={1}".format(attr, self.__dict__[attr]))
            except:
                pass
        
        return "<{0}>".format(";".join(reprPairs))

class GFF3:
    def __init__(self, file_location, strict_parse=True, fix_duplicated_ids=False, slim_index=False):
        self.fileLocation = file_location
        self.features = OrderedDict()
        self.types = {}
        self.contigs = set()
        self.parentTypes = set() # tells us which feature types to expect as being parents
        
        self.ncls = None
        self._nclsType = None
        self._nclsIndex = None
        
        self.isGFF3 = True
        self.parse_gff3(strict_parse=strict_parse, fix_duplicated_ids=fix_duplicated_ids, slim_index=slim_index)
    
    @staticmethod
    def make_feature_case_appropriate(featureType):
        if featureType.lower() == "gene":
            return "gene"
        elif featureType.lower() == "mrna":
            return "mRNA"
        elif featureType.lower() == "exon":
            return "exon"
        elif featureType.lower() == "cds":
            return "CDS"
        else:
            return featureType
    
    @staticmethod
    def longest_isoform(geneFeature):
        '''
        We pick out the representative gene based on length. If length is identical,
        we'll end up picking the entry listed first in the gff3 file since our > condition
        won't be met. I doubt this will happen much or at all though.
        '''
        assert hasattr(geneFeature, "mRNA"), \
            "Longest isoform finding can only occur on features that have .mRNA children"
        
        longestMrna = [None, 0]
        for mrnaFeature in geneFeature.mRNA:
            if hasattr(mrnaFeature, "CDS"):
                featType = "CDS"
            else:
                featType = "exon"
            
            mrnaLen = 0
            for subFeature in mrnaFeature.__dict__[featType]:
                mrnaLen += (subFeature.end - subFeature.start + 1)
                
            if mrnaLen > longestMrna[1]:
                longestMrna = [mrnaFeature, mrnaLen]
        return longestMrna[0]
    
    @staticmethod
    def _get_feature_coords(feature, exonOrCDS):
        '''
        Hidden function of retrieve_sequence_from_FASTA() to retrieve sorted
        coordinates lists for the exon or CDS features associated with what
        should be an mRNA feature. If it's not it'll probably crash, hence why
        this is a private method since I know exactly how it'll be used.
        '''
        coords = [f.coords for f in feature.__dict__[exonOrCDS]]
        frames = [f.frame for f in feature.__dict__[exonOrCDS]]
        forSorting = list(zip(coords, frames))
        
        if feature.strand == '+':
            forSorting.sort(key = lambda x: (int(x[0][0]), int(x[0][1])))
        else:
            forSorting.sort(key = lambda x: (-int(x[0][0]), -int(x[0][1])))
        
        coords = [c for c, f in forSorting]
        frames = [f for c, f in forSorting]
        
        return coords, frames
    
    def parse_gff3(self, strict_parse=True, fix_duplicated_ids=False, full_warning=False, slim_index=False):
        '''
        Parameters:
            strict_parse -- a boolean indicating whether this function should
                            die if the GFF3 doesn't meet strict format standards,
                            or if failing annotation values should simply be skipped
            full_warning -- a boolean indicating whether every single warning should
                            be printed, or just the first 10.
            slim_index -- a boolean indicating whether the GFF3 indexed should be fully
                          featured (slim_index == False) or if a slim index should be created
                          with minimal features detailed (slim_index == True)
        '''
        # Set up warning handling system
        warningContainer = { # this dict acts like a JSON for data storage
            "warningCount": 0,
            "warningLimit": 10,
            "hasHandledWarnings": False
        }
        def _handle_warning_message(warningContainer, message):
            if full_warning == True or warningContainer["warningCount"] < warningContainer["warningLimit"]:
                print(message)
                warningContainer["warningCount"] += 1
            if warningContainer["hasHandledWarnings"] == False and warningContainer["warningCount"] == warningContainer["warningLimit"]:
                print("Further warning messages will be suppressed.")
                warningContainer["hasHandledWarnings"] = True
        
        # Setup for slim parsing functionality
        def _format_attributes(attributes, slim_index):
            SLIM_ATTRIBUTES = ["id", "parent"]
            
            splitAttributes = []
            for a in attributes.split("="):
                if ";" in a:
                    splitAttributes += a.rsplit(";", maxsplit=1)
                else:
                    splitAttributes.append(a)
            
            if not slim_index:
                attributesDict = {splitAttributes[i]: splitAttributes[i+1] for i in range(0, len(splitAttributes)-(len(splitAttributes)%2), 2)}
            else:
                attributesDict = {
                    splitAttributes[i]: splitAttributes[i+1]
                    for i in range(0, len(splitAttributes)-(len(splitAttributes)%2), 2)
                    if splitAttributes[i].lower() in SLIM_ATTRIBUTES
                }
            return attributesDict
        
        # Gene object loop
        lineCount = 0
        with open(self.fileLocation, 'r') as fileIn:
            for line in fileIn:
                lineCount += 1
                line = line.replace('\r', '')
                
                # Skip filler and comment lines
                if line == "\n" or line.startswith("#"):
                    continue
                
                # Extract information from this line
                try:
                    contig, source, featureType, start, end, \
                        score, strand, frame, attributes \
                        = line.rstrip('\t\n').split('\t')
                except:
                    if strict_parse == True:
                        raise ValueError(
                            f"Error: Line #{lineCount} (\"{line}\") does not meet GFF3 standards; parsing failed"
                        )
                    else:
                        _handle_warning_message(warningContainer,
                            f"Warning: Line #{lineCount} (\"{line}\") does not meet GFF3 standards;" +
                            " strict parsing is disabled, so we'll just continue and hope for the best"
                        )
                        continue
                
                # Format attributes dictionary
                attributesDict = _format_attributes(attributes, slim_index)
                
                # Ensure case conformity
                featureType = GFF3.make_feature_case_appropriate(featureType)
                
                # Fix GFF3s which did not give exons an ID
                "I'm looking at you banana genome hub. You shouldn't do this."
                if 'ID' not in attributesDict and featureType.lower() == "exon":
                    parentID = attributesDict["Parent"].split(',')
                    assert len(parentID) == 1, \
                        ("I tried to fix missing exon IDs but found a sequence with >1 parent ", +
                         "i.e., {0}".format(attributesDict["Parent"]))
                    parentID = parentID[0]
                    
                    parentFeature = self.features[parentID]
                    try:
                        numExons = len(parentFeature.exon)
                    except:
                        numExons = 0
                    attributesDict["ID"] = f"{parentID}.exon{numExons+1}"
                
                # Skip un-indexable features
                if 'ID' not in attributesDict and featureType.lower() != "cds": # see the human genome GFF3 biological_region values for why this is necessary
                    continue
                self.contigs.add(contig) # we can index this contig now that we've skipped un-indexable features
                
                # Handle parent-level features
                if 'Parent' not in attributesDict: # If no Parent field this should BE the parent
                    featureID = attributesDict["ID"]
                    
                    # End parsing if duplicate ID is found (strict_parse is True)
                    if strict_parse == True:
                        assert featureID not in self.features, \
                            (f"Error: '{featureID}' feature occurs twice indicating poorly formatted file;" + 
                            f" for debugging, line #{lineCount} == {line}; parsing will stop now")
                    # Skip parsing if duplicate ID is found (strict_parse is False)
                    else:
                        if featureID in self.features:
                            _handle_warning_message(warningContainer,
                                f"Warning: '{featureID}' feature occurs more than once indicating poorly formatted file;"+
                                " strict parsing is disabled, so we will just continue and hope for the best"
                            )
                            continue
                    
                    # Create feature and populate it with details
                    feature = Feature()
                    feature.add_attributes(attributesDict)
                    if not slim_index:
                        feature.add_attributes({
                            "contig": contig, "source": source, "type": featureType,
                            "start": int(start), "end": int(end), "coords": [int(start), int(end)],
                            "score": score, "strand": strand, "frame": frame
                        })
                    else:
                        feature.add_attributes({
                            "contig": contig, "type": featureType,
                            "start": int(start), "end": int(end),
                            "strand": strand
                        })
                    
                    # Index feature
                    self.features[featureID] = feature
                    self.types.setdefault(featureType, [])
                    self.types[featureType].append(feature)
                    self.parentTypes.add(featureType)
                
                # Handle subfeatures
                else:
                    parents = attributesDict["Parent"].split(',')
                    
                    # Loop through parents and associate feature to them
                    for parentID in parents:
                        # Flexibly obtain a feature ID for normal features and for CDS that may lack IDs
                        try:
                            featureID = attributesDict["ID"]
                        except:
                            if featureType.lower() != "cds":
                                raise AttributeError(
                                    f"Error: '{featureType}' feature lacks an ID attribute and hence cannot be indexed;" +
                                    f" for debugging, line #{lineCount} == {line}; parsing will stop now"
                                )
                            else:
                                featureID = f"{parentID}.cds"
                        
                        # End parsing if parent doesn't exist (strict_parse is True)
                        if strict_parse == True:
                            assert parentID in self.features, \
                                (f"Error: '{featureID}' feature points to a non-existing parent '{parentID}' at the time of parsing;" + 
                                f" for debugging, line #{lineCount} == {line}; parsing will stop now")
                        else:
                            if parentID not in self.features:
                                _handle_warning_message(warningContainer,
                                    f"Warning: '{featureID}' feature points to a non-existing parent '{parentID}' at the time of parsing;" + 
                                    " strict parsing is disabled, so we will just continue and hope for the best"
                                )
                                continue
                        
                        # End parsing if duplicate ID is found and we aren't fixing it
                        if featureType.lower() != "cds": # CDS features are the ONLY feature allowed to have non-unique IDs
                            if self.features[parentID].retrieve_child(featureID) != None:
                                if fix_duplicated_ids:
                                    for i in range(2, len(self.features[parentID].retrieve_all_children())+3):
                                        if self.features[parentID].retrieve_child(f"{featureID}.{i}") == None:
                                            featureID = f"{featureID}.{i}"
                                            attributesDict["ID"] = featureID
                                            break
                                else:
                                    raise ValueError(
                                        (f"Error: '{featureID}' feature is associated to a parent '{parentID}' more than once;" + 
                                        f" for debugging, line #{lineCount} == {line}; parsing will stop now")
                                    )
                        
                        # Create feature and populate it with details
                        feature = Feature()
                        feature.add_attributes(attributesDict)
                        if not slim_index:
                            feature.add_attributes({
                                "contig": contig, "source": source, "type": featureType,
                                "start": int(start), "end": int(end), "coords": [int(start), int(end)],
                                "score": score, "strand": strand, "frame": frame
                            })
                        else:
                            feature.add_attributes({
                                "contig": contig, "type": featureType,
                                "start": int(start), "end": int(end),
                                "strand": strand
                            })
                        
                        # Index feature
                        if featureType.lower() != "cds": # since CDS features aren't guaranteed to have unique IDs, there's no point indexing them
                            self.features[featureID] = feature
                        self.types.setdefault(featureType, [])
                        self.types[featureType].append(feature)
                        
                        # Double-index for special child attributes
                        self._index_products(feature)
                        
                        # Add it as a child of the parent
                        "It's important to fully-specify the child before running add_child()"
                        self.features[parentID].add_child(feature)
            
            # Generate shortcut fields
            try:
                self.gene_values = self.types["gene"]
            except:
                print("WARNING: No gene features found in GFF3 file '{0}'".format(self.fileLocation))
                self.gene_values = None
            try:
                self.mrna_values = self.types["mRNA"]
            except:
                self.mrna_values = None
            
            # Sort contigs
            self.contigs = list(self.contigs)
            try:
                self.contigs.sort(key = lambda x: list(map(int, re.findall(r'\d+', x)))) # This should let us sort things like "contig1a2" and "contig1a1" and have the latter come first
            except:
                self.contigs.sort()
    
    def sort_CDS(self):
        '''
        This method will take any parent feature that indexes CDS features, and ensures
        that the CDS children are sorted in the way we would normally expect them to be.
        For example, a +ve stranded gene should have CDS sorted so that the first entry
        in a list is the "leftmost" on the chromosome (lesser coordinate). Alternatively,
        a -ve stranded gene should have CDS sorted so that the first entry in a list is
        the "rightmost" on the chromosome (greatest coordinate). The logic behind this
        is to have CDS features sorted in 5'->3' order, rather than focusing on the
        absolute chromosomal coordinates.
        '''
        assert "CDS" in self.types, \
            "There are no CDS features in this GFF3; sorting is irrelevant"
        
        for parentID in set([x.Parent for x in self.types["CDS"]]):
            parentFeature = self[parentID]
            if parentFeature.strand == "+":
                parentFeature.CDS.sort(key = lambda x: x.start)
            elif parentFeature.strand == "-":
                parentFeature.CDS.sort(key = lambda x: -x.end)
            else:
                raise ValueError(f"'{parentID} feature has CDS children but strand '{parentFeature.strand}' is unrecognised.")
    
    def sort_exon(self):
        '''
        Akin to sort_CDS(), this function does the same but for exon-indexing features.
        The differences comes in how we approach features without a proper strand
        being noted (+ve or -ve). In these cases, we'll just skip over the feature since,
        without strand info, we can't know how it's "supposed" to be sorted.
        '''
        assert "exon" in self.types, \
            "There are no exon features in this GFF3; sorting is irrelevant"
        
        for parentID in set([x.Parent for x in self.types["exon"]]):
            parentFeature = self[parentID]
            if parentFeature.strand == "+":
                parentFeature.exon.sort(key = lambda x: x.start)
            elif parentFeature.strand == "-":
                parentFeature.exon.sort(key = lambda x: -x.end)
            else:
                continue
    
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
            assert indexType in self.types, \
                "'{0}' not found as a feature type within the parsed GFF3 ('{1}')".format(indexType, self.fileLocation)
        
        nclsIndex = {}
        starts, ends, ids = [], [], []
        
        # Add features of the specified type to our NCLS structure
        ongoingCount = 0
        for indexType in typeToIndex:
            for feature in self.types[indexType]:
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
    
    def retrieve_coords(self, feature_or_featureID, sequenceType):
        '''
        Using this GFF3 instance, retrieve the exon or CDS coordinates for the corresponding
        sequenceID. Note that this function expects the provided sequenceID to directly
        point to a feature that contains CDS and/or exon values. It used to be able to accept
        gene features and return all the subfeature values, but this proved to be unwieldy.
        
        Parameters:
            feature_or_featureID -- a string corresponding to a feature within this GFF3 object.
                                    or just a feature object in general
            sequenceType -- a string corresponding to the type of sequence to retrieve
                            i.e., in the list ["CDS", "exon"]
        Returns:
            featureCoords -- a list of lists with format like:
                             [
                                 [start_1, end_1],
                                 [start_2, end_2],
                                 ...
                             ]
            startingFrames -- a list of integers indicating what the starting frame should
                              be in any translations (if applicable)
            featureTypes -- a list of strings indicating what type of feature's details
                            have been returned e.g., "mRNA" or "lnc_RNA".
        '''
        VALID_TYPES = ["cds", "exon"]
        assert sequenceType.lower() in VALID_TYPES, \
            "'{0}' is not recognised as a valid sequenceType; should be in list {1}".format(sequenceType.lower(), VALID_TYPES)
        
        # Figure out if we're handling a feature or featureID
        if isinstance(feature_or_featureID, str):
            featureID = feature_or_featureID
            assert featureID in self.features, \
                "'{0}' is not recognised as a feature within this GFF3".format(featureID)
            feature = self.features[featureID]
        elif type(feature_or_featureID).__name__ == "Feature" \
            or type(feature_or_featureID).__name__ == "ZS_GFF3IO.Feature" \
            or (hasattr(feature_or_featureID, "isFeature") and feature_or_featureID.isFeature is True):
                feature = feature_or_featureID
                featureID = feature.ID
        
        # Validate that the feature contains the relevant fields
        if sequenceType.lower() == "cds":
            assert hasattr(feature, "CDS"), \
                "CDS feature type is requested of feature '{0}' which lacks CDS".format(featureID)
        elif sequenceType.lower() == "exon":
            assert hasattr(feature, "exon"), \
                "exon feature type is requested of feature '{0}' which lacks exon".format(featureID)
        
        # Get the coordinates required for sequenceType retrieval
        if sequenceType.lower() == "exon":
            featureCoord, featureFrame = GFF3._get_feature_coords(feature, "exon")
            featureType = feature.type
            featureID = feature.ID
        elif sequenceType.lower() == "cds":
            featureCoord, featureFrame = GFF3._get_feature_coords(feature, "CDS")
            featureType = feature.type
            featureID = feature.ID
        
        # Reverse the coord lists if we're looking at a '-' model so we start at the 3' end of the gene model
        '''
        I truly have no idea why I do this. It was in my legacy code I wrote years back,
        and I have to assume I did it for a reason. I THINK it's because some GFF3s have
        truncated sequences, and when it's truncated at the 3' end of a -ve stranded gene,
        it won't accommodate that. So when we flip it, we need to take the final frame as
        our starting frame unlike the usual, non-truncated case where it'll always start
        in the 0-frame. Sorting beforehand is hence just a way of sorting the frames, not
        the coords.
        '''
        if feature.strand == '-':
            featureCoord.reverse()
        
        # Get the starting frame for the sequence
        startingFrame = featureFrame[0]
        return featureCoord, startingFrame, featureType
    
    def _index_products(self, feature):
        '''
        Hidden helper for indexing features with the .product attribute. These are found in
        some GFF3 files, and often any protein files generated from these will have the
        product ID associated to the translated features. This can prove problematic when
        trying to find their coordinates in the GFF3 since we'll find no matches.
        
        This method will index the product so it's discoverable within the GFF3 object.
        It won't be listed as a parentType (it should always be under a gene parent), and
        
        '''
        if hasattr(feature, "Product"):
            "We want the case of this to be predictable"
            feature.product = feature.Product
            delattr(feature, "Product")
        if hasattr(feature, "product"):
            self.features[feature.product] = feature
            self.types.setdefault("product", [])
            self.types["product"].append(feature)
    
    def __getitem__(self, key):
        return self.features[key]
    
    def __setitem__(self, key, item):
        self.features[key] = item
    
    def __len__(self):
        return len(self.features)
    
    def __iter__(self):
        return iter(self.features)
    
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
        return "<GFF3 object;file='{0}';num_contigs={1};{2}>".format(
            self.fileLocation,
            len(self.contigs),
            ";".join(["num_{0}={1}".format(key, len(self.types[key])) for key in self.types.keys()])
        )
