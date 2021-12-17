#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to quantify and assign labels to placed sequences
"""

import re
import json
from io import StringIO

from Bio import Phylo 

import phyloplacement.wrappers as wrappers
from phyloplacement.utils import (setDefaultOutputPath,
                                  TemporaryFilePath,
                                  readFromPickleFile)
from phyloplacement.database.parsers.mardb import MMPtaxonomyAssigner


class JplaceParser():
    """
    Methods to parse jplace files, as specified in 
    https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009
    """
    def __init__(self, path_to_jplace: str) -> None:
        with open(path_to_jplace, 'r') as JSON:
            jplace = json.load(JSON)
        self.jplace = jplace
    
    @property
    def meta(self):
        """
        Print metadata
        """
        return self.jplace['metadata']

    @property
    def fields(self): 
        """
        Print data fields
        """
        return self.jplace['fields']

    @property
    def tree(self, newick=False):
        """
        Return tree in original or newick format
        Original format contains branch labels
        in curly brackets. Newick format removes
        these labels.
        """
        if newick:
            return self.newickfyTree(self.jplace['tree'])
        else:
            return self.jplace['tree']
    
    @property
    def placements(self):
        """
        Return placement objects
        """
        return self.jplace['placements']
    
    @staticmethod
    def newickfyTree(tree_str: str) -> str:
        """
        Remove branch IDs from jplace tree string
        """
        subs_tree = re.sub("\{(\d+)\}", '', tree_str)
        return next(Phylo.parse(StringIO(subs_tree), 'newick'))
    
    def buildBranchDict(self) -> dict:
        """
        Build dictionary with edge/branch numbers as keys and 
        reference tree leaves as values
        """
        def get_id(s):
            return int(re.search("\{(\d+)\}", s).group(1))
        
        original_tree = self.jplace['tree']
        tree = self.newickfyTree(original_tree)
        leaves = tree.get_terminals()

        branches = {
            get_id(original_tree[original_tree.find(leaf.name):]): leaf.name
            for leaf in leaves
        }
        return branches

    def extractPlacementFields(self, pfielddata: list) -> dict:
        """
        Get dict with placement field values from list of values
        """
        fields = self.jplace['fields']
        return {field: pfielddata[i] for i, field in enumerate(fields)}

    def selectBestPlacement(self, placement_object: dict) -> dict:
        """
        Select placement with lowest likelihood
        """
        pdata = [
            self.extractPlacementFields(pfielddata)
            for pfielddata in placement_object['p']
            ]
        lowest_like_placement = sorted(pdata, key=lambda x: x['likelihood'])[0]
        return {'p': lowest_like_placement, 'n': placement_object['n']}

    def selectBestPlacements(self):
        """
        Select placement with lowest likelihood for 
        all placement objects in placements
        """
        best_placements = [
            self.selectBestPlacement(placement)
            for placement in self.jplace['placements']
        ]
        return best_placements

    def assignClusterCommonLabelsToQueries(self) -> dict:
        """
        Assign cluster function and lowest common taxonomy
        to query sequences placed within cluster
        """
        pass

    def labelPlacedQueries(self):
        """
        Assign reference label to best query placements
        """
        ref_dict = self.buildBranchDict()
        best_placements = self.selectBestPlacements()
        labeled_placements = {
            place_object['n'][0]: ref_dict[place_object['p']['edge_num']]
            for place_object in best_placements
            if place_object['p']['edge_num'] in ref_dict.keys()
        }
        return labeled_placements


def assignTaxonomyToPlacements(jplace: str, id_dict_pickle: str,
                               output_dir: str = None,
                               output_prefix: str = None) -> None:
    """
    Assign taxonomy to placed query sequences based on
    taxonomy assigned to tree reference sequences
    """
    if output_dir is None:
        output_dir = setDefaultOutputPath(jplace, only_dirname=True)
    if output_prefix is None:
        output_prefix = setDefaultOutputPath(jplace, only_filename=True)

    with TemporaryFilePath() as temptax:
        id_dict = readFromPickleFile(id_dict_pickle)
        taxonomy = MMPtaxonomyAssigner(
            complete='../data/taxonomy/CurrentComplete.tsv',
            partial='../data/taxonomy/CurrentPartial.tsv'
            )
        taxonomy.buildGappaTaxonomyTable(id_dict, output_file=temptax)
        wrappers.runGappaAssign(
            jplace=jplace,
            taxonomy_file=temptax,
            output_dir=output_dir,
            output_prefix=output_prefix,
            additional_args=None)