#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to perform phylogenetic tree reconstructions and
query sequence placements onto trees
"""

import os
import re
import json
from io import StringIO

from Bio import Phylo 

from phyloplacement.utils import setDefaultOutputPath
import phyloplacement.wrappers as wrappers
from phyloplacement.alignment import alignShortReadsToReferenceMSA
from phyloplacement.database.manipulation import getFastaRecordIDs, splitReferenceFromQueryAlignments


def inferTree(ref_aln: str,
              method: str = 'iqtree',
              substitution_model: str = 'TEST',
              output_dir: str = None,
              additional_args: str = None) -> None:
    """
    Infer tree from reference msa. Best substitution model
    selected by default.
    """
    if method.lower() in 'iqtree':
        wrappers.runIqTree(
        input_algns=ref_aln,
        output_dir=output_dir,
        output_prefix='ref_database',
        keep_recovery_files=True,
        substitution_model=substitution_model,
        additional_args=additional_args
        )
    elif method.lower() in 'fasttree':   
        wrappers.runFastTree(
            input_algns=ref_aln,
            output_file=os.path.join(output_dir, 'ref_database.fasttree'),
            additional_args=additional_args
        )
    else:
        raise ValueError('Wrong method, enter iqtree or fasttree')

def relabelTree(input_newick: str,
                label_dict: dict,
                output_file: str = None) -> None: 
    """
    Relabel tree leaves 
    """
    if output_file is None:
        output_file = setDefaultOutputPath(
            input_newick, tag='_relabel'
        )
    tree = next(Phylo.parse(input_newick, 'newick'))
    leaves = tree.get_terminals()
    for leaf in leaves:
        if leaf.name in label_dict.keys():
            leaf.name = label_dict[leaf.name]
    Phylo.write(tree, output_file, 'newick')

def getIqTreeModelFromLogFile(iqtree_log: str) -> str:
    """
    Parse iqtree log file and return best fit model

    If model supplied, search model in Command: iqtree ... -m 'model' 
    If not, then -m TEST or -m MFP
    If one of those, continue to line:
    Best-fit model: 'model' chosen according to BIC
    """
    with open(iqtree_log, 'r') as log:
        text = log.read()
        subtext = re.search('(?<=Command: iqtree)(.*)(?=\\n)', text).group(1)
        model = re.search('(?<=-m )(.*)(?= -bb)', subtext).group(1)
        if model.lower() in ['mfp', 'test']:
            model = re.search('(?<=Best-fit model: )(.*)(?= chosen)', text).group(1)
    return model

def placeReadsOntoTree(input_tree: str, 
                       tree_model: str,
                       ref_aln: str,
                       query_seqs: str,
                       aln_method: str = 'papara',
                       output_dir: str = None) -> None:
    """
    Performs short read placement onto phylogenetic tree
    tree_model: str, either the model name or path to log output by iqtree
    workflow example: https://github.com/Pbdas/epa-ng/wiki/Full-Stack-Example
    """
    if output_dir is None:
        output_dir = setDefaultOutputPath(query_seqs, only_dirname=True)
    else:
        output_dir = os.path.abspath(output_dir)
    
    if os.path.isfile(tree_model):
        tree_model = getIqTreeModelFromLogFile(tree_model)
        print(f'Running EPA-ng with inferred substitution model: {tree_model}')
    
    ref_ids = getFastaRecordIDs(ref_aln)
    ref_query_msa = os.path.join(
        output_dir, setDefaultOutputPath(query_seqs, extension='.faln',
                                         only_filename=True)
        )
    aln_ref_frac = os.path.splitext(ref_query_msa)[0] + '_ref_fraction.faln'
    aln_query_frac = os.path.splitext(ref_query_msa)[0] + '_query_fraction.faln'

    alignShortReadsToReferenceMSA(
        ref_msa=ref_aln,
        query_seqs=query_seqs,
        method=aln_method,
        tree_nwk=input_tree,
        output_dir=output_dir
    )
    
    splitReferenceFromQueryAlignments(
        ref_query_msa=ref_query_msa,
        ref_ids=ref_ids,
        out_dir=output_dir
    )

    wrappers.runEPAng(
        input_tree=input_tree,
        input_aln_ref=aln_ref_frac,
        input_aln_query=aln_query_frac,
        model=tree_model,
        output_dir=output_dir,
        n_threads=None,
        additional_args=None)

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