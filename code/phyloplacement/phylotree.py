#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to perform phylogenetic tree reconstructions and
query sequence placements onto trees
"""

import os
import shutil
import re

from Bio import Phylo 

from phyloplacement.utils import setDefaultOutputPath, easyPatternMatching
import phyloplacement.wrappers as wrappers


def inferTree(ref_aln: str,
              method: str = 'iqtree',
              substitution_model: str = 'modeltest',
              output_dir: str = None,
              additional_args: str = None) -> None:
    """
    Infer tree from reference msa. Best substitution model
    selected by default.
    """
    if method.lower() in 'iqtree':
        if 'modeltest' in substitution_model.lower():
            print('Selecting best subsitution model per modeltest-ng...')
            wrappers.runModelTest(
                input_algns=ref_aln,
                n_processes=None,
                output_dir=output_dir
            )
            best_model = getTreeModelFromModeltestLog(
                modeltest_log=os.path.join(output_dir, 'modeltest_result.log'),
                criterion='BIC'
            )
            modeltest_starting_tree = os.path.join(
                output_dir, 'modeltest_result.tree'
            )
        elif 'iqtest' in substitution_model.lower():
            best_model = 'TEST'
            modeltest_starting_tree = None
        else:
            best_model = substitution_model
            modeltest_starting_tree = None

        wrappers.runIqTree(
        input_algns=ref_aln,
        output_dir=output_dir,
        output_prefix='ref_database',
        keep_recovery_files=True,
        substitution_model=best_model,
        starting_tree=modeltest_starting_tree,
        additional_args=additional_args
        )
        shutil.move(
            os.path.join(output_dir, 'ref_database.contree'),
            os.path.join(output_dir, 'ref_database.newick')
            )

    elif method.lower() in 'fasttree':
        wrappers.runFastTree(
            input_algns=ref_aln,
            output_file=os.path.join(output_dir, 'ref_database.newick'),
            additional_args=additional_args
        )
    else:
        raise ValueError('Wrong method, enter iqtree or fasttree')

def sanityCheckForiTOL(label: str) -> str:
    """
    Reformat label to comply with iTOL requirements, remove:
    1. white spaces
    2. double underscores
    3. symbols outside english letters and numbers 
    """
    legal_chars = re.compile('[^a-zA-Z0-9]')
    itol_label = legal_chars.sub('_', label).replace('__', '_')
    return itol_label

def relabelTree(input_newick: str,
                label_dict: dict,
                output_file: str = None,
                iTOL=True) -> None: 
    """
    Relabel tree leaves with labels from 
    provided dictionary. If iTOL is set, then
    labels are checked for iTOL compatibility
    """
    if output_file is None:
        output_file = setDefaultOutputPath(
            input_newick, tag='_relabel'
        )
    if iTOL:
        sanity_check = sanityCheckForiTOL
    else:
        sanity_check = lambda x: x
    tree = next(Phylo.parse(input_newick, 'newick'))
    leaves = tree.get_terminals()
    for leaf in leaves:
        if leaf.name in label_dict.keys():
            leaf.name = sanity_check(label_dict[leaf.name])
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

def getTreeModelFromModeltestLog(modeltest_log: str, criterion: str = 'BIC') -> str:
    """
    Parse modeltest-ng log file and return best fit model
    according to selected criterion: BIC, AIC or AICc
    """
    with open(modeltest_log, 'r') as log:
        text = log.read()
        model = easyPatternMatching(
            easyPatternMatching(
                text, f'Best model according to {criterion}\n', '\nlnL'),
                left_pattern='Model:'
                ).strip()
        return model