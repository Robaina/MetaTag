"""
Tools to perform phylogenetic tree reconstructions and
query sequence placements onto trees
"""

import os
from Bio import Phylo 

from phyloplacement.utils import setDefaultOutputPath
import phyloplacement.wrappers as wrappers
from phyloplacement.alignment import alignShortReadsToReferenceMSA
from phyloplacement.database.manipulation import splitReferenceFromQueryAlignments

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
        leaf.name = label_dict[leaf.name]
    Phylo.write(tree, output_file, 'newick')

def placeReadsOntoTree(input_tree: str, 
                       tree_model: str,
                       ref_aln: str,
                       query_seqs: str,
                       aln_method: str = 'papara',
                       ref_prefix: str = '_ref',
                       output_dir: str = None) -> None:
    """
    Performs short read placement onto phylogenetic tree
    tree_model: str, either the model name or path to model file output by iqtree
    workflow example: https://github.com/Pbdas/epa-ng/wiki/Full-Stack-Example
    """
    if output_dir is None:
        output_dir = setDefaultOutputPath(query_seqs, only_dirname=True)
    ref_query_msa = os.path.join(
            output_dir, setDefaultOutputPath(query_seqs, extenion='.aln,',
                                             only_filename=True)
            ),
    aln_ref_frac = os.path.splitext(ref_query_msa)[0] + '_ref_fraction.fasta.aln'
    aln_query_frac = os.path.splitext(ref_query_msa)[0] + '_query_fraction.fasta.aln'

    alignShortReadsToReferenceMSA(
        ref_msa=ref_aln,
        query_seqs=query_seqs,
        method=aln_method,
        tree_nwk=input_tree,
        output_dir=output_dir
    )
    splitReferenceFromQueryAlignments(
        ref_query_msa=ref_query_msa,
        ref_prefix=ref_prefix,
        out_dir=output_dir
    )
    wrappers.runEPAng(
        input_tree=input_tree,
        input_aln=aln_ref_frac,
        input_query=aln_query_frac,
        model=tree_model,
        output_dir=output_dir,
        n_threads=None,
        additional_args=None)
