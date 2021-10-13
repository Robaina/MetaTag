"""
Tools to perform phylogenetic tree reconstructions and
query sequence placements onto trees
"""

import os
import shutil
import tempfile
from Bio import Phylo 

from phyloplacement.utils import setDefaultOutputPath
import phyloplacement.wrappers as wrappers
from phyloplacement.alignment import alignShortReadsToReferenceMSA, splitReferenceFromQueryAlignments


def placeReadsOntoTree(input_tree: str, 
                       tree_model: str,
                       ref_aln: str,
                       query_seqs: str,
                       output_dir: str = None) -> None:
    """
    Performs short read placement onto phylogenetic tree
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
        method='papara',
        tree_nwk=input_tree,
        output_dir=output_dir
    )
    splitReferenceFromQueryAlignments(
        ref_query_msa=ref_query_msa,
        ref_ids=set(),
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
