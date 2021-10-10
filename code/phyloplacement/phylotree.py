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
from phyloplacement.alignment import convertStockholmToFastaAln


path_to_papara_exec = '/home/robaina/Software/papara'

def alignPeptides():
    """
    Perform multiple alignment on a set of peptides

    NOTES:

    script: 0MusclePep.py, 0MusclePepTrimal.py
    """
    pass

def getPhyloTree():
    """
    Make phylogenetic tree out of peptide multiple aligment data

    NOTES:

    script: fasttree.sh, iqtree.sh
    """
    pass

def placeReadsOntoTree(input_tree: str, input_aln: str,
                       input_query: str,
                       query_already_aligned: bool = False,
                       output_file: str = None) -> None:
    """
    Performs short read placement onto phylogenetic tree

    NOTES:

    Runs hmmbuild, hmmalign, epa-ng

    workflow: https://github.com/Pbdas/epa-ng/wiki/Full-Stack-Example
    """
    if output_file is None:
        output_file = setDefaultOutputPath(input_query, tag='_placement',
                                           extension='.jplace')
    if not query_already_aligned:
         
        output_hmm = '/home/robaina/Documents/TRAITS/out.hmm'
        wrappers.runHMMbuild(input_aln=input_aln,
                    output_hmm=output_hmm,
                    additional_args=None)
        
        output_aln_seqs = '/home/robaina/Documents/TRAITS/aln_query.stk'
        wrappers.runHMMalign(input_aln=input_aln,
                    input_hmm=output_hmm,
                    input_seqs=input_query,
                    output_aln_seqs=output_aln_seqs,
                    additional_args=None)
        
        input_aln_query = '/home/robaina/Documents/TRAITS'
        convertStockholmToFastaAln(input_stockholm=output_aln_seqs,
                                output_fasta=input_aln_query)
    else:
        input_aln_query = input_query
    
    output_dir = '/home/robaina/Documents/TRAITS'
    wrappers.runEPAng(input_tree=input_tree, input_aln=input_aln,
             input_query=input_aln_query, output_dir=output_dir,
             n_threads=None, additional_args=None)

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

        
    
