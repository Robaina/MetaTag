#!/usr/bin/env python
# conda activate traits

import os
import argparse
import phyloplacement.wrappers as wrappers
from phyloplacement.phylotree import placeReadsOntoTree

"""
Placement:
1) Run papara or hmmalign to align query seqs to reference alignment
2) Run epa-ng to place query onto tree
3) Run gappa to obtain tree file with placed sequences
"""

parser = argparse.ArgumentParser(description='Place query sequences onto reference tree')
parser.add_argument('--aln', dest='aln', type=str,
                    help='Path to reference fasta alignment')
parser.add_argument('--tree', dest='tree', type=str,
                    help='Path to reference tree')
parser.add_argument('--query', dest='query', type=str,
                    help='Path to query peptide sequences. Query sequences should be already preprocessed to handle illegal symbols')
parser.add_argument('--out', dest='outdir', type=str,
                    help='Path to output directory')
parser.add_argument('--aln_method', dest='aln_method', type=str,
                    default='papara', choices=['papara', 'mafft'],
                    help='Choose method to align query sequences to reference alignment')
parser.add_argument('--tree_model', dest='tree_model', type=str,
                    default=None,
                    help='Provide subsitution model employed to infer tree. Can be either a valid model name or a path to the model.gz file returned by iqtree')
args = parser.parse_args()


output_placed_tree = os.path.join(args.outdir, 'ref_database.newick')
epa_jplace = os.path.join(args.outdir, 'epa_result.jplace')

def main():
    
    placeReadsOntoTree(
        input_tree=args.tree,
        tree_model=args.tree_model,
        ref_aln=args.aln,
        query_seqs=args.query,
        aln_method=args.aln_method,
        ref_prefix='_ref',
        output_dir=args.outdir
    )

    wrappers.runGappaHeatTree(
        input_jplace=epa_jplace,
        output_dir=args.outdir,
        output_prefix='epa_result',
        additional_args=None
    )

if __name__ == '__main__':
    main()