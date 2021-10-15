#!/usr/bin/env python
# conda activate traits

import os
import argparse
import phyloplacement.wrappers as wrappers
from phyloplacement.alignment import alignPeptides
from phyloplacement.phylotree import inferTree

"""
Reference tree:
1) Run muscle or mafft to perform msa on reference database
2) Run iqtree or fasttree to infer tree
3) Remove outlier branches in tree (optional)
"""

parser = argparse.ArgumentParser(description='MSA on reference database and infer reference tree')
parser.add_argument('--in', dest='data', type=str,
                    help='Path to reference database')
parser.add_argument('--out', dest='outdir', type=str,
                    help='Path to output directory')
parser.add_argument('--msa_method', dest='msa_method', type=str,
                    default='muscle', choices=['muscle', 'mafft'],
                    help='Choose method for msa')
parser.add_argument('--tree_method', dest='tree_method', type=str,
                    default='iqtree', choices=['iqtree', 'fasttree'],
                    help='Choose method for tree inference')
parser.add_argument('--tree_model', dest='tree_model', type=str,
                    default=None,
                    help='Choose substitution model for tree inference. Defaults to optimal.')
parser.add_argument('--shrinktree', dest='shrinktree', type=bool,
                    default=False,
                    help='Run treeshrink to remove branch outliers')
args = parser.parse_args()

output_aln = os.path.join(args.outdir, 'ref_database.faa.aln')
output_tree = os.path.join(args.outdir, 'ref_database.newick')

def main():
    
    alignPeptides(
        input_fasta=args.data,
        method=args.msa_method,
        output_file=output_aln,
        additional_args=None
    )

    inferTree(
        ref_aln=output_aln,
        method=args.tree_method,
        substitution_model=args.tree_model,
        output_dir=args.outdir,
        additional_args=''
    )

    # TODO: perhaps better to run in dedicated script
    if args.shrinktree:
        wrappers.runTreeShrink(
            input_tree='',
            input_aln='',
            output_dir='',
            additional_args=None
        )

if __name__ == '__main__':
    main()