#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Reference tree:
1) Run muscle or mafft to perform msa on reference database
2) Run iqtree or fasttree to infer tree
"""

import os
import argparse

from phyloplacement.alignment import alignPeptides
from phyloplacement.phylotree import inferTree


parser = argparse.ArgumentParser(description='MSA on reference database and infer reference tree')
parser.add_argument('--in', dest='data', type=str,
                    help='Path to reference database')
parser.add_argument('--outdir', dest='outdir', type=str,
                    help='Path to output directory')
parser.add_argument('--msa_method', dest='msa_method', type=str,
                    default='muscle', choices=['muscle', 'mafft'],
                    help='Choose method for msa')
parser.add_argument('--tree_method', dest='tree_method', type=str,
                    default='iqtree', choices=['iqtree', 'fasttree'],
                    help='Choose method for tree inference')
parser.add_argument('--tree_model', dest='tree_model', type=str,
                    default='TEST',
                    help='Choose substitution model for tree inference. Defaults to optimal.')

args = parser.parse_args()
output_aln = os.path.join(args.outdir, 'ref_database.faln')

def main():
    
    print('Aligning reference database...')
    alignPeptides(
        input_fasta=args.data,
        method=args.msa_method,
        output_file=output_aln,
        additional_args=None
    )
    
    print('Inferring reference tree...')
    inferTree(
        ref_aln=output_aln,
        method=args.tree_method,
        substitution_model=args.tree_model,
        output_dir=args.outdir,
        additional_args=''
    )

    print('Finished!')

if __name__ == '__main__':
    main()