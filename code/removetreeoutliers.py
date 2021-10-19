#!/usr/bin/env python
# conda activate traits

import argparse
import phyloplacement.wrappers as wrappers

"""
Shrink reference tree:
1) Run shrinkTree to remove outlier branches from reference tree and msa
"""

parser = argparse.ArgumentParser(description='Detect and remove outlier branches from tree and msa')
parser.add_argument('--tree', dest='tree', type=str,
                    help='Path to tree in newick format')
parser.add_argument('--aln', dest='aln', type=str,
                    help='Path to reference fasta alignment')
parser.add_argument('--outdir', dest='outdir', type=str,
                    help='Path to output directory')

args = parser.parse_args()

def main():
    
    print('Removing tree branch outliers')
    wrappers.runTreeShrink(
        input_tree=args.tree,
        input_aln=args.aln,
        output_dir=args.outdir,
        additional_args=None
    )

if __name__ == '__main__':
    main()