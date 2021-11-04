#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Shrink reference tree:
1) Run shrinkTree to remove outlier branches from reference tree and msa
"""

import argparse

import phyloplacement.wrappers as wrappers
from phyloplacement.utils import setDefaultOutputPath


parser = argparse.ArgumentParser(
    description='Detect and remove outlier branches from tree and msa',
    epilog='Semidán Robaina Estévez (srobaina@ull.edu.es), 2021'
    )

optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
parser._action_groups.append(optional)

required.add_argument('--tree', dest='tree', type=str, required=True,
                      help='path to tree in newick format')
required.add_argument('--aln', dest='aln', type=str, required=True,
                      help='path to reference fasta alignment')
optional.add_argument('--outdir', dest='outdir', type=str,
                      help='path to output directory')

args = parser.parse_args()
if args.outdir is None:
    args.outdir = setDefaultOutputPath(args.tree, only_dirname=True)

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