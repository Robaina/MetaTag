#!/usr/bin/env python
# conda activate traits

import argparse
from phyloplacement.utils import readFromPickleFile, setDefaultOutputPath
from phyloplacement.phylotree import relabelTree

"""
Shrink reference tree:
1) Run shrinkTree to remove outlier branches from reference tree and msa
"""

parser = argparse.ArgumentParser(description='Relabel tree based on input label dictionaries')
parser.add_argument('--tree', dest='tree', type=str,
                    help='Path to tree in newick format')
parser.add_argument('--labels', dest='labels', type=str,
                    help='Path to label dict in pickle format. More than one comma-separated path can be input')
parser.add_argument('--outfile', dest='outfile', type=str,
                    default=None,
                    help='Path to output file directory')

args = parser.parse_args()
if args.outfile is None:
    args.outfile = setDefaultOutputPath(args.tree, tag='_relabel')
label_dicts = [readFromPickleFile(label_path.strip()) for label_path in args.labels.split(',')]
label_dict = {k: k.split('_')[1] + '_' + v for d in label_dicts for k, v in d.items()}

    
print('Relabelling tree...')
relabelTree(
    input_newick=args.tree,
    label_dict=label_dict,
    output_file=args.outfile
)