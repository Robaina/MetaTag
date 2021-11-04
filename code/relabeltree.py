#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Relabel tree from label dict as pickle file
"""

import argparse

from phyloplacement.utils import readFromPickleFile, setDefaultOutputPath
from phyloplacement.phylotree import relabelTree


parser = argparse.ArgumentParser(
    description='Relabel tree based on input label dictionaries',
    epilog='Semidán Robaina Estévez (srobaina@ull.edu.es), 2021'
    )

optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
parser._action_groups.append(optional)

required.add_argument('--tree', dest='tree', type=str, required=True,
                      help='path to tree in newick format')
required.add_argument('--labels', dest='labels', type=str, required=True,
                      help='path to label dict in pickle format. More than one comma-separated path can be input')
optional.add_argument('--outfile', dest='outfile', type=str,
                      help='path to output file directory')

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