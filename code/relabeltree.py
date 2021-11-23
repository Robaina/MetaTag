#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Relabel tree and msa from label dicts as pickle files
"""

import os
import argparse

from phyloplacement.utils import readFromPickleFile, setDefaultOutputPath
from phyloplacement.database.preprocessing import setOriginalRecordIDsInFASTA
from phyloplacement.phylotree import relabelTree


parser = argparse.ArgumentParser(
    description='Relabel tree and msa based on input label dictionaries',
    epilog='Semidán Robaina Estévez (srobaina@ull.edu.es), 2021'
    )

optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
parser._action_groups.append(optional)

required.add_argument('--tree', dest='tree', type=str, required=True,
                      help='path to tree in newick format')
required.add_argument('--labels', dest='labels', type=str, required=True,
                      nargs='+',
                      help=(
                          'path to label dict in pickle format. '
                          'More than one space-separated path can be input')
                          )
optional.add_argument('--label_prefixes', dest='labelprefixes', type=str,
                      nargs='+',
                      help=(
                          'prefix(es) to be added to sequences in each label dict,'
                          'input in same order as labels.'
                          'More than one space-separated prefix can be input')
                          )
optional.add_argument('--aln', dest='aln', type=str,
                      help='path to fasta alignment file to be relabelled')
optional.add_argument('--outdir', dest='outdir', type=str,
                      help='path to output directory')

args = parser.parse_args()
if args.outdir is None:
    args.outdir = setDefaultOutputPath(args.tree, only_dirname=True)

treeout = os.path.join(args.outdir, setDefaultOutputPath(args.tree, tag='_relabel'))
if args.aln is not None:
    alnout = os.path.join(args.outdir, setDefaultOutputPath(args.aln, tag='_relabel'))

if args.labelprefixes is None:
    label_pre = ['' for _ in args.labels]
else:
    label_pre = args.labelprefixes

label_dicts = [readFromPickleFile(label_path.strip()) for label_path in args.labels]
label_dict = {
    k: prefix + v 
    for prefix, labels in zip(label_pre, label_dicts) 
    for (k, v) in labels.items()
    }

def main():
    print('* Relabelling tree...')
    relabelTree(
        input_newick=args.tree,
        label_dict=label_dict,
        output_file=treeout
    )

    if args.aln is not None:
        setOriginalRecordIDsInFASTA(
            input_fasta=args.aln, 
            label_dict=label_dict,
            output_fasta=alnout
        )

if __name__ == '__main__':
    main()