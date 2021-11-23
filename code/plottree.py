#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
open empress tree with placements in browser. Only if script runs locally.
"""

import os
import argparse

from phyloplacement.utils import setDefaultOutputPath
from phyloplacement.visualization import makeFeatureMetadataTable, plotTreeInBrowser


parser = argparse.ArgumentParser(
    description='Place query sequences onto reference tree',
    epilog='Semidán Robaina Estévez (srobaina@ull.edu.es), 2021'
    )

optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
parser._action_groups.append(optional)

required.add_argument('--tree', dest='tree', type=str, required=True,
                      help='path to tree file')
optional.add_argument('--outdir', dest='outdir', type=str,
                      help='path to output directory')

args = parser.parse_args()
if args.outdir is None:
    args.outdir = setDefaultOutputPath(args.tree, only_dirname=True)

def main():

    print('* Drawing tree in browser...')
    # makeFeatureMetadataTable(
    #     label_dict=label_dict,
    #     output_tsv=os.path.join(args.outdir, 'empress_metadata.tsv'),
    #     original_labels=True
    # )

    plotTreeInBrowser(
        input_tree=args.tree,
        output_dir=None) #,
        #feature_metadata=os.path.join(args.outdir, 'empress_metadata.tsv')
    #)

if __name__ == '__main__':
    main()