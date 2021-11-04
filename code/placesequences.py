#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Placement:
1) Run papara or hmmalign to align query seqs to reference alignment
2) Run epa-ng to place query onto tree
3) Run gappa to obtain tree file with placed sequences
"""

import os
import argparse

from phyloplacement.utils import readFromPickleFile, setDefaultOutputPath
import phyloplacement.wrappers as wrappers
from phyloplacement.phylotree import placeReadsOntoTree, relabelTree
from phyloplacement.visualization import makeFeatureMetadataTable, plotTreeInBrowser


parser = argparse.ArgumentParser(
    description='Place query sequences onto reference tree',
    epilog='Semidán Robaina Estévez (srobaina@ull.edu.es), 2021'
    )

optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
parser._action_groups.append(optional)

required.add_argument('--aln', dest='aln', type=str, required=True,
                      help='path to reference fasta alignment')
required.add_argument('--tree', dest='tree', type=str, required=True,
                      help='path to reference tree')
required.add_argument('--query', dest='query', type=str, required=True,
                      help=(
                        'path to query peptide sequences. \n'
                        'Query sequences should be already preprocessed to handle illegal symbols')
                        )
optional.add_argument('--outdir', dest='outdir', type=str,
                      help='path to output directory')
optional.add_argument('--aln_method', dest='aln_method', type=str,
                      default='papara', choices=['papara', 'hmmalign'],
                      help='choose method to align query sequences to reference alignment')
optional.add_argument('--tree_model', dest='tree_model', type=str,
                      default=None,
                      help=(
                          'provide subsitution model employed to infer tree. '
                          'Can be: 1) a valid model name or 2) a path to the log file returned by iqtree')
                          )
optional.add_argument('--plot_placements', dest='plot_placements', action='store_true',
                      default=False,
                      help='open empress tree with placements in browser. Only if script runs locally.'
                        )

args = parser.parse_args()
if args.outdir is None:
    args.outdir = setDefaultOutputPath(args.aln, only_dirname=True)
if args.tree_model is None:
    raise ValueError('Missing tree model.')
epa_jplace = os.path.join(args.outdir, 'epa_result.jplace')

def main():
    
    print('Placing reads on tree...')
    placeReadsOntoTree(
        input_tree=args.tree,
        tree_model=args.tree_model,
        ref_aln=args.aln,
        query_seqs=args.query,
        aln_method=args.aln_method,
        ref_prefix='ref_',
        output_dir=args.outdir
    )
    
    print('Writing tree with placements...')
    wrappers.runGappaGraft(
        input_jplace=epa_jplace,
        output_dir=args.outdir,
        output_prefix=None,
        additional_args=None
    )

    print('Relabelling final tree...')
    ref_dict = readFromPickleFile(
        path_to_file=os.path.join(args.outdir, 'ref_database_id_dict.pickle')
    )
    query_dict = readFromPickleFile(
        path_to_file=os.path.join(args.outdir, 'query_cleaned_id_dict.pickle')
    )
    label_dict = {**ref_dict, **query_dict}
    relabelTree(
        input_newick=os.path.join(args.outdir, 'epa_result.newick'),
        label_dict=label_dict,
        output_file=os.path.join(args.outdir, 'epa_result_relabel.newick')
    )
    
    if args.plot_placements:
        print('Drawing tree in browser...')
        makeFeatureMetadataTable(
            label_dict=label_dict,
            output_tsv=os.path.join(args.outdir, 'empress_metadata.tsv'),
            original_labels=True
        )

        plotTreeInBrowser(
            input_tree=os.path.join(args.outdir, 'epa_result_relabel.newick'),
            output_dir=None,
            feature_metadata=os.path.join(args.outdir, 'empress_metadata.tsv')
        )

    print('Finished!')

if __name__ == '__main__':
    main()