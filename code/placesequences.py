#!/usr/bin/env python
# conda activate traits

import os
import argparse
from phyloplacement.utils import readFromPickleFile

import phyloplacement.wrappers as wrappers
from phyloplacement.phylotree import placeReadsOntoTree, relabelTree
from phyloplacement.visualization import makeFeatureMetadataTable, plotTreeInBrowser

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
                    help=(
                        'Path to query peptide sequences. '
                        'Query sequences should be already preprocessed to handle illegal symbols')
                        )
parser.add_argument('--outdir', dest='outdir', type=str,
                    help='Path to output directory')
parser.add_argument('--aln_method', dest='aln_method', type=str,
                    default='papara', choices=['papara', 'hmmalign'],
                    help='Choose method to align query sequences to reference alignment')
parser.add_argument('--tree_model', dest='tree_model', type=str,
                    default=None,
                    help=(
                        'Provide subsitution model employed to infer tree. '
                        'Can be: 1) a valid model name or 2) a path to the log file returned by iqtree')
                        )

args = parser.parse_args()
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

    print('Relabeling final tree...')
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

    makeFeatureMetadataTable(
        label_dict=label_dict,
        output_tsv=os.path.join(args.outdir, 'empress_metadata.tsv'),
        original_labels=True
    )

    print('Drawing tree in browser...')
    plotTreeInBrowser(
        input_tree=os.path.join(args.outdir, 'epa_result_relabel.newick'),
        output_dir=None,
        feature_metadata=os.path.join(args.outdir, 'empress_metadata.tsv')
    )

    print('Finished!')

if __name__ == '__main__':
    main()