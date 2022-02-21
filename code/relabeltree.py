#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Relabel tree and msa from label dicts as pickle files
"""

import os
import argparse
from typing import Tuple, Dict

from phyloplacement.utils import readFromPickleFile, setDefaultOutputPath
from phyloplacement.database.preprocessing import setOriginalRecordIDsInFASTA
from phyloplacement.taxonomy import TaxonomyAssigner
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
optional.add_argument('--taxonomy', dest='taxonomy', action='store_true',
                      default=False, help=(
                          'assign GTDB taxonomy to labels containing MMP identifiers')
                          )
optional.add_argument('--aln', dest='aln', type=str,
                      help='path to fasta alignment file to be relabelled')
optional.add_argument('--outdir', dest='outdir', type=str,
                      help='path to output directory')

args = parser.parse_args()
if args.outdir is None:
    args.outdir = setDefaultOutputPath(args.tree, only_dirname=True)

treeout = os.path.join(args.outdir, setDefaultOutputPath(args.tree, tag='_relabel'))
taxoout = os.path.join(args.outdir, setDefaultOutputPath(args.tree, tag='_taxonomy', extension='.tsv'))
if args.aln is not None:
    alnout = os.path.join(args.outdir, setDefaultOutputPath(args.aln, tag='_relabel'))


def initializeLabelDict(args) -> dict:
    """
    Initialize label dictionary for tree relabelling
    """
    if args.labelprefixes is None:
        label_pre = ['' for _ in args.labels]
    else:
        label_pre = args.labelprefixes

    label_dicts = [
        readFromPickleFile(label_path.strip()) for label_path in args.labels
        ]
    label_dict = {
        k: prefix + label
        for prefix, labels in zip(label_pre, label_dicts) 
        for (k, label) in labels.items()
        }
    return label_dict

def assignTaxonomyToLabels(label_dict: dict) -> Tuple[Dict]:
    """
    Assign GTDB taxonomy to tree labels
    """
    taxo_dict, export_label_dict, tree_label_dict = {}, {}, {}
    taxonomy = TaxonomyAssigner(
        taxo_file='./data/taxonomy/merged_taxonomy.tsv'
    )
    for k, label in label_dict.items():
        taxopath = taxonomy.assignTaxonomyToLabel(label)
        taxo_dict[k] = taxopath
        export_label_dict[label] = taxopath
    
    for k, label in label_dict.items():
        tree_label_dict[k] = f'{label}_{taxo_dict[k]}'

    return tree_label_dict, export_label_dict

def exportTaxonomyTable(export_label_dict: dict, outfile: str) -> None:
    """
    Build and export table containing assigned taxonomy
    """
    lines = ['label\ttaxopath\n']
    with open(outfile, 'w') as file:
        for label, taxopath in export_label_dict.items():
            line = label + '\t' + taxopath + '\n'
            lines.append(line)
        file.writelines(lines)


def main():

    print('* Relabelling tree...')

    label_dict = initializeLabelDict(args)
    if args.taxonomy:
        tree_label_dict, export_label_dict = assignTaxonomyToLabels(label_dict)
        exportTaxonomyTable(export_label_dict, taxoout)
    else:
        tree_label_dict = label_dict

    relabelTree(
        input_newick=args.tree,
        label_dict=tree_label_dict,
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