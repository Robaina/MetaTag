#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Evaluation of placed sequences:
1) Run gappa examine assign to infer taxonomy of placed sequences from reference tree taxonomy
2) Assign function to placed sequences
3) Count placed sequences
"""

import argparse

from phyloplacement.utils import setDefaultOutputPath, readFromPickleFile
from phyloplacement.database.preprocessing import is_fasta, writeRecordNamesToFile
from phyloplacement.placement import assignTaxonomyToPlacements


parser = argparse.ArgumentParser(
    description='Evaluate placed sequences: assgin taxonomy, function and quantify',
    epilog='Semidán Robaina Estévez (srobaina@ull.edu.es), 2021'
    )

optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
parser._action_groups.append(optional)

required.add_argument('--jplace', dest='jplace', type=str, required=True,
                      help='path to placements jplace file')
required.add_argument('--labels', dest='labels', type=str, required=True,
                      nargs='+',
                      help=(
                          'path to label dict in pickle format. '
                          'More than one space-separated path can be input')
                          )
optional.add_argument('--outgroup', dest='outgroup', type=str,
                      help=(
                          'path to text file containing IDs of sequences to be considered '
                          'as an outgroup to root the tree. It can also be a fasta file from '
                          'which sequence names will be extracted. The outgroup will be used to '
                          'recover missing taxonomic infomation by gappa examine assign. '
                          'IMPORTANT: outgroup sequences have to be fully taxonomically resolved')
                          )
optional.add_argument('--prefix', dest='prefix', type=str,
                      default='',
                      help='prefix to be added to output files')
optional.add_argument('--outdir', dest='outdir', type=str,
                      help='path to output directory')

args = parser.parse_args()
if args.outdir is None:
    args.outdir = setDefaultOutputPath(args.jplace, only_dirname=True)


def main():

    if args.outgroup is not None:
        if is_fasta(args.outgroup):
            output_file = setDefaultOutputPath(args.outgroup,
                                               tag='_record_ids', extension='.txt')
            writeRecordNamesToFile(args.outgroup, output_file)
        else:
            output_file = args.outgroups
        args_str = f'--resolve-missing-paths --root-outgroup {output_file}'
    else:
        args_str = '--resolve-missing-paths'

    args_str = None
    label_dicts = [
        readFromPickleFile(label_path.strip()) for label_path in args.labels
    ]
    label_dict = {
        k: label 
        for labels in label_dicts 
        for (k, label) in labels.items()
        }

    assignTaxonomyToPlacements(
        jplace=args.jplace,
        id_dict=label_dict,
        output_dir=args.outdir,
        output_prefix=args.prefix,
        only_best_hit=True,
        additional_args=args_str
    )

if __name__ == '__main__':
    main()