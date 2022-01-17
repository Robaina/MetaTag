#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Evaluation of placed sequences:
1) Count placed sequences
"""

import argparse
from phyloplacement.utils import setDefaultOutputPath
from phyloplacement.placement import TaxAssignParser


parser = argparse.ArgumentParser(
    description='Count placed sequences based on taxon level, function and quality score',
    epilog='Semidán Robaina Estévez (srobaina@ull.edu.es), 2021'
    )

optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
parser._action_groups.append(optional)

required.add_argument('--taxtable', dest='taxtable', type=str, required=True,
                      help='path to placements taxonomy table file')

required.add_argument('--taxlevel', dest='taxlevel', type=str, required=True,
                      help='specify tax level to count hits')
required.add_argument('--cluster_ids', dest='cluster_ids', type=str, required=True,
                      nargs='+',
                      help=(
                          'list of space-separated target cluster ids of the reference tree '
                          'corresponding to the selecte function to filter counts'
                          )
                    )
optional.add_argument('--outfile', dest='outfile', type=str,
                      help='path to output file')


args = parser.parse_args()
if args.outfile is None:
    args.outfile = setDefaultOutputPath(args.taxtable, tag='_counts')


taxparser = TaxAssignParser(args.taxtable)
counts = taxparser.countHits(
    cluster_ids=args.cluster_ids, taxlevel=args.taxlevel,
    normalize=True, taxopath_type='taxopath'
    )
counts.to_csv(args.outfile, header=None, sep='\t')