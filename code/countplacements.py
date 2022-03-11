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
optional.add_argument('--cluster_ids', dest='cluster_ids', type=str, default=None,
                      nargs='+',
                      help=(
                          'list of space-separated target cluster ids of the reference tree '
                          'corresponding to the selected function to filter counts. If not '
                          'provided, then all clusters in the tree are considered for counting.'
                          )
                    )
optional.add_argument('--score_threshold', dest='score_threshold', type=float, default=None,
                      help='cluster score threshold value to filter placement results')
optional.add_argument('--outfile', dest='outfile', type=str, default=None,
                      help='path to output results file')
optional.add_argument('--outpdf', dest='outpdf', type=str, default=None,
                      help='path to output results figure')


args = parser.parse_args()
if args.outfile is None:
    args.outfile = setDefaultOutputPath(args.taxtable, tag='_counts')


taxparser = TaxAssignParser(args.taxtable)
counts = taxparser.countHits(
    cluster_ids=args.cluster_ids,
    score_threshold=args.score_threshold,
    taxlevel=args.taxlevel,
    normalize=True, taxopath_type='taxopath'
    )
    
column_id = 'frequency'
counts.to_csv(args.outfile, header=[column_id], index=True, sep='\t')

if args.outpdf is not None:
    fig = counts.plot.pie(figsize=(15,15), title=f"Represented {args.taxlevel}").get_figure()
    fig.savefig(args.outpdf, format='pdf')