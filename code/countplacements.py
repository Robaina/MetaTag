#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Evaluation of placed sequences:
1) Count placed sequences
"""

import os
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
required.add_argument('--taxlevels', dest='taxlevels', type=str, required=True,
                      nargs='+',
                      help='specify space-separated tax levels to count hits')
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
optional.add_argument('--outdir', dest='outdir', type=str, default=None,
                      help='path to output directory')
optional.add_argument('--prefix', dest='outprefix', type=str, default=None,
                      help='prefix to be added to output files')


args = parser.parse_args()
if args.outdir is None:
    args.outdir = setDefaultOutputPath(args.taxtable, only_dirname=True)
if args.outprefix is None:
    args.outprefix = 'placed_'

taxparser = TaxAssignParser(args.taxtable)
for taxlevel in args.taxlevels:
    outfile = os.path.join(args.outdir, f"{args.outprefix}{taxlevel}_counts.tsv")
    outpdf = os.path.join(args.outdir, f"{args.outprefix}{taxlevel}_counts.pdf")

    counts = taxparser.countHits(
        cluster_ids=args.cluster_ids,
        score_threshold=args.score_threshold,
        taxlevel=taxlevel,
        taxopath_type='taxopath'
        )
        
    counts.to_csv(outfile, sep='\t')
    
    # .value_counts(normalize=True)
    # fig = counts.plot.pie(figsize=(15,15), title=f"Represented {taxlevel}", rotatelabels=True).get_figure()
    fig = counts.plot.bar(figsize=(15,15), title=f"Represented {taxlevel}", rot=90).get_figure()
    fig.savefig(outpdf, format='pdf')