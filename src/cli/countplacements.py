#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Evaluation of placed sequences:
1) Count placed sequences
"""

import os
import argparse
from metatag.utils import set_default_output_path
from metatag.placement import TaxAssignParser


parser = argparse.ArgumentParser(
    description="Count placed sequences based on taxon level, function and quality score",
    epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2021",
)

optional = parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
parser._action_groups.append(optional)

required.add_argument(
    "--taxtable",
    dest="taxtable",
    type=str,
    required=True,
    help="path to placements taxonomy table file",
)
required.add_argument(
    "--taxlevels",
    dest="taxlevels",
    type=str,
    required=True,
    nargs="+",
    help="specify space-separated tax levels to count hits",
)
optional.add_argument(
    "--cluster_ids",
    dest="cluster_ids",
    type=str,
    default=None,
    nargs="+",
    help=(
        "list of space-separated target cluster ids of the reference tree "
        "corresponding to the selected function to filter counts. If not "
        "provided, then all clusters in the tree are considered for counting."
    ),
)
optional.add_argument(
    "--score_threshold",
    dest="score_threshold",
    type=float,
    default=None,
    help="cluster score threshold value to filter placement results",
)
optional.add_argument(
    "--outdir", dest="outdir", type=str, default=None, help="path to output directory"
)
optional.add_argument(
    "--prefix",
    dest="outprefix",
    type=str,
    default=None,
    help="prefix to be added to output files",
)
optional.add_argument(
    "--export_right_queries",
    dest="export_right_queries",
    default=False,
    action="store_true",
    help=(
        "export a tsv file containing queries with cluster assginments "
        "that passed provided filters (--score_threshold and/or --cluster_ids)"
    ),
)


args = parser.parse_args()
if args.outdir is None:
    args.outdir = set_default_output_path(args.taxtable, only_dirname=True)
if args.outprefix is None:
    args.outprefix = "placed_"
if args.export_right_queries:
    path_to_query_list = os.path.join(
        args.outdir, f"{args.outprefix}rightly_cassified_queries.tsv"
    )
else:
    path_to_query_list = None

taxparser = TaxAssignParser(args.taxtable)
for taxlevel in args.taxlevels:
    outfile = os.path.join(args.outdir, f"{args.outprefix}{taxlevel}_counts.tsv")
    outpdf = os.path.join(args.outdir, f"{args.outprefix}{taxlevel}_counts.pdf")

    taxlevel_counter = taxparser.count_hits(
        cluster_ids=args.cluster_ids,
        score_threshold=args.score_threshold,
        taxopath_type="taxopath",
        path_to_query_list=path_to_query_list,
    )

    counts, fig = taxlevel_counter.get_counts(
        taxlevel=taxlevel, output_tsv=outfile, plot_type="bar", output_pdf=outpdf
    )
