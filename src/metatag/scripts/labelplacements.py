#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Taxonomic and function al labelling of placed sequences:
1) Run gappa examine assign to infer taxonomy of placed sequences from reference tree taxonomy
2) Assign function to placed sequences
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from metatag.database.preprocessing import (
    is_fasta,
    write_record_names_to_file,
)
from metatag.placement import (
    JplaceParser,
    add_duplicates_to_assignment_table,
    assign_labels_to_placements,
)
from metatag.utils import DictMerger, set_default_output_path, init_logger



def initialize_parser() -> argparse.ArgumentParser:
    """_summary_
    Returns:
        argparse.ArgumentParser: _description_
    """
    parser = argparse.ArgumentParser(
        description="Assgin taxonomy and function to placed query sequences",
        usage="metatag assign [-h] [args] \n",
        epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2021",
    )

    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    parser._action_groups.append(optional)

    required.add_argument(
        "--jplace",
        dest="jplace",
        type=Path,
        required=True,
        help="path to placements jplace file",
    )
    required.add_argument(
        "--labels",
        dest="labels",
        type=Path,
        required=True,
        nargs="+",
        help=(
            "path to label dict in pickle format. "
            "More than one space-separated path can be input"
        ),
    )
    optional.add_argument(
        "--query_labels",
        dest="query_labels",
        type=Path,
        default=None,
        nargs="+",
        help=(
            "path to query label dict in pickle format. "
            "More than one space-separated path can be input"
        ),
    )
    optional.add_argument(
        "--ref_clusters",
        dest="ref_clusters",
        type=Path,
        default=None,
        help=(
            "path to tsv file containing cluster assignment to each reference "
            'sequence id. Must contain one column named "id" and another '
            '(tab-separated) column named "cluster"'
        ),
    )
    optional.add_argument(
        "--ref_cluster_scores",
        dest="ref_cluster_scores",
        type=Path,
        default=None,
        help=(
            "path to tsv file containing cluster quality scores assigned to each "
            'cluster in the reference tree. Contains one column named "cluster" '
            'and another (tab-separated) column named "score"'
        ),
    )
    optional.add_argument(
        "--outgroup",
        dest="outgroup",
        type=Path,
        default=None,
        help=(
            "path to text file containing IDs of sequences to be considered "
            "as an outgroup to root the tree. It can also be a fasta file from "
            "which sequence names will be extracted. It can also be a string containing "
            "a tag to filter record labels by it. The outgroup will be used to "
            "recover missing taxonomic infomation by gappa examine assign. "
        ),
    )
    optional.add_argument(
        "--prefix",
        dest="prefix",
        type=str,
        default="placed_tax_",
        help="prefix to be added to output files",
    )
    optional.add_argument(
        "--outdir", dest="outdir", type=Path, help="path to output directory"
    )
    optional.add_argument(
        "--max_placement_distance",
        dest="max_distance",
        type=float,
        default=None,
        help=(
            "Maximum allowed pendant distance to consider a placement as valid. "
            'Change distance measure with parameter: "distance_measure" (defaults to pendant length)'
        ),
    )
    optional.add_argument(
        "--distance_measure",
        dest="distance_measure",
        type=str,
        default="pendant",
        choices=["pendant", "pendant_distal_ratio", "pendant_diameter_ratio"],
        help=(
            "Choose distance measure to remove placements with distance larger than "
            '"max_placement_distance". Choose among: '
            '1. "pendant": corresponding to pendant length of placement '
            '2. "pendant_distal_ratio": ratio between pendant and distal distances '
            '3. "pendant_diameter_ratio": ratio between pendant and tree diameter (largest pairwise distance) ratio. '
            "See https://github.com/lczech/gappa/wiki for a description of distal and pendant lengths."
        ),
    )
    optional.add_argument(
        "--min_placement_lwr",
        dest="minimum_lwr",
        type=float,
        default=None,
        help=(
            "Minimum allowed placement LWR to consider a placement as valid. Values between 0 and 1."
        ),
    )
    optional.add_argument(
        "--duplicated_query_ids",
        dest="duplicated_query_ids",
        type=Path,
        default=None,
        help="path to text file containing duplicated query ids as output by seqkit rmdup",
    )
    optional.add_argument(
        "--taxonomy_file",
        dest="taxofile",
        type=Path,
        default=None,
        help=(
            "path to tsv containing taxonomy, formated like GTDB taxopaths, for each genome ID in reference database. "
            "Defaults to None, in which case a custom GTDB taxonomy database of marine prokaryotes is used."
        ),
    )
    optional.add_argument(
        "-l",
        "--log",
        dest="logfile",
        type=Path,
        default=None,
        metavar="",
        required=False,
        help="path to log file. Log not written by default.",
    )
    return parser


def run(args: argparse.ArgumentParser) -> None:
    """_summary_

    Args:
        args (argparse.ArgumentParser): _description_
    """
    logger = init_logger(args)
    if args.outdir is None:
        args.outdir = set_default_output_path(args.jplace, only_dirname=True)
    else:
        args.outdir = Path(args.outdir).resolve()

    taxtable = args.outdir / f"{args.prefix}assignments.tsv"
    ref_labels = DictMerger.from_pickle_paths(args.labels).merge()
    if args.query_labels is not None:
        query_labels = DictMerger.from_pickle_paths(args.query_labels).merge()
    else:
        query_labels = None
    outgroup_file_generated = False

    if args.outgroup is not None:
        args.outgroup = Path(args.outgroup).resolve()
        if args.outgroup.exists():
            if is_fasta(args.outgroup):
                outgroup_file = set_default_output_path(
                    args.jplace, tag="_outgroup_ids", extension=".txt"
                )
                write_record_names_to_file(args.outgroup, output_file=outgroup_file)
                outgroup_file_generated = True
            else:
                outgroup_file = args.outgroup
        else:
            matched_labels = [
                f"{ref}\n" for ref in ref_labels.keys() if args.outgroup in ref
            ]
            if not matched_labels:
                logger.error("No matched labels for given outgroup pattern")
                sys.exit(1)
            outgroup_file = set_default_output_path(
                args.jplace, tag="_outgroup_ids", extension=".txt"
            )
            with open(outgroup_file, "w") as outfile:
                outfile.writelines(matched_labels)
            outgroup_file_generated = True

        args_str = f"--resolve-missing-paths --root-outgroup {outgroup_file}"
    else:
        args_str = ""

    if args.max_distance is not None:
        logger.info(
            f'Filtering placements by maximum distance: "{args.distance_measure}" of {args.max_distance}'
        )
        filtered_jplace = set_default_output_path(args.jplace, tag="_distance_filtered")
        jplaceparser = JplaceParser(args.jplace)
        if args.distance_measure.lower() == "pendant":
            jplaceparser.filter_placements_by_max_pendant_length(
                max_pendant_length=args.max_distance, output_file=filtered_jplace
            )
        elif "pendant_distal_ratio" in args.distance_measure.lower():
            jplaceparser.filter_placements_by_max_pendant_to_distal_length_ratio(
                max_pendant_ratio=args.max_distance, output_file=filtered_jplace
            )
        elif "pendant_diameter_ratio" in args.distance_measure.lower():
            jplaceparser.filter_placements_by_max_pendant_to_tree_diameter_ratio(
                max_pendant_ratio=args.max_distance, output_file=filtered_jplace
            )
        else:
            logger.error("Distance measure unavailable. Please choose a valid one.")
            sys.exit(1)
        args.jplace = filtered_jplace

    if args.minimum_lwr is not None:
        logger.info(f"Filtering placements by minimum LWR of: {args.minimum_lwr}")
        filtered_jplace = set_default_output_path(args.jplace, tag="_lwr_filtered")
        jplaceparser = JplaceParser(args.jplace)
        jplaceparser.filter_placements_by_minimum_lwr(
            minimum_lwr=args.minimum_lwr, output_file=filtered_jplace
        )
        args.jplace = filtered_jplace

    assign_labels_to_placements(
        jplace=args.jplace,
        ref_labels=ref_labels,
        query_labels=query_labels,
        output_dir=args.outdir,
        output_prefix=args.prefix,
        only_best_hit=False,
        ref_clusters_file=args.ref_clusters,
        ref_cluster_scores_file=args.ref_cluster_scores,
        gappa_additional_args=args_str,
        only_unique_cluster=True,
        taxo_file=args.taxofile,
    )

    if args.duplicated_query_ids is not None:
        add_duplicates_to_assignment_table(taxtable, args.duplicated_query_ids)

    if outgroup_file_generated:
        outgroup_file.unlink(missing_ok=True)
    logger.info("Done!")
    logging.shutdown()


if __name__ == "__main__":
    args = initialize_parser().parse_args()
    run(args)
