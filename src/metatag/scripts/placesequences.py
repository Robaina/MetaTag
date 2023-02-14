#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Placement:
1) Run papara or hmmalign to align query seqs to reference alignment
2) Run epa-ng to place query onto tree
3) Run gappa to obtain tree file with placed sequences
"""
from __future__ import annotations
import os
import argparse

import metatag.wrappers as wrappers
from metatag.utils import set_default_output_path
from metatag.placement import place_reads_onto_tree


def initialize_parser(arg_list: list[str] = None) -> argparse.ArgumentParser:
    """_summary_

    Returns:
        _type_: _description_
    """
    parser = argparse.ArgumentParser(
        description="Place query sequences onto reference tree",
        usage="metatag place [-h] [args] \n",
        epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2021",
    )

    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    parser._action_groups.append(optional)

    required.add_argument(
        "--aln",
        dest="aln",
        type=str,
        required=True,
        help="path to reference fasta alignment",
    )
    required.add_argument(
        "--tree", dest="tree", type=str, required=True, help="path to reference tree"
    )
    required.add_argument(
        "--query",
        dest="query",
        type=str,
        required=True,
        help=(
            "path to query peptide sequences. \n"
            "Query sequences should be already preprocessed to handle illegal symbols"
        ),
    )
    required.add_argument(
        "--tree_model",
        dest="tree_model",
        type=str,
        required=True,
        help=(
            "provide subsitution model employed to infer tree. "
            "Can be: 1) a valid model name or 2) a path to the log file returned by iqtree"
        ),
    )
    optional.add_argument(
        "--outdir", dest="outdir", type=str, help="path to output directory"
    )
    optional.add_argument(
        "--aln_method",
        dest="aln_method",
        type=str,
        default="papara",
        choices=["papara", "hmmalign"],
        help="choose method to align query sequences to reference alignment",
    )

    if arg_list is None:
        return parser.parse_args()
    else:
        return parser.parse_args(arg_list)


def run(args: argparse.ArgumentParser) -> None:
    """_summary_

    Raises:
        ValueError: _description_
    """
    if args.outdir is None:
        args.outdir = set_default_output_path(args.aln, only_dirname=True)
    if args.tree_model is None:
        raise ValueError("Missing tree model.")
    epa_jplace = os.path.join(args.outdir, "epa_result.jplace")

    print("* Placing reads on tree...")
    place_reads_onto_tree(
        input_tree=args.tree,
        tree_model=args.tree_model,
        ref_aln=args.aln,
        query_seqs=args.query,
        aln_method=args.aln_method,
        output_dir=args.outdir,
    )

    print("* Writing tree with placements...")
    wrappers.run_gappa_graft(
        input_jplace=epa_jplace,
        output_dir=args.outdir,
        output_prefix=None,
        additional_args="--fully-resolve",
    )

    print("Finished!")


if __name__ == "__main__":
    args = initialize_parser()
    run(args)
