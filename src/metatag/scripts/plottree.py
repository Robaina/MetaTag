#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
open empress tree with placements in browser. 
Only if script runs locally.
"""
from __future__ import annotations
import os
import argparse

from metatag.utils import set_default_output_path, DictMerger
from metatag.visualization import (
    make_feature_metadata_table,
    plot_tree_in_browser,
)


def initialize_parser(arg_list: list[str] = None) -> argparse.ArgumentParser:
    """_summary_

    Returns:
        argparse.ArgumentParser: _description_
    """
    parser = argparse.ArgumentParser(
        description="Place query sequences onto reference tree",
        epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2021",
    )

    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    parser._action_groups.append(optional)

    required.add_argument(
        "--tree", dest="tree", type=str, required=True, help="path to tree file"
    )
    optional.add_argument(
        "--labels",
        dest="labels",
        type=str,
        required=False,
        nargs="+",
        help=(
            "path to label dict in pickle format. "
            "More than one space-separated path can be input"
        ),
    )
    optional.add_argument(
        "--outdir", dest="outdir", type=str, help="path to output directory"
    )

    if arg_list is None:
        return parser.parse_args()
    else:
        return parser.parse_args(arg_list)


def run(args: argparse.ArgumentParser) -> None:
    """_summary_"""
    args = initialize_parser()
    if args.outdir is None:
        args.outdir = set_default_output_path(args.tree, only_dirname=True)
    os.makedirs(args.outdir, exist_ok=True)

    print("* Drawing tree in browser...")
    if args.labels is not None:
        label_dict = DictMerger.from_pickle_paths(args.labels).merge()
        make_feature_metadata_table(
            label_dict=label_dict,
            output_tsv=os.path.join(args.outdir, "empress_metadata.tsv"),
            original_labels=False,
        )
        feature_metadata = os.path.join(args.outdir, "empress_metadata.tsv")
    else:
        feature_metadata = None

    plot_tree_in_browser(
        input_tree=args.tree,
        output_dir=os.path.join(args.outdir, "empress-plot"),
        feature_metadata=feature_metadata,
    )


if __name__ == "__main__":
    args = initialize_parser()
    run(args)
