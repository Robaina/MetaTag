#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Reference tree:
1) Run muscle or mafft to perform msa on reference database
2) Run iqtree or fasttree to infer tree
"""
from __future__ import annotations

import argparse
import logging
import os

from metatag.alignment import align_peptides
from metatag.phylotree import infer_tree
from metatag.utils import set_default_output_path

logger = logging.getLogger(__name__)


def initialize_parser(arg_list: list[str] = None) -> argparse.ArgumentParser:
    """_summary_

    Args:
        arg_list (list[str], optional): _description_. Defaults to None.

    Returns:
        argparse.ArgumentParser: _description_
    """
    parser = argparse.ArgumentParser(
        description="MSA on reference database and infer reference tree",
        usage="metatag tree [-h] [args] \n",
        epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2021",
    )

    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    parser._action_groups.append(optional)

    required.add_argument(
        "--in", dest="data", type=str, required=True, help="path to reference database"
    )
    optional.add_argument(
        "--outdir", dest="outdir", type=str, help="path to output directory"
    )
    optional.add_argument(
        "--msa_method",
        dest="msa_method",
        type=str,
        default="muscle",
        choices=["muscle", "mafft"],
        help="choose method for msa",
    )
    optional.add_argument(
        "--tree_method",
        dest="tree_method",
        type=str,
        default="iqtree",
        choices=["iqtree", "fasttree"],
        help="choose method for tree inference",
    )
    optional.add_argument(
        "--tree_model",
        dest="tree_model",
        type=str,
        default="modeltest",
        help=(
            "choose substitution model for iqtree inference. "
            'Choices=["iqtest", "modeltest", "a valid model name"]. '
            "Defaults to optimal per modeltest-ng."
        ),
    )

    if arg_list is None:
        return parser.parse_args()
    else:
        return parser.parse_args(arg_list)


def run(args: argparse.ArgumentParser) -> None:
    """_summary_"""
    if args.outdir is None:
        args.outdir = set_default_output_path(args.data, only_dirname=True)
    output_aln = os.path.join(args.outdir, "ref_database.faln")

    logger.info("Aligning reference database...")
    align_peptides(
        input_fasta=args.data,
        method=args.msa_method,
        output_file=output_aln,
        additional_args=None,
    )

    logger.info("Inferring reference tree...")
    infer_tree(
        ref_aln=output_aln,
        method=args.tree_method,
        substitution_model=args.tree_model,
        output_dir=args.outdir,
        additional_args="",
    )

    logger.info("Done!")


if __name__ == "__main__":
    args = initialize_parser()
    run(args)
