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
from pathlib import Path

from metatag.alignment import align_peptides
from metatag.phylotree import infer_tree
from metatag.utils import set_default_output_path, init_logger


def initialize_parser() -> argparse.ArgumentParser:
    """_summary_
    Returns:
        argparse.ArgumentParser: _description_
    """
    parser = argparse.ArgumentParser(
        description="MSA on reference database and infer reference tree",
        usage="metatag tree [-h] [args] \n",
        epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2021",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    parser._action_groups.append(optional)

    required.add_argument(
        "--in", dest="data", type=Path, required=True, help="path to reference database"
    )
    optional.add_argument(
        "--outdir", dest="outdir", type=Path, help="path to output directory"
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
            "choose substitution model for iqtree inference. \n"
            'Choices=["iqtest", "modeltest", "a valid model name"]. \n'
            "Defaults to optimal per modeltest-ng."
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
        args.outdir = set_default_output_path(args.data, only_dirname=True)
    else:
        args.outdir = Path(args.outdir).resolve()
    output_aln = args.outdir / "ref_database.faln"

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
    logging.shutdown()


if __name__ == "__main__":
    args = initialize_parser().parse_args()
    run(args)
