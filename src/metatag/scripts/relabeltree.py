#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Relabel tree and msa from label dicts as pickle files
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from metatag.database.preprocessing import set_original_record_ids_in_fasta
from metatag.phylotree import relabel_tree
from metatag.taxonomy import TaxonomyAssigner
from metatag.utils import DictMerger, set_default_output_path

package_dir = Path(__file__).parent.parent

logger = logging.getLogger(__name__)


def initialize_label_dict(args) -> dict:
    """
    Initialize label dictionary for tree relabelling
    """
    if args.labelprefixes is None:
        label_pre = ["" for _ in args.labels]
    else:
        label_pre = args.labelprefixes
    label_dict = DictMerger.from_pickle_paths(args.labels).merge()
    prefix_label_dict = DictMerger.from_pickle_paths(args.labels).merge(
        dict_prefixes=label_pre
    )
    return label_dict, prefix_label_dict


def assign_taxonomy_to_labels(
    prefix_label_dict: dict, label_dict: dict, taxo_file: Path = None
) -> tuple[dict]:
    """
    Assign GTDB taxonomy to tree labels
    """
    if taxo_file is None:
        taxo_file = package_dir / "data" / "merged_taxonomy.tsv"
    else:
        taxo_file = Path(taxo_file).resolve()
    taxo_dict, export_label_dict, tree_label_dict = {}, {}, {}
    taxonomy = TaxonomyAssigner(taxo_file=taxo_file)
    for k, label in label_dict.items():
        taxopath = taxonomy.assign_taxonomy_to_label(label)
        taxo_dict[k] = taxopath
        export_label_dict[prefix_label_dict[k]] = taxopath

    for k, label in prefix_label_dict.items():
        tree_label_dict[k] = f"{label}_{taxo_dict[k]}"
    return tree_label_dict, export_label_dict


def export_taxonomy_table(export_label_dict: dict, output_file: Path) -> None:
    """
    Build and export table containing assigned taxonomy
    """
    lines = ["label\ttaxopath\n"]
    with open(output_file, "w") as file:
        for label, taxopath in export_label_dict.items():
            line = label + "\t" + taxopath + "\n"
            lines.append(line)
        file.writelines(lines)


def initialize_parser() -> argparse.ArgumentParser:
    """_summary_

    Returns:
        argparse.ArgumentParser: _description_
    """
    parser = argparse.ArgumentParser(
        description="Relabel tree and msa based on input label dictionaries",
        usage="metatag relabel [-h] [args] \n",
        epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2021",
    )

    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    parser._action_groups.append(optional)

    required.add_argument(
        "--tree",
        dest="tree",
        type=Path,
        required=True,
        help="path to tree in newick format",
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
        "--label_prefixes",
        dest="labelprefixes",
        type=str,
        nargs="+",
        help=(
            "prefix(es) to be added to sequences in each label dict,"
            "input in same order as labels."
            "More than one space-separated prefix can be input"
        ),
    )
    optional.add_argument(
        "--taxonomy",
        dest="taxonomy",
        action="store_true",
        default=False,
        help=("assign GTDB taxonomy to labels containing MMP identifiers"),
    )
    optional.add_argument(
        "--aln",
        dest="aln",
        type=Path,
        help="path to fasta alignment file to be relabelled",
    )
    optional.add_argument(
        "--outdir", dest="outdir", type=Path, help="path to output directory"
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
    return parser


def run(args: argparse.ArgumentParser) -> None:
    """_summary_"""
    if args.outdir is None:
        args.outdir = set_default_output_path(args.tree, only_dirname=True)
    else:
        args.outdir = Path(args.outdir).resolve()

    treeout = args.outdir / set_default_output_path(args.tree, tag="_relabel")
    taxoout = args.outdir / set_default_output_path(
        args.tree, tag="_taxonomy", extension=".tsv"
    )
    if args.aln is not None:
        alnout = args.outdir / set_default_output_path(args.aln, tag="_relabel")
    else:
        alnout = None

    logger.info("Relabelling tree...")
    label_dict, prefix_label_dict = initialize_label_dict(args)
    if args.taxonomy:
        tree_label_dict, export_label_dict = assign_taxonomy_to_labels(
            prefix_label_dict, label_dict, taxo_file=args.taxofile
        )
        export_taxonomy_table(export_label_dict, taxoout)
    else:
        tree_label_dict = prefix_label_dict

    relabel_tree(
        input_newick=args.tree, label_dict=tree_label_dict, output_file=treeout
    )

    if args.aln is not None:
        set_original_record_ids_in_fasta(
            input_fasta=args.aln, label_dict=label_dict, output_fasta=alnout
        )


if __name__ == "__main__":
    args = initialize_parser().parse_args()
    run(args)
