#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Shrink reference tree:
1) Run shrinkTree to remove outlier branches from reference tree and msa
"""

import argparse

import pandas as pd

import metatag.wrappers as wrappers
from metatag.utils import set_default_output_path


parser = argparse.ArgumentParser(
    description="Detect and remove outlier branches from tree and msa",
    epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2021",
)

optional = parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
parser._action_groups.append(optional)

required.add_argument(
    "--tree", dest="tree", type=str, required=True, help="path to tree in newick format"
)
optional.add_argument(
    "--ref_database",
    dest="infasta",
    type=str,
    default=None,
    help=(
        "path to reference database in fasta format, if provided, "
        "a new fasta will be returned with outlier sequences removed"
    ),
)
optional.add_argument(
    "--aln",
    dest="aln",
    type=str,
    default=None,
    help="path to reference fasta alignment",
)
optional.add_argument(
    "--outdir", dest="outdir", type=str, help="path to output directory"
)
optional.add_argument(
    "--additional_args",
    dest="addargs",
    type=str,
    default=None,
    help="addtional arguments to treeshrink passed as a string",
)

args = parser.parse_args()
if args.outdir is None:
    args.outdir = set_default_output_path(args.tree, only_dirname=True)


def read_ref_sequences_to_remove(
    treeshrink_deleted: str, except_prefix: str = None
) -> list:
    """
    Get list of reference sequences output by treeshrink (as a tsv)
    @params:
    threshrink_deleted: path to threshrink output file containing deleted nodes
    except_prefix: prefix of nodes to be excluded from the list
    """
    df = pd.read_csv(treeshrink_deleted, sep="\t", header=None)
    outliers = [str(node_id) for node_id in df.values.flatten() if not pd.isna(node_id)]
    if except_prefix is not None:
        return [
            str(node_id)
            for node_id in outliers
            if not node_id.startswith(except_prefix)
        ]
    else:
        return outliers


def main():

    print("* Removing tree branch outliers...")
    wrappers.run_tree_shrink(
        input_tree=args.tree,
        input_aln=args.aln,
        output_dir=args.outdir,
        output_deleted_nodes=True,
        additional_args=args.addargs,
    )

    # if args.infasta is not None:
    #     print('* Removing outlier sequences from FASTA...')
    #     outfasta = os.path.join(
    #         args.outdir,
    #         set_default_output_path(args.infasta, tag='_removed_outliers', only_filename=True)
    #         )
    #     seqs_to_remove = read_ref_sequences_to_remove(
    #         treeshrink_deleted='',
    #         except_prefix=''
    #     )


if __name__ == "__main__":
    main()
