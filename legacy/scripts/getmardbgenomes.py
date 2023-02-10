#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Obtain MAR db genomes given a set of query sequences
"""

import os
import argparse

import phyloplacement.utils as utils
from phyloplacement.database.parsers.mardb import (
    getMARdbGenomeByEntryCode,
    getMarDBentryCode,
)
from phyloplacement.database.manipulation import is_empty_fasta


parser = argparse.ArgumentParser(
    description="Get MarDB genomes given set of sequences",
    epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2021",
)

optional = parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
parser._action_groups.append(optional)

required.add_argument(
    "--data",
    dest="data",
    type=str,
    required=True,
    help="path mardb assembly fasta file",
)
required.add_argument(
    "--ids",
    dest="ids",
    type=str,
    required=True,
    help="path to pickle with original query sequence ids",
)
optional.add_argument("--out", dest="outdir", type=str, help="path to output directory")
optional.add_argument(
    "--proc", dest="proc", type=int, help="number of processes to use"
)

args = parser.parse_args()
if args.outdir is None:
    args.outdir = utils.set_default_output_path(args.data, only_dirname=True)
if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)


def parallel_genome(input_id, input_fasta: str, output_dir: str):
    output_fasta = os.path.join(output_dir, f"{input_id}.fa")
    getMARdbGenomeByEntryCode(
        input_id, input_fasta=input_fasta, output_fasta=output_fasta, clean_seqs=True
    )


def main():
    # Retrieve mardb genomes
    label_dict = utils.read_from_pickle_file(args.ids)
    nxr_entry_codes = {getMarDBentryCode(v) for v in label_dict.values()}

    utils.parallelizeOverInputFiles(
        parallel_genome,
        input_list=nxr_entry_codes,
        n_processes=args.proc,
        input_fasta=args.data,
        output_dir=args.outdir,
    )

    # Remove fasta files without records
    for fasta in utils.full_path_list_dir(args.outdir):
        if is_empty_fasta(fasta):
            os.remove(fasta)

    print("Finished!")


if __name__ == "__main__":
    main()
