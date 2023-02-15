#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Reference database:
1) Run hmmer to extract peptides of interest
2) Reduce redundancy: cd-hit and/or repset
3) Relabel entries with temporary ids to avoid donwstream conflicts
"""
from __future__ import annotations

import argparse
import logging
import os
import shutil

from metatag.database.manipulation import (
    filter_fasta_by_hmm,
    filter_fasta_by_sequence_length,
)
from metatag.database.preprocessing import (
    merge_fastas,
    remove_duplicates_from_fasta,
    set_temp_record_ids_in_fasta,
)
from metatag.database.reduction import reduce_database_redundancy
from metatag.utils import (
    DictMerger,
    TemporaryDirectoryPath,
    TemporaryFilePath,
    full_path_list_dir,
    set_default_output_path,
)

logger = logging.getLogger(__name__)


def initialize_parser(arg_list: list[str] = None) -> argparse.ArgumentParser:
    """_summary_

    Args:
        arg_list (list[str], optional): _description_. Defaults to None.

    Returns:
        argparse.ArgumentParser: _description_
    """
    parser = argparse.ArgumentParser(
        description="Build peptide reference database",
        usage="metatag database [-h] [args] \n",
        epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2021",
    )

    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    parser._action_groups.append(optional)

    required.add_argument(
        "--hmms",
        dest="hmms",
        type=str,
        nargs="+",
        required=True,
        help="a single or a list of space-separated paths to tigrfam or pfam models",
    )
    required.add_argument(
        "--in", dest="data", type=str, required=True, help="path to peptide database"
    )
    optional.add_argument(
        "--outdir",
        dest="outdir",
        type=str,
        default=None,
        help="path to output directory",
    )
    optional.add_argument(
        "--prefix",
        dest="prefix",
        type=str,
        default="",
        help="prefix to be added to output files",
    )
    optional.add_argument(
        "--relabel_prefixes",
        dest="relabel_prefixes",
        type=str,
        nargs="+",
        default=None,
        help=(
            "List of space-separated prefixes to be added to each hmm-derived set of "
            'sequences after relabelling. Only used if "--relabel" is set. '
            'Label prefixes set to "None" are assigned "ref_" by default.'
        ),
    )
    optional.add_argument(
        "--max_sizes",
        dest="maxsizes",
        default=None,
        type=str,
        nargs="+",
        help=(
            "maximum size of representative set of sequences for each hmm model. "
            'Each (space-separated) integer corresponds to a hmm model inputed in "--hmms", '
            'thus, sorted in the same order. A value of "None" may be given to a hmm model '
            "in the list, in which case the maximum number of sequences is unlimited "
            "for that hmm."
            "Defaults to full set of sequences for all hmm modells inputed."
        ),
    )
    optional.add_argument(
        "--min_seq_length",
        dest="minseqlength",
        default=None,
        type=int,
        help=("minimum sequence length in reference database. " "Defaults to zero"),
    )
    optional.add_argument(
        "--max_seq_length",
        dest="maxseqlength",
        default=None,
        type=int,
        required=False,
        help=("maximum sequence length in reference database. " "Defaults to inf"),
    )
    optional.add_argument(
        "--relabel",
        dest="relabel",
        action="store_true",
        required=False,
        default=False,
        help=(
            "relabel record IDs with numerical ids. "
            "Unrequired to build database, but highly recommended "
            "to avoid possible conflicts downstream the pipeline."
        ),
    )
    optional.add_argument(
        "--nocdhit",
        dest="nocdhit",
        action="store_true",
        required=False,
        default=False,
        help=("do not run cd-hit on peptide database"),
    )
    optional.add_argument(
        "--remove_duplicates",
        dest="noduplicates",
        action="store_true",
        required=False,
        default=False,
        help=("remove duplicated sequences from final (merged) database"),
    )
    optional.add_argument(
        "--hmmsearch_args",
        dest="hmmsearch_args",
        type=str,
        default=None,
        required=False,
        help=(
            "a string of comma-separated additional arguments for each hmm passed to hmmsearch. "
            'e.g. inputing 3 hmms: " --cut_ga --cpu 4, --cut_nc, None". '
            "IMPORTANT: the string must be preceded by a white space. "
            "A single string may be provided, in which case the same additinal arguments will be passed for each hmm. "
            'Defaults to additional arguments string: "--cut_nc". If no additional arguments are needed, provide the value "None"'
        ),
    )

    if arg_list is None:
        return parser.parse_args()
    else:
        return parser.parse_args(arg_list)


def run(args: argparse.ArgumentParser) -> None:
    """_summary_"""
    if args.maxsizes is None:
        args.maxsizes = [None for _ in args.hmms]
    else:
        args.maxsizes = [int(arg) if arg.isdigit() else None for arg in args.maxsizes]
    if args.relabel_prefixes is None:
        args.relabel_prefixes = [None for _ in args.hmms]
    if args.outdir is None:
        args.outdir = set_default_output_path(args.data, only_dirname=True)
    args.outdir = os.path.abspath(args.outdir)
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    output_fasta = os.path.join(args.outdir, f"{args.prefix}ref_database.faa")
    output_pickle_short_ids = os.path.join(
        args.outdir, f"{args.prefix}ref_database_id_dict.pickle"
    )
    if args.hmmsearch_args is None:
        args.hmmsearch_args = ",".join(["None" for _ in args.hmms])
    hmmsearch_args_list = list(map(lambda x: x.strip(), args.hmmsearch_args.split(",")))
    hmmsearch_args_list = list(
        map(lambda x: "--cut_nc" if x == "None" else x, hmmsearch_args_list)
    )
    if len(hmmsearch_args_list) < len(args.hmms):
        hmmsearch_args_list = [hmmsearch_args_list[0] for _ in args.hmms]

    logger.info("Making peptide-specific reference database...")
    with TemporaryDirectoryPath() as tempdir1, TemporaryDirectoryPath() as tempdir2:
        for n, (hmm, maxsize, prefix, hmmsearch_args) in enumerate(
            zip(args.hmms, args.maxsizes, args.relabel_prefixes, hmmsearch_args_list)
        ):
            hmm_name = os.path.basename(hmm)
            if prefix is None:
                prefix = f"ref_{n}_"
            logger.info(
                f"Processing hmm {hmm_name} with additional arguments: {hmmsearch_args}"
            )
            hmmer_output = os.path.join(args.outdir, f"hmmer_output_{hmm_name}.txt")

            with TemporaryFilePath() as tempfasta, TemporaryFilePath() as tempfasta2, TemporaryFilePath() as tempfasta3:
                filter_fasta_by_hmm(
                    hmm_model=hmm,
                    input_fasta=args.data,
                    output_fasta=tempfasta,
                    hmmer_output=hmmer_output,
                    additional_args=hmmsearch_args,
                )

                if (args.minseqlength is not None) or (args.maxseqlength is not None):
                    logger.info("Filtering sequences by established length bounds...")
                    filter_fasta_by_sequence_length(
                        input_fasta=tempfasta,
                        min_length=args.minseqlength,
                        max_length=args.maxseqlength,
                        output_fasta=tempfasta2,
                    )
                    shutil.move(tempfasta2, tempfasta)

                reduce_database_redundancy(
                    input_fasta=tempfasta,
                    output_fasta=tempfasta3,
                    cdhit=(not args.nocdhit),
                    cdhit_args=None,
                    maxsize=maxsize,
                )

                if args.relabel:
                    logger.info("Relabelling records in reference database...")
                    output_fasta_short = os.path.join(
                        tempdir2, f"{tempfasta3}_short_ids"
                    )
                    set_temp_record_ids_in_fasta(
                        input_fasta=tempfasta3, output_dir=tempdir2, prefix=prefix
                    )
                    shutil.move(output_fasta_short, tempdir1)
                else:
                    shutil.move(tempfasta3, tempdir1)

        merge_fastas(input_fastas_dir=tempdir1, output_fasta=output_fasta)

        if args.noduplicates:
            with TemporaryFilePath() as tmp_file_path:
                logger.info("Removing duplicates...")
                remove_duplicates_from_fasta(
                    input_fasta=output_fasta,
                    output_fasta=tmp_file_path,
                    method="seqkit",
                    export_duplicates=False,
                )
                shutil.move(tmp_file_path, output_fasta)

        pickle_dict_paths = [
            file for file in full_path_list_dir(tempdir2) if file.endswith(".pickle")
        ]
        if pickle_dict_paths:
            DictMerger.from_pickle_paths(pickle_dict_paths).merge(
                save_pickle_path=output_pickle_short_ids
            )

    logger.info("Done!")


if __name__ == "__main__":
    args = initialize_parser()
    run(args)
