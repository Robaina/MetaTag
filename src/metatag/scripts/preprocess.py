#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Preprocessing:
0) Merge fasta files if more than one
1) Remove duplicates
2) Assert correct sequence format for downstream analysis
3) Translate if nucleotide
4) Relabel if asked for
"""
from __future__ import annotations

import argparse
import logging
import shutil
from pathlib import Path

from metatag.database.preprocessing import (
    assert_correct_sequence_format,
    fasta_contains_nucleotide_sequences,
    merge_fastas,
    remove_duplicates_from_fasta,
    set_temp_record_ids_in_fasta,
)
from metatag.utils import (
    TemporaryDirectoryPath,
    TemporaryFilePath,
    set_default_output_path,
)
from metatag.wrappers import run_prodigal

logger = logging.getLogger(__name__)


def initialize_parser(arg_list: list[str] = None) -> argparse.ArgumentParser:
    """_summary_

    Args:
        arg_list (list[str], optional): _description_. Defaults to None.

    Returns:
        argparse.ArgumentParser: _description_
    """
    parser = argparse.ArgumentParser(
        description=(
            "Database preprocessing: removal of duplicated sequences and of sequences with illegal symbols. "
            "To preferentially keep one duplicate sequence over another, place preferred sequences first."
        ),
        usage="metatag preprocess [-h] [args] \n",
        epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2021",
    )

    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    parser._action_groups.append(optional)

    required.add_argument(
        "--in",
        dest="data",
        type=Path,
        required=True,
        help="path to fasta file or directory containing fasta files",
    )
    optional.add_argument(
        "--dna",
        dest="dna",
        action="store_true",
        default=False,
        help="declare if sequences are nucleotides. Defaults to peptide sequences.",
    )
    optional.add_argument(
        "--translate",
        dest="translate",
        action="store_true",
        default=False,
        help="choose whether nucleotide sequences are translated with prodigal",
    )
    optional.add_argument(
        "--export-duplicates",
        dest="export_dup",
        action="store_true",
        default=False,
        help="choose whether duplicated sequences are exported to file (same directory than outfile)",
    )
    optional.add_argument(
        "--outfile", dest="outfile", type=Path, help="path to output fasta file"
    )
    optional.add_argument(
        "--relabel",
        dest="relabel",
        default=False,
        action="store_true",
        help="relabel record IDs with numeral ids",
    )
    optional.add_argument(
        "--idprefix",
        "--relabel_prefix",
        dest="idprefix",
        type=str,
        default=None,
        help="prefix to be added to sequence IDs",
    )
    if arg_list is None:
        return parser.parse_args()
    else:
        return parser.parse_args(arg_list)


def run(args: argparse.ArgumentParser) -> None:
    """_summary_

    Args:
        args (argparse.ArgumentParser): _description_
    """
    args.data = Path(args.data).resolve()
    if args.outfile is None:
        outfasta = set_default_output_path(args.data, tag="_cleaned")
    else:
        outfasta = Path(args.outfile).resolve()
    if args.idprefix is None:
        args.idprefix = "label_"
    output_dir = outfasta.parent

    if args.data.is_dir():
        logger.info("Merging input files...")
        file_ext = next(args.data.iterdir()).suffix
        data_path = output_dir / f"merged_data{file_ext}"
        merge_fastas(input_fastas_dir=args.data, output_fasta=data_path)
    else:
        data_path = args.data

    if args.dna:
        is_peptide = False
    if not args.dna:
        if fasta_contains_nucleotide_sequences(data_path):
            logger.info("Inferred data contain nucleotide sequences")
            is_peptide = False
        else:
            is_peptide = True

    with TemporaryFilePath() as tmp_file_path:
        logger.info("Removing duplicates...")
        duplicates_file = set_default_output_path(
            outfasta, tag="_duplicates", extension=".txt"
        )
        remove_duplicates_from_fasta(
            input_fasta=data_path,
            output_fasta=tmp_file_path,
            export_duplicates=args.export_dup,
            duplicates_file=duplicates_file,
        )
        logger.info("Asserting correct sequence format...")
        assert_correct_sequence_format(
            fasta_file=tmp_file_path, output_file=outfasta, is_peptide=is_peptide
        )

    if args.translate and not is_peptide:
        outprefix = set_default_output_path(outfasta, only_filename=True)
        logger.info("Translating nucleotide sequences...")
        with TemporaryDirectoryPath() as tempdir:
            run_prodigal(
                input_file=outfasta,
                output_prefix=outprefix,
                output_dir=tempdir,
                metagenome=True,
                additional_args=None,
            )
            outfaa = tempdir / f"{outprefix}.faa"
            outgbk = tempdir / f"{outprefix}.gbk"

            assert_correct_sequence_format(
                fasta_file=outfaa, output_file=outfasta, is_peptide=True
            )
            shutil.move(outgbk, output_dir)
    elif args.translate and is_peptide:
        logger.info("Data already translated!")

    if args.relabel:
        logger.info("Relabelling records...")
        outfasta_short = set_default_output_path(outfasta, tag="_short_ids")
        set_temp_record_ids_in_fasta(
            input_fasta=outfasta,
            output_dir=output_dir,
            prefix=args.idprefix,
        )
        shutil.move(outfasta_short, outfasta)

    if args.data.is_dir():
        data_path.unlink(missing_ok=True)

    logger.info("Done!")


if __name__ == "__main__":
    args = initialize_parser()
    run(args)
