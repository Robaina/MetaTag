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

import os
import shutil
import argparse

from metatag.utils import (
    set_default_output_path,
    TemporaryFilePath,
    TemporaryDirectoryPath,
)
from metatag.wrappers import run_prodigal
from metatag.database.preprocessing import (
    assert_correct_sequence_format,
    remove_duplicates_from_fasta,
    merge_fastas,
    set_temp_record_ids_in_fasta,
    fasta_contains_nucleotide_sequences,
)


def initialize_parser() -> argparse.ArgumentParser:
    """_summary_

    Returns:
        argparse.ArgumentParser: _description_
    """
    parser = argparse.ArgumentParser(
        description=(
            "Database preprocessing: removal of duplicated sequences and of sequences with illegal symbols. "
            "To preferentially keep one duplicate sequence over another, place preferred sequences first."
        ),
        epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2021",
    )

    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    parser._action_groups.append(optional)

    required.add_argument(
        "--in",
        dest="data",
        type=str,
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
        "--outfile", dest="outfile", type=str, help="path to output fasta file"
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
    optional.add_argument(
        "--duplicate_method",
        dest="duplicate_method",
        type=str,
        default="seqkit",
        help="choose method to remove duplicates: seqkit or biopython",
    )

    args = parser.parse_args()
    return args


def run():
    """_summary_"""
    args = initialize_parser()
    if args.idprefix is None:
        args.idprefix = "label_"
    if args.outfile is None:
        outfasta = set_default_output_path(args.data, tag="_cleaned")
    else:
        outfasta = os.path.abspath(args.outfile)
    output_dir = os.path.abspath(os.path.dirname(args.outfile))

    if os.path.isdir(args.data):
        print("Merging input files...")
        _, file_ext = os.path.splitext(os.listdir(args.data)[0])
        data_path = os.path.abspath(os.path.join(output_dir, f"merged_data{file_ext}"))
        merge_fastas(input_fastas_dir=args.data, output_fasta=data_path)
    else:
        data_path = os.path.abspath(args.data)

    if args.dna:
        is_peptide = False
    if not args.dna:
        if fasta_contains_nucleotide_sequences(data_path):
            print("Inferred data contain nucleotide sequences")
            is_peptide = False
        else:
            is_peptide = True

    with TemporaryFilePath() as tmp_file_path:
        print("* Removing duplicates...")
        duplicates_file = set_default_output_path(
            outfasta, tag="_duplicates", extension=".txt"
        )
        remove_duplicates_from_fasta(
            input_fasta=data_path,
            output_fasta=tmp_file_path,
            method=args.duplicate_method,
            export_duplicates=args.export_dup,
            duplicates_file=duplicates_file,
        )
        print("* Asserting correct sequence format...")
        assert_correct_sequence_format(
            fasta_file=tmp_file_path, output_file=outfasta, is_peptide=is_peptide
        )

    if args.translate and not is_peptide:
        outprefix = set_default_output_path(outfasta, only_filename=True)
        print("* Translating nucleotide sequences...")
        with TemporaryDirectoryPath() as tempdir:
            run_prodigal(
                input_file=outfasta,
                output_prefix=outprefix,
                output_dir=tempdir,
                metagenome=True,
                additional_args=None,
            )
            outfaa = os.path.join(tempdir, outprefix + ".faa")
            outgbk = os.path.join(tempdir, outprefix + ".gbk")

            assert_correct_sequence_format(
                fasta_file=outfaa, output_file=outfasta, is_peptide=True
            )
            shutil.move(outgbk, output_dir)
    elif args.translate and is_peptide:
        print("Data already translated!")

    if args.relabel:
        print("* Relabelling records...")
        outfasta_short = set_default_output_path(outfasta, tag="_short_ids")
        set_temp_record_ids_in_fasta(
            input_fasta=outfasta,
            output_dir=os.path.dirname(args.outfile),
            prefix=args.idprefix,
        )
        shutil.move(outfasta_short, outfasta)

    # Remove temporary merged fasta
    if os.path.isdir(args.data):
        os.remove(data_path)

    print("Finished!")


if __name__ == "__main__":
    run()
