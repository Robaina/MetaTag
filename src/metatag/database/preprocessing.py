#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to preprocess sequence databases

1. Remove illegal characters from peptide sequences
2. Remove illegal symbols from file paths
3. Relabel fasta records and make dictionary with old labels
"""

import logging
import re
from pathlib import Path

import pyfastx
from Bio import SeqIO

import metatag.wrappers as wrappers
from metatag.utils import (
    save_to_pickle_file,
    set_default_output_path,
    terminal_execute,
)

logger = logging.getLogger(__name__)


def remove_duplicates_from_fasta(
    input_fasta: Path,
    output_fasta: Path = None,
    export_duplicates: bool = False,
    duplicates_file: Path = None,
) -> None:
    """
    Removes duplicate entries (either by sequence or ID) from fasta.
    """
    if output_fasta is None:
        output_fasta = set_default_output_path(input_fasta, "_noduplicates")
    else:
        output_fasta = Path(output_fasta).resolve()

    if duplicates_file is None:
        duplicates_file = set_default_output_path(
            input_fasta, "_duplicates", extension=".txt"
        )
    wrappers.run_seqkit_nodup(
        input_fasta=input_fasta,
        output_fasta=output_fasta,
        export_duplicates=export_duplicates,
        duplicates_file=duplicates_file,
    )


def merge_fastas(input_fastas_dir: Path, output_fasta: Path = None) -> None:
    """
    Merge input fasta files into a single fasta
    """
    input_fastas_dir = Path(input_fastas_dir).resolve()
    if output_fasta is None:
        output_fasta = input_fastas_dir / "merged.fasta"
    else:
        output_fasta = Path(output_fasta).resolve()
    file = open(output_fasta, "w", encoding="UTF-8")
    cmd_str = f"printf '%s\\0' * | xargs -0 cat > {output_fasta}"
    terminal_execute(cmd_str, work_dir=input_fastas_dir, suppress_shell_output=False)
    file.close()


def assert_correct_file_path(file: Path) -> None:
    """
    Remove illegal symbols from file path
    """
    file = Path(file).resolve()
    upper_lower_digits = re.compile("[^a-zA-Z0-9]")
    fdir = file.parent
    fname, ext = file.stem, file.suffix
    clean_fname = upper_lower_digits.sub("_", fname).replace("__", "_").strip("_")
    return fdir / f"{clean_fname}{ext}"


def fasta_contains_nucleotide_sequences(fasta_file: Path) -> bool:
    """
    Check whether fasta file contains nucleotide sequences
    """
    ffile = open(fasta_file, "r")
    seq_1 = next(SeqIO.parse(ffile, "fasta"))
    seq = seq_1.seq
    ffile.close()
    return is_legit_dna_sequence(seq)


def is_legit_peptide_sequence(record_seq: str) -> bool:
    """
    Assert that peptide sequence only contains valid symbols
    """
    aas = {
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    }
    seq_symbols = {s.upper() for s in record_seq}
    return seq_symbols.issubset(aas)


def is_legit_dna_sequence(record_seq: str) -> bool:
    """
    Assert that DNA sequence only contains valid symbols
    """
    nts = {"A", "G", "T", "C"}
    seq_symbols = {s.upper() for s in record_seq}
    return seq_symbols.issubset(nts)


def assert_correct_sequence_format(
    fasta_file: Path,
    output_file: Path = None,
    is_peptide: bool = True,
) -> None:
    """
    Filter out (DNA or peptide) sequences containing illegal characters
    """
    fasta_file = Path(fasta_file).resolve()
    dirname = fasta_file.parent
    fname, ext = fasta_file.stem, fasta_file.suffix

    def remove_stop_codon_signals(record_seq: str) -> str:
        return record_seq.replace("*", "")

    if output_file is None:
        output_file = dirname / f"{fname}_modified{ext}"
    else:
        output_file = Path(output_file).resolve()
    if is_peptide:
        is_legit_sequence = is_legit_peptide_sequence
    else:
        is_legit_sequence = is_legit_dna_sequence

    with open(output_file, "w") as outfile, open(fasta_file, "r") as ffile:
        fasta = pyfastx.Fasta(ffile.name, build_index=False, full_name=True)
        for record_name, record_seq in fasta:
            if is_peptide:
                record_seq = remove_stop_codon_signals(record_seq)
            if is_legit_sequence(record_seq):
                outfile.write(f">{record_name}\n{record_seq}\n")


def set_temp_record_ids_in_fasta(
    input_fasta: Path,
    output_dir: Path = None,
    prefix: str = None,
    output_fasta: Path = None,
    output_dict: Path = None,
) -> None:
    """
    Change record ids for numbers and store then in a dictionary
    """
    input_fasta = Path(input_fasta).resolve()
    if output_dir is None:
        output_dir = input_fasta.parent
    else:
        output_dir = Path(output_dir).resolve()
    if prefix is not None:
        prefix_str = prefix
    else:
        prefix_str = ""

    fasta_file = set_default_output_path(
        input_fasta, tag="_short_ids", only_filename=True
    )
    dict_file = set_default_output_path(
        input_fasta, tag="_id_dict", extension=".pickle", only_filename=True
    )
    if output_fasta is None:
        output_fasta = output_dir / fasta_file
    if output_dict is None:
        output_dict = output_dir / dict_file

    id_dict = dict()
    with open(output_fasta, "w", encoding="UTF-8") as outfasta, open(
        input_fasta, "r", encoding="UTF-8"
    ) as ffile:
        fasta = pyfastx.Fasta(ffile.name, build_index=False, full_name=True)
        for n, (record_name, record_seq) in enumerate(fasta):
            new_id = f"{prefix_str}{n}"
            id_dict[new_id] = record_name
            outfasta.write(f">{new_id}\n{record_seq}\n")
    save_to_pickle_file(id_dict, output_dict)


def set_original_record_ids_in_fasta(
    input_fasta: Path, label_dict: dict = None, output_fasta: Path = None
):
    """
    Relabel temporary record ID by original IDs
    """
    input_fasta = Path(input_fasta).resolve()
    if output_fasta is None:
        output_fasta = set_default_output_path(input_fasta, tag="_original_ids")
    else:
        output_fasta = Path(output_fasta).resolve()

    def relabel_records():
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.name in label_dict.keys():
                name = label_dict[record.name]
                record.name = name
                record.id = name
                record.description = name
            yield record

    SeqIO.write(relabel_records(), output_fasta, "fasta")


def write_record_names_to_file(
    input_fasta: Path, filter_by_tag: str = None, output_file: Path = None
):
    """
    Write a txt file containing a list of record IDs in fasta
    @params:
    filter_by_tag: set to str containing a pattern to match
    in record labels. In this case, only matched record labels
    are returned.
    """
    input_fasta = Path(input_fasta).resolve()
    if output_file is None:
        output_file = set_default_output_path(input_fasta, extension=".txt")
    else:
        output_file = Path(output_file).resolve()
    if filter_by_tag is not None:
        pattern = filter_by_tag
    else:
        pattern = ">"
    cmd_str = f"grep '{pattern}' {input_fasta} | cut -c 2- > {output_file}"
    terminal_execute(cmd_str, suppress_shell_output=False)


def fastq2fasta(
    input_fastq: Path, output_fasta: Path = None, unzip: bool = True
) -> None:
    """
    Convert Fastq to FASTA format via sed
    """
    input_fastq = Path(input_fastq).resolve()
    if output_fasta is None:
        output_fasta = set_default_output_path(input_fastq, extension=".fasta")
    else:
        output_fasta = Path(output_fasta).resolve()
    if unzip:
        input_uncompressed = input_fastq.as_posix().strip(".gz")
        cmd_str = f"gzip -d {input_fastq} > {input_uncompressed}"
        terminal_execute(cmd_str, suppress_shell_output=False)
        input_fastq = input_uncompressed

    cmd_str = f"sed -n '1~4s/^@/>/p;2~4p' {input_fastq} > {output_fasta}"
    terminal_execute(cmd_str, suppress_shell_output=False)


def is_fasta(filename: Path):
    filename = Path(filename).resolve()
    if filename.exists():
        with open(filename, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            return any(fasta)
    else:
        return False
