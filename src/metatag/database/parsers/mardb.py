#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to process MARdb data
"""

import os
import re
import shutil

import pandas as pd
import pyfastx

from metatag.utils import (
    create_temporary_file_path,
    set_default_output_path,
    terminal_execute,
)


class MARdbLabelParser:
    """
    Parse MARdb entry label to extract coded info
    """

    def __init__(self) -> None:
        pass

    @staticmethod
    def extract_mmp_id(label: str) -> str:
        """
        Extract mardb mmp id from reference label
        """
        db_entry = re.compile("_MMP\d+")  # {7}|_MMP\d{8}|_MMP\d{9}')
        try:
            return re.search(db_entry, label).group(0).strip("_")
        except Exception:
            return ""

    def parse(self, label: str) -> dict:
        """
        Parse MarDB sequence labels to obtain contig and locus info
        """
        parsed_dict = {
            "full": label,
            "species": "",
            "mmp_id": "",
            "contig": "",
            "gene_pos": None,
            "locus_pos": None,
            "strand": "",
        }
        try:
            entry = label.split("__")[0]
            mmp_id = self.extract_mmp_id(label)
            species = entry.strip(mmp_id).strip("_")
            meta = label.split("__")[1]
            strand = meta.split("_")[-1]
            locus_pos = tuple([int(pos) for pos in meta.split("_")[-3:-1]])
            gene_pos = int(meta.split("_")[-4])
            contig = "_".join(meta.split("_")[:-4])
            parsed_dict["species"] = species
            parsed_dict["mmp_id"] = mmp_id
            parsed_dict["contig"] = contig
            parsed_dict["gene_pos"] = gene_pos
            parsed_dict["locus_pos"] = locus_pos
            parsed_dict["strand"] = strand
        except Exception:
            pass
        return parsed_dict

    def parse_from_list(self, labels=list) -> pd.DataFrame:
        """
        Parse labels in list of labels and return DataFrame
        """
        return pd.DataFrame([self.parse(label) for label in labels])


# *** Tagged as possibly trash code ***


def get_mardb_entry_code(label: str) -> str:
    db_entry = re.compile("\[mmp_id=(.*)\] ")
    return re.search(db_entry, label).group(1)


def filter_mardb_records_by_entry_codes(
    input_fasta: str, entry_codes: set, output_fasta: str = None
) -> None:
    """
    Filter records in mardb fasta file matching provided entry codes
    """
    if output_fasta is None:
        output_fasta = set_default_output_path(input_fasta, "_fitered")

    fasta = pyfastx.Fasta(input_fasta, build_index=False, full_name=True)
    with open(output_fasta, "w") as outfile:
        for record_name, record_seq in fasta:
            entry_code = get_mardb_entry_code(record_name)
            if entry_code in entry_codes:
                outfile.write(f">{record_name}\n{record_seq}\n")


def get_mardb_genome_by_entry_code(
    entry_code: str, input_fasta: str, output_fasta: str = None, clean_seqs: bool = True
) -> None:
    """
    Get full or partial genomes with given MARdb entry code.
    If clean = True, remove characters which are not letters
    """
    if output_fasta is None:
        output_fasta = set_default_output_path(
            input_fasta, tag=f"_genome_{entry_code}", extension=".fa"
        )

    def clean_output_fasta(output_fasta: str) -> None:
        """
        Check if illegal symbols in sequences,
        then remove and tag file as cleaned
        """
        not_capital_letters = re.compile("[^A-Z]")
        fname, ext = os.path.splitext(output_fasta)
        was_cleaned = False
        cleaned_fasta = f"{fname}_cleaned{ext}"
        temp_file_path = create_temporary_file_path()
        with open(output_fasta, "r") as fasta, open(temp_file_path, "a+") as tfasta:
            for line in fasta.readlines():
                if (">" not in line) and (not_capital_letters.search(line)):
                    line = not_capital_letters.sub("", line)
                    was_cleaned = True
                tfasta.write(line)
        if was_cleaned:
            shutil.move(fasta.name, cleaned_fasta)
            shutil.move(temp_file_path, cleaned_fasta)
        else:
            os.remove(temp_file_path)

    cmd_str = f"grep -A1 {entry_code} {input_fasta} > {output_fasta}"
    terminal_execute(cmd_str, suppress_shell_output=False)
    if clean_seqs:
        clean_output_fasta(output_fasta)


def relabel_mardb(label_dict: dict) -> dict:
    """
    Convert mardb long labels into short labels
    displaying mardb id and species (if present)
    """
    db_code_pattern = re.compile("\[mmp_(.*)\]")
    species_pattern = re.compile("\[(.*?)\]")

    def edit_mardb_label(label: str) -> str:
        try:
            species = re.search(
                species_pattern, re.sub(db_code_pattern, "", label)
            ).group(1)
        except Exception:
            species = "Undetermined"
        mar_id = label.split(" ")[0]
        return f"{mar_id}_{species}"

    return {k: edit_mardb_label(v) for k, v in label_dict.items()}
