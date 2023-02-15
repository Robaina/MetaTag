#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to reduce the size of the peptide-specific
reference database

Currently based on:
1. CD-HIT
2. Repset: https://onlinelibrary.wiley.com/doi/10.1002/prot.25461
"""

import os
import shutil
import warnings
from pathlib import Path

import metatag.wrappers as wrappers
from metatag.database.manipulation import filter_fasta_by_ids
from metatag.utils import (
    TemporaryDirectoryPath,
    TemporaryFilePath,
    set_default_output_path,
    terminal_execute,
)


def get_representative_set(
    input_seqs: str, input_pi: str, max_size: int = None, outfile: str = None
) -> None:
    """
    Runs repset.py to obtain a representative
    set of size equal to max_size (or smaller if less sequences than max_size)
    or an ordered list (by 'representativeness') of representative sequences
    if max_size set to None.
    """
    input_seqs = os.path.abspath(input_seqs)
    input_pi = os.path.abspath(input_pi)
    repset_exe = Path(__file__).parent.parent / "vendor" / "repset_min.py"

    if outfile is None:
        outfile = set_default_output_path(input_seqs, tag="_repset")

    with TemporaryDirectoryPath() as tempdir:
        cmd_str = (
            f"python {repset_exe} --seqs {input_seqs} --pi {input_pi} "
            f"--outdir {tempdir} --size {max_size}"
        )
        terminal_execute(cmd_str, suppress_shell_output=True)

        with open(os.path.join(tempdir, "repset.txt")) as repset:
            rep_ids = [rep_id.strip("\n") for rep_id in repset.readlines()]

    if (max_size is not None) and (max_size < len(rep_ids)):
        rep_ids = rep_ids[:max_size]

    filter_fasta_by_ids(
        input_fasta=input_seqs, record_ids=rep_ids, output_fasta=outfile
    )


def reduce_database_redundancy(
    input_fasta: str,
    output_fasta: str = None,
    cdhit: bool = True,
    maxsize: int = None,
    cdhit_args: str = None,
) -> None:
    """
    Reduce redundancy of peptide datatabase.
    Runs cd-hit, if selected, additional arguments to cdhit
    may be passed as a string (cdhit_args).
    Runs repset to obtain a final database size no larger
    (number of sequences) than selected maxsize.
    If maxsize = None, repset is not run.
    """
    if (not cdhit) and (maxsize is None):
        warnings.warn("No reduction algorithm has been selected.")

    if output_fasta is None:
        output_fasta = set_default_output_path(input_fasta, tag="_reduced")

    with TemporaryFilePath() as tempaln, TemporaryFilePath() as tempfasta, TemporaryFilePath() as tempfasta2, TemporaryFilePath() as tempident:

        if cdhit:
            wrappers.run_cdhit(
                input_fasta=input_fasta,
                output_fasta=tempfasta,
                additional_args=cdhit_args,
            )
            os.remove(tempfasta + ".clstr")
        else:
            shutil.move(input_fasta, tempfasta)

        if maxsize is not None:
            wrappers.run_mafft(
                input_fasta=tempfasta,
                output_file=tempaln,
                n_threads=-1,
                parallel=True,
                additional_args="--retree 1 --maxiterate 0",
            )

            wrappers.get_percent_identity_from_msa(
                input_msa=tempaln, output_file=tempident
            )

            print("Finding representative sequences for reference database...")
            get_representative_set(
                input_seqs=tempfasta,
                input_pi=tempident,
                max_size=maxsize,
                outfile=tempfasta2,
            )
            shutil.move(tempfasta2, tempfasta)

        shutil.move(tempfasta, output_fasta)
