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

import phyloplacement.wrappers as wrappers
from phyloplacement.utils import (
    setDefaultOutputPath,
    terminalExecute,
    TemporaryFilePath,
    TemporaryDirectoryPath,
)
from phyloplacement.database.manipulation import filterFASTAbyIDs


def getRepresentativeSet(
    input_seqs: str, input_PI: str, max_size: int = None, outfile: str = None
) -> None:
    """
    Runs repset.py to obtain a representative
    set of size equal to max_size (or smaller if less sequences than max_size)
    or an ordered list (by 'representativeness') of representative sequences
    if max_size set to None.
    """
    input_seqs = os.path.abspath(input_seqs)
    input_PI = os.path.abspath(input_PI)
    repset_exe = os.path.abspath("code/vendor/repset_min.py")

    if outfile is None:
        outfile = setDefaultOutputPath(input_seqs, tag="_repset")

    with TemporaryDirectoryPath() as tempdir:
        cmd_str = (
            f"python {repset_exe} --seqs {input_seqs} --pi {input_PI} "
            f"--outdir {tempdir} --size {max_size}"
        )
        terminalExecute(cmd_str, suppress_shell_output=True)

        with open(os.path.join(tempdir, "repset.txt")) as repset:
            rep_ids = [rep_id.strip("\n") for rep_id in repset.readlines()]

    if (max_size is not None) and (max_size < len(rep_ids)):
        rep_ids = rep_ids[:max_size]

    filterFASTAbyIDs(input_fasta=input_seqs, record_ids=rep_ids, output_fasta=outfile)


def reduceDatabaseRedundancy(
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
        output_fasta = setDefaultOutputPath(input_fasta, tag="_reduced")

    with TemporaryFilePath() as tempaln, TemporaryFilePath() as tempfasta, TemporaryFilePath() as tempfasta2, TemporaryFilePath() as tempident:

        if cdhit:
            wrappers.runCDHIT(
                input_fasta=input_fasta,
                output_fasta=tempfasta,
                additional_args=cdhit_args,
            )
            os.remove(tempfasta + ".clstr")
        else:
            shutil.move(input_fasta, tempfasta)

        if maxsize is not None:
            wrappers.runMAFFT(
                input_fasta=tempfasta,
                output_file=tempaln,
                n_threads=-1,
                parallel=True,
                additional_args="--retree 1 --maxiterate 0",
            )

            wrappers.getPercentIdentityFromMSA(input_msa=tempaln, output_file=tempident)

            print("Finding representative sequences for reference database...")
            getRepresentativeSet(
                input_seqs=tempfasta,
                input_PI=tempident,
                max_size=maxsize,
                outfile=tempfasta2,
            )
            shutil.move(tempfasta2, tempfasta)

        shutil.move(tempfasta, output_fasta)
