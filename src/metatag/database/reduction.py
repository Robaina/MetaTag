#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to reduce the size of the peptide-specific
reference database

Currently based on:
1. CD-HIT
2. Repset: https://onlinelibrary.wiley.com/doi/10.1002/prot.25461
"""

import logging
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

logger = logging.getLogger(__name__)


def get_representative_set(
    input_seqs: Path, input_pi: Path, max_size: int = None, output_file: Path = None
) -> None:
    """Runs repset.py to obtain a representative
       set of size equal to max_size (or smaller if less sequences than max_size)
       or an ordered list (by 'representativeness') of representative sequences
       if max_size set to None.

    Args:
        input_seqs (Path): path to input FASTA file containing sequences
        input_pi (Path): path to input file containing pairwise identity
        max_size (int, optional): maximum number of sequences in reduced database.
            Defaults to None.
        output_file (Path, optional): path to output file. Defaults to None.
    """
    input_seqs = Path(input_seqs).resolve()
    input_pi = Path(input_pi).resolve()
    repset_exe = Path(__file__).parent.parent / "vendor" / "repset_min.py"

    if output_file is None:
        output_file = set_default_output_path(input_seqs, tag="_repset")
    else:
        output_file = Path(output_file).resolve()

    with TemporaryDirectoryPath() as tempdir:
        cmd_str = (
            f"python {repset_exe} --seqs {input_seqs} --pi {input_pi} "
            f"--outdir {tempdir} --size {max_size}"
        )
        terminal_execute(cmd_str, suppress_shell_output=False)

        with open(Path(tempdir) / "repset.txt") as repset:
            rep_ids = [rep_id.strip("\n") for rep_id in repset.readlines()]

    if (max_size is not None) and (max_size < len(rep_ids)):
        rep_ids = rep_ids[:max_size]

    filter_fasta_by_ids(
        input_fasta=input_seqs, record_ids=rep_ids, output_fasta=output_file
    )


def reduce_database_redundancy(
    input_fasta: Path,
    output_fasta: Path = None,
    cdhit: bool = True,
    maxsize: int = None,
    cdhit_args: str = None,
) -> None:
    """Reduce redundancy of peptide datatabase.
       Runs cd-hit, if selected, additional arguments to cdhit
       may be passed as a string (cdhit_args).
       Runs repset to obtain a final database size no larger
       (number of sequences) than selected maxsize.
       If maxsize = None, repset is not run.

    Args:
        input_fasta (Path): path to input FASTA file.
        output_fasta (Path, optional): path to output, reduced fasta.
            Defaults to None.
        cdhit (bool, optional): whether to use CD-HIT alongside repset. Defaults to True.
        maxsize (int, optional): maximum number of sequences in final database. Defaults to None.
        cdhit_args (str, optional): additional arguments to CD-HIT. Defaults to None.
    """
    input_fasta = Path(input_fasta).resolve()
    if output_fasta is None:
        output_fasta = set_default_output_path(input_fasta, tag="_reduced")
    else:
        output_fasta = Path(output_fasta).resolve()
    if (not cdhit) and (maxsize is None):
        warnings.warn("No reduction algorithm has been selected.")

    with TemporaryFilePath() as tempaln, TemporaryFilePath() as tempfasta, TemporaryFilePath() as tempfasta2, TemporaryFilePath() as tempident:
        if cdhit:
            wrappers.run_cdhit(
                input_fasta=input_fasta,
                output_fasta=tempfasta,
                additional_args=cdhit_args,
            )
            Path(tempfasta.as_posix() + ".clstr").unlink(missing_ok=True)
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

            logger.info("Finding representative sequences for reference database...")
            get_representative_set(
                input_seqs=tempfasta,
                input_pi=tempident,
                max_size=maxsize,
                output_file=tempfasta2,
            )
            shutil.move(tempfasta2, tempfasta)

        shutil.move(tempfasta, output_fasta)
