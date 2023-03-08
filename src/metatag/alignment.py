#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to perform multiple sequence alignments
"""

import logging
import sys
from pathlib import Path

import metatag.wrappers as wrappers
from metatag.database.manipulation import (
    convert_fasta_aln_to_phylip,
    convert_phylip_to_fasta_aln,
    convert_stockholm_to_fasta_aln,
)
from metatag.utils import TemporaryFilePath, set_default_output_path

logger = logging.getLogger(__name__)


def align_peptides(
    input_fasta: Path,
    method: str = "muscle",
    output_file: Path = None,
    additional_args: str = None,
):
    """
    Perform MSA on reference peptide sequences.
    Outputs in format fasta.aln
    """
    if output_file is None:
        output_file = set_default_output_path(
            input_fasta, extension=".fasta.aln", only_filename=True
        )
    input_fasta = Path(input_fasta)
    if method.lower() in "muscle":
        wrappers.run_muscle(
            input_fasta=input_fasta,
            output_file=output_file,
            additional_args=additional_args,
        )
    elif method.lower() in "mafft":
        wrappers.run_mafft(
            input_fasta=input_fasta,
            output_file=output_file,
            additional_args=additional_args,
        )
    else:
        logger.error("Invalid method. Valid methods: muscle or mafft")
        sys.exit(1)


def align_short_reads_to_reference_msa(
    ref_msa: Path,
    query_seqs: Path,
    method: str = "papara",
    tree_nwk: Path = None,
    output_dir: Path = None,
) -> None:
    """
    Align short read query sequences to reference MSA (fasta format).
    Outputs fasta msa alignment between query and reference sequences
    """
    ref_msa = Path(ref_msa)
    query_seqs = Path(query_seqs)
    tree_nwk = Path(tree_nwk)

    if output_dir is None:
        output_dir = set_default_output_path(ref_msa, only_dirname=True)
    output_hmm = output_dir / set_default_output_path(
        ref_msa, extension=".hmm", only_filename=True
    )
    output_aln_seqs = output_dir / set_default_output_path(
        query_seqs, extension=".faln", only_filename=True
    )

    if method.lower() in "hmmalign":
        with TemporaryFilePath() as temp_file_path:
            wrappers.run_hmmbuild(input_aln=ref_msa, output_hmm=output_hmm)

            wrappers.run_hmmalign(
                input_aln=ref_msa,
                input_hmm=output_hmm,
                input_seqs=query_seqs,
                output_aln_seqs=temp_file_path,
            )
            convert_stockholm_to_fasta_aln(
                input_stockholm=temp_file_path, output_fasta=output_aln_seqs
            )
    elif method.lower() in "papara":
        with TemporaryFilePath() as temp_phy_path, TemporaryFilePath() as temp_aln_path:
            convert_fasta_aln_to_phylip(
                input_fasta_aln=ref_msa, output_phylip=temp_phy_path
            )
            wrappers.run_papara(
                tree_nwk=tree_nwk,
                msa_phy=temp_phy_path,
                query_fasta=query_seqs,
                output_aln=temp_aln_path,
            )
            convert_phylip_to_fasta_aln(
                input_phylip=temp_aln_path, output_file=output_aln_seqs
            )
    else:
        logger.error("Alignment method not implemented")
        sys.exit(1)
