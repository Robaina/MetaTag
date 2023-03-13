#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to create peptide-specific sequence databases

1. Implement hmmr
2. Filter fasta files based on query sequences
"""

from __future__ import annotations

import logging
import sys
import tempfile
import warnings
from collections import defaultdict
from pathlib import Path

import pandas as pd
import pyfastx
from Bio import AlignIO, SearchIO

import metatag.wrappers as wrappers
from metatag.utils import set_default_output_path, terminal_execute

logger = logging.getLogger(__name__)


def filter_fasta_by_sequence_length(
    input_fasta: Path,
    min_length: int = None,
    max_length: int = None,
    output_fasta: Path = None,
) -> None:
    """
    Filter sequences by length in fasta file
    """
    if (min_length is None) and (max_length is None):
        warnings.warn("Missing boundary values for sequence length")
        return
    input_fasta = Path(input_fasta).resolve()
    with open(input_fasta, "r") as ffile:
        fasta = pyfastx.Fasta(ffile.name)
        record_ids = fasta.keys()
        if min_length is None:
            min_length = 0
        if max_length is not None:
            max_tag = str(max_length)
            record_ids.filter(record_ids >= min_length, record_ids <= max_length)
        else:
            max_tag = ""
            record_ids.filter(record_ids >= min_length)
        if output_fasta is None:
            output_fasta = set_default_output_path(
                input_fasta, f"_length_{min_length}_{max_tag}"
            )
        else:
            output_fasta = Path(output_fasta).resolve()
        if not record_ids:
            logger.error("No records found with given sequence length bounds")
            sys.exit(1)
        with open(output_fasta, "w") as fp:
            for record_id in record_ids:
                record_obj = fasta[record_id]
                fp.write(record_obj.raw)
    Path(input_fasta.as_posix() + ".fxi").unlink(missing_ok=True)


def parse_hmmsearch_output(hmmer_output: Path) -> pd.DataFrame:
    """
    Parse hmmsearch or hmmscan summary table output file
    """
    attribs = ["id", "bias", "bitscore", "description"]
    hits = defaultdict(list)
    with open(hmmer_output) as handle:
        for queryresult in SearchIO.parse(handle, "hmmer3-tab"):
            for hit in queryresult.hits:
                for attrib in attribs:
                    hits[attrib].append(getattr(hit, attrib))
    return pd.DataFrame.from_dict(hits)


def filter_fasta_by_ids(
    input_fasta: Path, record_ids: list, output_fasta: Path = None
) -> None:
    """
    Filter records in fasta file matching provided IDs
    """
    input_fasta = Path(input_fasta).resolve()
    if output_fasta is None:
        output_fasta = set_default_output_path(input_fasta, "_fitered")
    else:
        output_fasta = Path(output_fasta).resolve()
    index_file = Path(input_fasta.as_posix() + ".fxi").resolve()
    index_file.unlink(missing_ok=True)
    record_ids = set(record_ids)
    with tempfile.NamedTemporaryFile(mode="w+t") as tmp_ids:
        tmp_ids.writelines("\n".join(record_ids))
        tmp_ids.flush()
        tmp_ids_path = tmp_ids.name
        cmd_str = f"seqkit grep -i -f {tmp_ids_path} {input_fasta} -o {output_fasta}"
        terminal_execute(cmd_str, suppress_shell_output=True)


def filter_fasta_by_hmm(
    hmm_model: Path,
    input_fasta: Path,
    output_fasta: Path = None,
    hmmer_output: Path = None,
    method: str = "hmmsearch",
    additional_args: str = None,
) -> None:
    """
    Generate protein-specific database by filtering
    sequence database to only contain sequences
    corresponing to protein of interest

    @Arguments:
    additional_args: additional arguments to hmmsearch or hmmscan
    """
    hmm_model = Path(hmm_model).resolve()
    input_fasta = Path(input_fasta).resolve()
    if hmmer_output is None:
        hmmer_output = set_default_output_path(
            input_fasta, tag=f"_{hmm_model.stem}", extension=".txt"
        )
    else:
        hmmer_output = Path(hmmer_output).resolve()
    if output_fasta is None:
        output_fasta = set_default_output_path(
            input_fasta, tag=f"filtered_{hmm_model.stem}"
        )
    else:
        output_fasta = Path(output_fasta).resolve()

    logger.info("Running Hmmer...")
    wrappers.run_hmmsearch(
        hmm_model=hmm_model,
        input_fasta=input_fasta,
        output_file=hmmer_output,
        method=method,
        additional_args=additional_args,
    )
    logger.info("Parsing Hmmer output file...")
    hmmer_hits = parse_hmmsearch_output(hmmer_output)
    if not hmmer_hits.id.values.tolist():
        logger.error("No records found in database matching provided hmm")
        sys.exit(1)
    logger.info("Filtering Fasta...")
    filter_fasta_by_ids(
        input_fasta, record_ids=hmmer_hits.id.values, output_fasta=output_fasta
    )


def convert_fasta_aln_to_phylip(
    input_fasta_aln: Path, output_phylip: Path = None
) -> None:
    """
    Convert alignments in Fasta to Phylip.
    """
    input_fasta_aln = Path(input_fasta_aln).resolve()
    if output_phylip is None:
        output_phylip = set_default_output_path(input_fasta_aln, extension=".phylip")
    else:
        output_phylip = Path(output_phylip).resolve()
    with open(input_fasta_aln, "r") as input_handle, open(
        output_phylip, "w"
    ) as output_handle:
        alignments = AlignIO.parse(input_handle, "fasta")
        AlignIO.write(alignments, output_handle, "phylip-relaxed")


def convert_phylip_to_fasta_aln(input_phylip: Path, output_file: Path = None) -> None:
    """
    Convert alignments in Phylip to Fasta format
    """
    input_phylip = Path(input_phylip).resolve()
    if output_file is None:
        output_file = set_default_output_path(input_phylip, extension=".faln")
    else:
        output_file = Path(output_file).resolve()
    alignments = AlignIO.parse(input_phylip, "phylip-relaxed")
    AlignIO.write(alignments, output_file, "fasta")


def convert_stockholm_to_fasta_aln(
    input_stockholm: Path, output_fasta: Path = None
) -> None:
    """
    Convert alignment file in Stockholm format to fasta
    """
    input_stockholm = Path(input_stockholm).resolve()
    if output_fasta is None:
        output_fasta = set_default_output_path(input_stockholm, extension=".faln")
    else:
        output_fasta = Path(output_fasta).resolve()
    alignments = AlignIO.read(input_stockholm, "stockholm")
    AlignIO.write(alignments, output_fasta, "fasta")


def split_reference_from_query_alignments(
    ref_query_msa: Path,
    ref_ids: set = None,
    ref_prefix: str = None,
    output_dir: Path = None,
) -> None:
    """
    Separate reference sequences from query sequences in msa fasta file
    """
    ref_query_msa = Path(ref_query_msa).resolve()
    if output_dir is None:
        output_dir = ref_query_msa.parent
    else:
        output_dir = Path(output_dir).resolve()
    if (ref_ids is not None) and (ref_prefix is None):

        def is_reference(record_name):
            return record_name in ref_ids

    elif (ref_ids is None) and (ref_prefix is not None):

        def is_reference(record_name):
            return record_name.lower().startswith(ref_prefix.lower())

    else:
        logger.error("Provide either set of ref ids or ref prefix")
        sys.exit(1)
    out_ref_msa = set_default_output_path(ref_query_msa, tag="_ref_fraction")
    out_query_msa = set_default_output_path(ref_query_msa, tag="_query_fraction")

    with open(out_ref_msa, "w") as outref, open(out_query_msa, "w") as outquery, open(
        ref_query_msa, "r"
    ) as msa:
        fasta = pyfastx.Fasta(msa.name, build_index=False, full_name=True)
        for record_name, record_seq in fasta:
            if is_reference(record_name):
                outref.write(f">{record_name}\n{record_seq}\n")
            else:
                outquery.write(f">{record_name}\n{record_seq}\n")


def get_fasta_record_ids(fasta_file: Path) -> set:
    """
    Extract record ids from fasta
    """
    fasta_file = Path(fasta_file).resolve()
    with open(fasta_file, "r") as ffile:
        fasta = pyfastx.Fasta(ffile.name, full_name=True)
        fasta_keys = set(fasta.keys())
    return fasta_keys
