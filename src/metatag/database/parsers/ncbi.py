#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Obtain FASTA with gene sequences from list of NCBI accession ids
Dependencies: ncbi-acc-download
"""

import logging
import os

from Bio import SeqIO

from metatag.utils import terminal_execute

logger = logging.getLogger(__name__)


def download_gbk_from_ncbi(entry_ids: list, output_dir: str = None) -> None:
    """
    Download genbank files from NCBI from given list of entry IDs
    """
    for n, entry_id in enumerate(entry_ids):
        logger.info(f"Downloading entry: {entry_id} ({n + 1} / {len(entry_ids)})")
        outfasta = os.path.join(output_dir, f"{entry_id}.gbk")
        cmd_str = f"ncbi-acc-download -o {outfasta} {entry_id}"
        terminal_execute(cmd_str)


def get_protein_sequence_from_gbk(
    gbk: str, cds_keywords: dict, case_insensitive: bool = True
) -> dict:
    """
    Extract cds record matching keywords from gbk file.
      @Arguments:
      gbk: path to genbank file
      keywords is a dictionary in which keys correspond to gbk cds fields
      and values to keywords to find in each field. For instace,
      keywords = {
          'gene': ['ureC'],
          'product': ['urease', 'alpha']
      }
      case_insensity: whether or not to care for case when matching keywords
    """

    def contains_keywords(text: str, keywords: list) -> bool:
        if case_insensitive:
            return all([key.lower() in text.lower() for key in keywords])
        else:
            return all([key in text for key in keywords])

    def is_a_match(cds: dict) -> bool:
        return all(
            [
                field in cds.keys() and contains_keywords(cds[field][0], keywords)
                for field, keywords in cds_keywords.items()
            ]
        )

    gbk = list(SeqIO.parse(gbk, "genbank"))[0]
    cds_records = [f.qualifiers for f in gbk.features[1:] if "CDS" in f.type]
    return [cds for cds in cds_records if is_a_match(cds)]


def write_fasta_from_cds_qualifiers(records: dict, output_fasta: str = None) -> None:
    """
    Write FASTA file from dict of gbk record qualifiers (Ordered dict)
    """
    with open(output_fasta, "w") as file:
        for record_id, record in records.items():
            if record:
                record = record[0]
                ref_id = f'{record_id}_{record["protein_id"][0]}_{"_".join(record["product"][0].split())}'
                file.write(f">{ref_id}\n")
                file.write(f'{record["translation"][0]}\n')


def get_fasta_for_gene(
    gbk_dir: str,
    gene_keywords: dict,
    output_fasta: str = None,
    case_insensitive: bool = True,
) -> None:
    """
    Write FASTA from list of GenBank files and cds entry field keywords
    """
    gbk_dir = os.path.abspath(gbk_dir)
    records_dict = {
        gbk_file.split(".gbk")[0]: get_protein_sequence_from_gbk(
            os.path.join(gbk_dir, gbk_file),
            cds_keywords=gene_keywords,
            case_insensitive=case_insensitive,
        )
        for gbk_file in os.listdir(gbk_dir)
    }
    write_fasta_from_cds_qualifiers(records_dict, output_fasta=output_fasta)
