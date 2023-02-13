#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to process Paoli et al. 2019 data
https://doi.org/10.1101/2021.03.24.436479
"""

from metatag.utils import list_tar_dir


def get_genome_ids(paoli_tar: str) -> list:
    """
    Get Paoli database genome IDs from tar file
    """
    return [
        file.split("aPaoli/")[1].split(".fasta")[0]
        for file in list_tar_dir(paoli_tar)[1:]
    ]
