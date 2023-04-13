#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to parse sequence labels from different databases
"""

import re
import warnings


class LabelParser:
    def __init__(self, label: str) -> None:
        """Parse labels to extract genome ID and metadata

        Args:
            label (str): label to be parsed
        """
        self._label = label

    def extract_genome_id(self) -> str:
        """Extract genome ID from sequence label

        Returns:
            str: Genome ID
        """
        mmp_id = self.extract_mmp_id()
        taxid = self.extract_taxid()
        if mmp_id and taxid:
            warnings.warn("Label contains conflicting genome IDs")
            return ""
        elif mmp_id:
            genome_id = mmp_id
        elif taxid:
            genome_id = taxid
        else:
            genome_id = self.extract_oceanmicrobiome_id()
        return genome_id

    def extract_mmp_id(self) -> str:
        """Extract MMP ID from sequence label

        Returns:
            str: MMP ID
        """
        db_entry = re.compile("_MMP\d+")
        try:
            return re.search(db_entry, self._label).group(0).strip("_")
        except Exception:
            return ""

    def extract_taxid(self) -> str:
        """Extract NCBI taxon ID from sequence label

        Returns:
            str: NCBI taxon ID
        """
        db_entry = re.compile("(OX=)(\d+)")
        try:
            return f"taxid_{re.search(db_entry, self._label).group(2)}"
        except Exception:
            return ""

    def extract_oceanmicrobiome_id(self) -> str:
        """Extract OceanMicrobiome ID from sequence label

        Returns:
            str: OceanMicrobiome ID
        """
        return "_".join(self._label.split("_")[:4])
