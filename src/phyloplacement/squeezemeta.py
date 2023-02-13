#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to parse and modify SqueezeMeta's output table
"""
from pathlib import Path

import pandas as pd


class SqueezeMetaOutputParser:
    def __init__(self, squeeze_output: Path) -> None:
        """ """
        self._squeeze_path = squeeze_output
        self._squeeze_out = pd.read_csv(
            squeeze_output, sep="\t", header=1
        )  # .dropna(subset=["Tax"])
        with open(squeeze_output) as f:
            firstline = f.readline()
        self._metaheader = firstline

    @staticmethod
    def taxo_to_squeeze_format(taxopath: str) -> str:
        """
        Convert GTDB taxopath format to SqueezeMeta's format
        """
        taxolevels = {
            "d__": "k_",
            "p__": "p_",
            "c__": "c_",
            "o__": "o_",
            "f__": "f_",
            "g__": "g_",
            "s__": "s_",
        }
        for level_gtdb, level_squeeze in taxolevels.items():
            taxopath = taxopath.replace(level_gtdb, level_squeeze)
        return taxopath

    def replace_by_placement_taxopath(
        self, placement_assignments: Path, output_file: Path = None
    ) -> None:
        """
        Replace SqueezeMeta's taxopath by placement taxopath in right format
        """
        if output_file is None:
            output_file = (
                self._squeeze_path.parent
                / f"{self._squeeze_path.stem}_replaced{self._squeeze_path.suffix}"
            )
        placed_df = pd.read_csv(placement_assignments, sep="\t")
        # placed_df = placed_df[placed_df.taxopath != "Unspecified"]
        placed_df.taxopath = placed_df.taxopath.apply(
            lambda p: "" if p == "Unspecified" else p
        )
        filtered_squeeze = self._squeeze_out[
            self._squeeze_out["ORF ID"].isin(placed_df["query_name"])
        ].copy()
        for i, row in filtered_squeeze.iterrows():
            gtdb_tax = placed_df.loc[
                placed_df.query_name == row["ORF ID"], "taxopath"
            ].item()
            filtered_squeeze.loc[i, "Tax"] = self.taxo_to_squeeze_format(gtdb_tax)
        filtered_squeeze.to_csv(output_file, sep="\t", index=False)
        # Prepend meta header
        with open(output_file, "r+") as file:
            content = file.read()
            file.seek(0)
            file.write(self._metaheader + content)


class SqueezeMetaTaxonomyParser:
    def __init__(self, squeeze_output: Path) -> None:
        """ """
        self._squeeze_path = squeeze_output
        self._squeeze_out = pd.read_csv(
            squeeze_output, sep="\t", header=None, skiprows=1
        ).rename({0: "ORF ID", 1: "Tax"}, axis=1)
        with open(squeeze_output) as f:
            firstline = f.readline()
        self._metaheader = firstline

    @staticmethod
    def taxo_to_squeeze_format(taxopath: str) -> str:
        """
        Convert GTDB taxopath format to SqueezeMeta's format
        """
        taxolevels = {
            "d__": "k_",
            "p__": "p_",
            "c__": "c_",
            "o__": "o_",
            "f__": "f_",
            "g__": "g_",
            "s__": "s_",
        }
        for level_gtdb, level_squeeze in taxolevels.items():
            taxopath = taxopath.replace(level_gtdb, level_squeeze)
        return taxopath

    def replace_by_placement_taxopath(
        self, placement_assignments: Path, output_file: Path = None
    ) -> None:
        """
        Replace SqueezeMeta's taxopath by placement taxopath in right format
        """
        if output_file is None:
            output_file = (
                self._squeeze_path.parent
                / f"{self._squeeze_path.stem}_replaced{self._squeeze_path.suffix}"
            )
        placed_df = pd.read_csv(placement_assignments, sep="\t")
        placed_df.taxopath = placed_df.taxopath.apply(
            lambda p: "" if p == "Unspecified" else p
        )
        filtered_squeeze = self._squeeze_out[
            self._squeeze_out["ORF ID"].isin(placed_df["query_name"])
        ].copy()
        for i, row in filtered_squeeze.iterrows():
            gtdb_tax = placed_df.loc[
                placed_df.query_name == row["ORF ID"], "taxopath"
            ].item()
            filtered_squeeze.loc[i, "Tax"] = self.taxo_to_squeeze_format(gtdb_tax)
        filtered_squeeze.to_csv(output_file, sep="\t", index=False, header=False)
        # Prepend meta header
        with open(output_file, "r+") as file:
            content = file.read()
            file.seek(0)
            file.write(self._metaheader + content)
