#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Pre-set pipelines to infer trees and place and label env sequences
"""

from __future__ import annotations

from pathlib import Path

from cli import MetaTag


class ReferenceTreeBuilder:
    """_summary_"""

    def __init__(
        self,
        input_database: Path,
        hmms: list[Path],
        max_hmm_reference_size: list[int] = None,
        min_sequence_length: int = 10,
        max_sequence_length: int = 1000,
        output_directory: Path = None,
        translate: bool = False,
        relabel: bool = True,
        remove_duplicates: bool = True,
        relabel_prefixes: list[str] = None,
        hmmsearch_args: str = None,
        tree_method: str = "fasttree",
    ) -> None:
        """_summary_

        Args:
            input_database (Path): _description_
            hmms (list[Path]): _description_
            max_hmm_reference_size (list[int], optional): _description_. Defaults to None.
            min_sequence_length (int, optional): _description_. Defaults to 10.
            max_sequence_length (int, optional): _description_. Defaults to 1000.
            output_directory (Path, optional): _description_. Defaults to None.
            translate (bool, optional): _description_. Defaults to False.
            relabel (bool, optional): _description_. Defaults to True.
            remove_duplicates (bool, optional): _description_. Defaults to True.
            relabel_prefixes (list[str], optional): _description_. Defaults to None.
            hmmsearch_args (str, optional): _description_. Defaults to None.
            tree_method (str, optional): _description_. Defaults to "fasttree".
        """
        if max_hmm_reference_size is None:
            max_hmm_reference_size = [1000 for _ in hmms]

        self.input_database = Path(input_database)
        self.output_directory = Path(output_directory)

    def call_subcommand(self, subcommand: str, args: list[str]) -> None:
        """_summary_

        Args:
            subcommand (str): _description_
            args (list[str]): _description_
        """
        MetaTag(subcommand, args, silent=True)._call_subcommand(subcommand)

    def run(self):
        """_summary_"""

        self.call_subcommand(
            "preprocess",
            [
                "--in",
                self.input_database,
                "--outfile",
                Path(self.output_directory / f"{self.input_database.stem}_cleaned.faa"),
                "--translate",
            ],
        )

        self.call_subcommand(
            "database",
            [
                "--in",
                "",
                "--outdir",
                "",
                "--hmms",
                "",
                "--max_sizes",
                "",
                "--min_seq_length",
                "",
                "--max_seq_length",
                "",
                "--relabel_prefixes",
                "",
                "--hmmsearch_args",
                "",
            ],
        )

        self.call_subcommand(
            "tree",
            [
                "--in",
                "",
                "--outdir",
                "",
                "--msa_method",
                "",
                "--tree_model",
                "",
                "--tree_method",
                "",
            ],
        )

        self.call_subcommand("relabel", ["--tree", "", "--labels", "", "--taxonomy"])


class QueryLabeller:
    """_summary_"""

    def __init__(self, input_query: Path) -> None:
        self.input_query = Path(input_query)

    def call_subcommand(self, subcommand: str, args: list[str]) -> None:
        """_summary_

        Args:
            subcommand (str): _description_
            args (list[str]): _description_
        """
        MetaTag(subcommand, args, silent=True)._call_subcommand(subcommand)

    def run(self):
        """_summary_"""
        self.call_subcommand(
            "preprocess",
            [
                "--in",
                self.input_query,
                "--outfile",
                Path(self.output_directory / f"{self.input_query.stem}_cleaned.faa"),
                "--translate",
            ],
        )

        self.call_subcommand(
            "place",
            [
                "--aln",
                "",
                "--tree",
                "",
                "--query",
                "",
                "--outdir",
                "",
                "--aln_method",
                "",
                "--tree_model",
                "",
            ],
        )

        self.call_subcommand(
            "assign",
            [
                "--jplace",
                "",
                "--labels",
                "",
                "--query_labels",
                "",
                "--ref_clusters",
                "",
                "--ref_cluster_scores",
                "",
                "--prefix",
                "placed_tax_",
                "--outdir",
                "",
                "--max_placement_distance",
                "1.0",
                "--distance_measure",
                "pendant_diameter_ratio",
                "--min_placement_lwr",
                "0.8",
            ],
        )

        self.call_subcommand(
            "count",
            [
                "--taxtable",
                "",
                "--taxlevels",
                "",
                "--cluster_ids",
                "",
                "--score_threshold",
                "",
                "--outdir",
                "",
                "--export_right_queries",
            ],
        )

        self.call_subcommand(
            "relabel",
            ["--tree", "", "--labels", "", "--label_prefixes", "", "--taxonomy"],
        )
