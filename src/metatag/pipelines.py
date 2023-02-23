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
    """
    Place queries onto reference tree and assign function and taxonomy
    TODO:
    1) Subdirectories in output folder
    2) Enable re-run for count without doing alignnment and placement
    3) Track output files
    """

    def __init__(
        self,
        input_query: Path,
        reference_alignment: Path,
        reference_tree: Path,
        tree_model: str,
        tree_clusters: Path,
        tree_cluster_scores: Path,
        alignment_method: str = "papara",
        output_directory: Path = None,
        reference_labels: Path = None,
        maximum_placement_distance: float = 1.0,
        distance_measure: str = "pendant_diameter_ratio",
        minimum_placement_lwr: float = 0.8,
    ) -> None:
        self.input_query = Path(input_query)
        self.reference_alignment = Path(reference_alignment)
        self.reference_tree = Path(reference_tree)
        self.tree_model = Path(tree_model) if Path(tree_model).is_file() else tree_model
        self.tree_clusters = Path(tree_clusters)
        self.tree_cluster_scores = Path(tree_cluster_scores)
        self.aligment_method = alignment_method
        if output_directory is None:
            self.output_directory = self.input_query.parent / "placement_results"
            self.output_directory.mkdir(exist_ok=False)
        else:
            self.output_directory = Path(output_directory)
        if reference_labels is not None:
            self.reference_labels = Path(reference_labels)
        self.maximum_placement_distance = maximum_placement_distance
        self.distance_measure = distance_measure
        self.minimum_placement_lwr = minimum_placement_lwr

        self.out_cleaned_query = Path(
            self.output_directory / f"{self.input_query.stem}_cleaned.faa"
        )
        self.out_query_labels = Path()
        self.out_jplace = Path()
        self.out_taxtable = Path()
        self.out_placements_tree = Path()

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
                self.out_cleaned_query,
                "--translate",
            ],
        )

        self.call_subcommand(
            "place",
            [
                "--aln",
                self.reference_alignment,
                "--tree",
                self.reference_tree,
                "--query",
                self.out_cleaned_query,
                "--outdir",
                self.output_directory,
                "--aln_method",
                self.aligment_method,
                "--tree_model",
                self.tree_model,
            ],
        )

        self.call_subcommand(
            "assign",
            [
                "--jplace",
                self.out_jplace,
                "--labels",
                self.reference_labels,
                "--query_labels",
                self.out_query_labels,
                "--ref_clusters",
                self.tree_clusters,
                "--ref_cluster_scores",
                self.tree_cluster_scores,
                "--prefix",
                "placed_tax_",
                "--outdir",
                self.output_directory,
                "--max_placement_distance",
                self.maximum_placement_distance,
                "--distance_measure",
                self.distance_measure,
                "--min_placement_lwr",
                self.minimum_placement_lwr,
            ],
        )

        self.call_subcommand(
            "count",
            [
                "--taxtable",
                self.out_taxtable,
                "--taxlevels",
                "genus",
                "family",
                "order",
                "class",
                "phylum",
                "--cluster_ids",
                self.tree_cluster_ids,
                "--score_threshold",
                self.tree_cluster_scores,
                "--outdir",
                self.output_directory,
                "--export_right_queries",
            ],
        )

        self.call_subcommand(
            "relabel",
            [
                "--tree",
                self.out_placements_tree,
                "--labels",
                self.reference_labels,
                self.out_query_labels,
                "--label_prefixes",
                "ref_",
                "query_",
                "--taxonomy",
            ],
        )
