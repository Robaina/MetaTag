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
        maximum_hmm_reference_sizes: list[int] = None,
        minimum_sequence_length: int = 10,
        maximum_sequence_length: int = 1000,
        output_directory: Path = None,
        translate: bool = False,
        relabel: bool = True,
        remove_duplicates: bool = True,
        relabel_prefixes: list[str] = None,
        hmmsearch_args: str = None,
        tree_method: str = "fasttree",
        tree_model: str = "iqtest",
        msa_method: str = "muscle",
    ):
        """Reconstruct reference phylogenetic tree from sequence database and hmms

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
            tree_model (str, optional): choose method to select substitution model: "iqtest", "modeltest",
                or a valid substitution model name (compatible with EPA-ng). Defaults to "iqtest".
            msa_method (str, optional): choose msa method for reference database: "muscle" or "mafft".
                Defaults to "muscle".
        """
        self.input_database = Path(input_database)
        self.hmms = hmms
        if maximum_hmm_reference_sizes is None:
            self.max_hmm_reference_sizes = [1000 for _ in hmms]
        else:
            self.max_hmm_reference_sizes = maximum_hmm_reference_sizes
        self.minimum_sequence_length = minimum_sequence_length
        self.maximum_sequence_length = maximum_sequence_length
        self.translate = translate
        self.relabel = relabel
        self.remove_duplicates = remove_duplicates
        self.relabel_prefixes = relabel_prefixes
        self.hmmsearch_args = hmmsearch_args
        self.msa_method = msa_method
        self.tree_method = tree_method
        self.tree_model = tree_model
        self.output_directory = Path(output_directory)

        self.out_cleaned_database = Path(
            self.output_directory / f"{self.input_database.stem}_cleaned.faa"
        )
        self.out_reference_database = self.output_directory / "ref_database.faa"
        self.out_reference_tree = self.output_directory / "ref_database.newick"
        self.out_reference_alignment = self.output_directory / "ref_database.faln"
        self.out_reference_labels = (
            self.output_directory / "ref_database_id_dict.pickle"
        )
        self.out_tree_model = None

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
                self.out_cleaned_database,
                "--translate",
            ],
        )

        self.call_subcommand(
            "database",
            [
                "--in",
                self.out_cleaned_database,
                "--outdir",
                self.output_directory,
                "--hmms",
                self.hmms,
                "--max_sizes",
                self.max_hmm_reference_sizes,
                "--min_seq_length",
                self.minimum_sequence_length,
                "--max_seq_length",
                self.maximum_sequence_length,
                "--relabel_prefixes",
                self.relabel_prefixes,
                "--hmmsearch_args",
                self.hmmsearch_args,
            ],
        )

        self.call_subcommand(
            "tree",
            [
                "--in",
                self.out_reference_database,
                "--outdir",
                self.output_directory,
                "--msa_method",
                self.msa_method,
                "--tree_model",
                self.tree_model,
                "--tree_method",
                self.tree_method,
            ],
        )

        self.call_subcommand(
            "relabel",
            [
                "--tree",
                self.out_reference_tree,
                "--labels",
                self.out_reference_labels,
                "--taxonomy",
            ],
        )


class QueryLabeller:
    """
    Place queries onto reference tree and assign function and taxonomy
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
    ):
        """_summary_

        Args:
            input_query (Path): _description_
            reference_alignment (Path): _description_
            reference_tree (Path): _description_
            tree_model (str): _description_
            tree_clusters (Path): _description_
            tree_cluster_scores (Path): _description_
            alignment_method (str, optional): _description_. Defaults to "papara".
            output_directory (Path, optional): _description_. Defaults to None.
            reference_labels (Path, optional): _description_. Defaults to None.
            maximum_placement_distance (float, optional): _description_. Defaults to 1.0.
            distance_measure (str, optional): _description_. Defaults to "pendant_diameter_ratio".
            minimum_placement_lwr (float, optional): _description_. Defaults to 0.8.
        """
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

        self.out_cleaned_query = Path(self.output_directory / "query_cleaned.faa")
        self.out_query_labels = self.output_directory / "query_cleaned_id_dict.pickle"

        self.place_outdir = self.output_directory / "placements"
        self.place_outdir.mkdir(parents=True, exist_ok=True)

        self.assign_outdir = self.output_directory / "assignments"
        self.assign_outdir.mkdir(parents=True, exist_ok=True)

        self.count_outdir = self.output_directory / "counts"
        self.count_outdir.mkdir(parents=True, exist_ok=True)

        if (self.maximum_placement_distance is None) and (
            self.minimum_placement_lwr is None
        ):
            self.out_jplace = self.place_outdir / "epa_result.jplace"
        elif (self.maximum_placement_distance is None) and (
            self.minimum_placement_lwr is not None
        ):
            self.out_jplace = self.place_outdir / "epa_result_lwr_filtered.jplace"
        elif (self.maximum_placement_distance is not None) and (
            self.minimum_placement_lwr is None
        ):
            self.out_jplace = self.place_outdir / "epa_result_distance_filtered.jplace"
        else:
            self.out_jplace = (
                self.place_outdir / "epa_result_distance_filtered_lwr_filtered.jplace"
            )

        self.out_placements_tree = self.place_outdir / "epa_result.newick"
        self.out_taxtable = self.assign_outdir / "placed_taxonomy_assignments.tsv"

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
                self.place_outdir,
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
                self.assign_outdir,
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
                self.count_outdir,
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
