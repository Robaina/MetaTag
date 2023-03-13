#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Pre-set pipelines to infer trees and place and label env sequences
"""

from __future__ import annotations

from pathlib import Path

from metatag.scripts import (
    buildtree,
    countplacements,
    labelplacements,
    makedatabase,
    placesequences,
    preprocess,
    relabeltree,
)


class CommandArgs:
    """Base class to hold command line arguments."""

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


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
        self.input_database = Path(input_database).resolve()
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
        if isinstance(hmmsearch_args, str) and (not hmmsearch_args.startswith(" ")):
            hmmsearch_args = " " + hmmsearch_args
        self.hmmsearch_args = hmmsearch_args
        self.msa_method = msa_method
        self.tree_method = tree_method
        self.tree_model = tree_model
        self.output_directory = Path(output_directory).resolve()

        self.out_cleaned_database = Path(
            self.output_directory / f"{self.input_database.stem}_cleaned.faa"
        ).resolve()
        self.out_reference_database = self.output_directory / "ref_database.faa"
        self.out_reference_tree = self.output_directory / "ref_database.newick"
        self.out_reference_alignment = self.output_directory / "ref_database.faln"
        self.out_reference_labels = (
            self.output_directory / "ref_database_id_dict.pickle"
        )
        self.out_tree_model = None

    def run(self) -> None:
        """Run pipeline to build reference tree."""

        preprocess_args = CommandArgs(
            data=self.input_database.as_posix(),
            outfile=self.out_cleaned_database.as_posix(),
            translate=True,
            dna=False,
            export_dup=False,
            relabel=False,
            idprefix=None,
            duplicate_method="seqkit",
        )
        preprocess.run(preprocess_args)

        database_args = CommandArgs(
            data=self.out_cleaned_database.as_posix(),
            outdir=self.output_directory.as_posix(),
            hmms=self.hmms,
            prefix="",
            relabel=True,
            nocdhit=False,
            noduplicates=False,
            relabel_prefixes=self.relabel_prefixes,
            maxsizes=self.max_hmm_reference_sizes,
            minseqlength=self.minimum_sequence_length,
            maxseqlength=self.maximum_sequence_length,
            hmmsearch_args=self.hmmsearch_args,
        )
        makedatabase.run(database_args)

        tree_args = CommandArgs(
            data=self.out_reference_database.as_posix(),
            outdir=self.output_directory.as_posix(),
            msa_method=self.msa_method,
            tree_model=self.tree_model,
            tree_method=self.tree_method,
        )
        buildtree.run(tree_args)

        relabel_args = CommandArgs(
            tree=self.out_reference_tree.as_posix(),
            labels=[self.out_reference_labels.as_posix()],
            labelprefixes=None,
            taxonomy=True,
            aln=None,
            outdir=None,
            taxofile=None,
        )
        relabeltree.run(relabel_args)


class QueryLabeller:
    """
    Place queries onto reference tree and assign function and taxonomy
    """

    def __init__(
        self,
        input_query: Path,
        reference_alignment: Path,
        reference_tree: Path,
        reference_labels: Path,
        tree_model: str,
        tree_clusters: Path,
        tree_cluster_scores: Path,
        tree_cluster_score_threshold: float = None,
        alignment_method: str = "papara",
        output_directory: Path = None,
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
            reference_labels (Path, optional): _description_.
            alignment_method (str, optional): _description_. Defaults to "papara".
            output_directory (Path, optional): _description_. Defaults to None.
            maximum_placement_distance (float, optional): _description_. Defaults to 1.0.
            distance_measure (str, optional): _description_. Defaults to "pendant_diameter_ratio".
            minimum_placement_lwr (float, optional): _description_. Defaults to 0.8.
        """
        self.input_query = Path(input_query).resolve()
        self.reference_alignment = Path(reference_alignment).resolve()
        self.reference_tree = Path(reference_tree).resolve()
        self.tree_model = (
            Path(tree_model).resolve() if Path(tree_model).is_file() else tree_model
        )
        self.tree_clusters = Path(tree_clusters).resolve()
        self.tree_cluster_scores = Path(tree_cluster_scores).resolve()
        self.tree_cluster_score_threshold = tree_cluster_score_threshold
        self.alignment_method = alignment_method
        self.reference_labels = reference_labels
        if output_directory is None:
            self.output_directory = self.input_query.parent / "placement_results"
            self.output_directory.mkdir(exist_ok=False)
        else:
            self.output_directory = Path(output_directory).resolve()
        self.maximum_placement_distance = maximum_placement_distance
        self.distance_measure = distance_measure
        self.minimum_placement_lwr = minimum_placement_lwr

        self.out_cleaned_query = Path(
            self.output_directory / "query_cleaned.faa"
        ).resolve()
        self.out_query_labels = self.output_directory / "query_cleaned_id_dict.pickle"

        self.place_outdir = self.output_directory / "placements"
        self.place_outdir.mkdir(parents=True, exist_ok=True)

        self.assign_outdir = self.output_directory / "assignments"
        self.assign_outdir.mkdir(parents=True, exist_ok=True)

        self.count_outdir = self.output_directory / "counts"
        self.count_outdir.mkdir(parents=True, exist_ok=True)

        self.jplace = self.filtered_jplace = self.place_outdir / "epa_result.jplace"
        if (self.maximum_placement_distance is None) and (
            self.minimum_placement_lwr is None
        ):
            self.filtered_jplace = self.place_outdir / "epa_result.jplace"
        elif (self.maximum_placement_distance is None) and (
            self.minimum_placement_lwr is not None
        ):
            self.filtered_jplace = self.place_outdir / "epa_result_lwr_filtered.jplace"
        elif (self.maximum_placement_distance is not None) and (
            self.minimum_placement_lwr is None
        ):
            self.filtered_jplace = (
                self.place_outdir / "epa_result_distance_filtered.jplace"
            )
        else:
            self.filtered_jplace = (
                self.place_outdir / "epa_result_distance_filtered_lwr_filtered.jplace"
            )

        self.out_placements_tree = self.place_outdir / "epa_result.newick"
        self.out_taxtable = self.assign_outdir / "placed_tax_assignments.tsv"

    def run(self) -> None:
        """Run pipeline to annotate query sequences through evolutionary placement."""
        preprocess_args = CommandArgs(
            data=self.input_query.as_posix(),
            outfile=self.out_cleaned_query.as_posix(),
            translate=True,
            dna=False,
            relabel=True,
            idprefix="query_",
            export_dup=True,
            duplicate_method="seqkit",
        )
        preprocess.run(preprocess_args)

        place_args = CommandArgs(
            aln=self.reference_alignment.as_posix(),
            tree=self.reference_tree.as_posix(),
            query=self.out_cleaned_query.as_posix(),
            outdir=self.place_outdir.as_posix(),
            aln_method=self.alignment_method,
            tree_model=self.tree_model,
        )
        placesequences.run(place_args)

        assign_args = CommandArgs(
            jplace=self.jplace.as_posix(),
            labels=self.reference_labels,
            query_labels=[self.out_query_labels.as_posix()],
            ref_clusters=self.tree_clusters.as_posix(),
            ref_cluster_scores=self.tree_cluster_scores.as_posix(),
            outgroup=None,
            prefix="placed_tax_",
            outdir=self.assign_outdir.as_posix(),
            max_distance=self.maximum_placement_distance,
            distance_measure=self.distance_measure,
            minimum_lwr=self.minimum_placement_lwr,
            duplicated_query_ids=None,
            taxofile=None,
        )
        labelplacements.run(assign_args)

        count_args = CommandArgs(
            taxtable=self.out_taxtable,
            taxlevels=["genus", "family", "order", "class", "phylum"],
            cluster_ids=None,
            score_threshold=self.tree_cluster_score_threshold,
            outdir=self.count_outdir.as_posix(),
            outprefix=None,
            export_right_queries=True,
        )
        countplacements.run(count_args)

        relabel_args = CommandArgs(
            tree=self.out_placements_tree.as_posix(),
            labels=self.reference_labels + [self.out_query_labels.as_posix()],
            labelprefixes=["ref_", "query_"],
            taxonomy=True,
            aln=None,
            outdir=None,
            taxofile=None,
        )
        relabeltree.run(relabel_args)
