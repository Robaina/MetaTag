#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Integration tests for MetaTag pipeline
"""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from metatag.pipelines import ReferenceTreeBuilder, QueryLabeller

tests_dir = Path(__file__).parent


class TestReferenceTreeBuilder(unittest.TestCase):
    def test_run(self):
        pass

class TestQueryLabeller(unittest.TestCase):
    def test_run(self):
        pass

class TestPipeline(unittest.TestCase):
    def test_pipeline(self):
        with TemporaryDirectory() as tempdir:

            tree_builder = ReferenceTreeBuilder(
                input_database=tests_dir / "test_data" / "database",
                hmms=[
                (tests_dir / "test_data" / "TIGR01287.1.HMM").as_posix(),
                (tests_dir / "test_data" / "TIGR02016.1.HMM").as_posix()
                ],
                maximum_hmm_reference_sizes=[20, 5],
                relabel_prefixes=["ref_", "out_"],
                relabel=True,
                remove_duplicates=True,
                hmmsearch_args="None, --cut_ga",
                output_directory=tempdir,
                msa_method="muscle",
                tree_method="fasttree",
                tree_model="iqtest"
            )
            tree_builder.run()
            labeller = QueryLabeller(
                input_query=tests_dir / "test_data" / "query.faa",
                reference_alignment=tree_builder.out_reference_alignment,
                reference_tree=tree_builder.out_reference_tree,
                reference_labels=[(Path(tempdir) / "ref_database_id_dict.pickle").as_posix()],
                tree_model="JTT",
                tree_clusters=tests_dir / "test_data" / "clusters.tsv",
                tree_cluster_scores=tests_dir / "test_data" / "cluster_scores.tsv",
                tree_cluster_score_threshold=0.6,
                alignment_method="papara",
                output_directory=tempdir,
                maximum_placement_distance=1.0,
                distance_measure="pendant_diameter_ratio",
                minimum_placement_lwr=0.8
            )
            labeller.run()
            family_counts = pd.read_csv(Path(tempdir) / "counts" / "placed_family_counts.tsv", sep="\t")
        self.assertGreater(
            len(family_counts.family.values),
            0,
            "Failed to retrieve correct taxonomy labels")



if __name__ == "__main__":
    unittest.main()