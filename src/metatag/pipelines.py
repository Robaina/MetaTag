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
        pass

    def run(self):
        """_summary_"""
        MetaTag.preprocess()
        MetaTag.database()
        MetaTag.tree()
        MetaTag.relabel()


class QueryLabeller:
    """_summary_"""

    def __init__(self) -> None:
        pass

    def run(self):
        """_summary_"""
        MetaTag.preprocess()
        MetaTag.place()
        MetaTag.assign()
        MetaTag.count()
