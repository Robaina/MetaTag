#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to visualize phylogenetic trees with and
without short read placement.
"""

import webbrowser
from pathlib import Path

import pandas as pd

from metatag.utils import terminal_execute


def make_feature_metadata_table(
    label_dict: dict, output_tsv: Path, original_labels: bool = True
) -> None:
    """Construct feature metadata tsv classifiying reference and
    query sequences for empress
    https://github.com/biocore/empress/issues/548

    Args:
        label_dict (dict): dictionary of sequence short IDs and labels
        output_tsv (Path): path to output tsv file
        original_labels (bool, optional): whether to include original
            long labels in tree. Defaults to True.
    """
    feature_dict = {
        seq_name if original_labels else seq_id: "ref" if "ref_" in seq_id else "query"
        for seq_id, seq_name in label_dict.items()
    }
    df = pd.DataFrame.from_dict(feature_dict, orient="index", columns=["Sequence type"])
    df.index.name = "Feature ID"
    df.to_csv(output_tsv, sep="\t")


def make_tree_html(
    input_tree: Path, output_dir: Path = None, feature_metadata: Path = None
) -> None:
    """Runs empress tree-plot
    empress:  https://github.com/biocore/empress

    Args:
    input_tree (Path): path to tree in newick format
    output_dir (Path, optional): path to output directory. Defaults to None.
    feature_metadata (Path, optional): path to fieature metadata
        table as output by make_feature_metadata_table. Defaults to None.
    """
    input_tree = Path(input_tree).resolve()
    if output_dir is None:
        output_dir = input_tree.parent / "empress-plot"
    else:
        output_dir = Path(output_dir).resolve()
    if feature_metadata is not None:
        meta_str = f"-fm {feature_metadata}"
    else:
        meta_str = ""
    cmd_str = f"empress tree-plot -t {input_tree} {meta_str} -o {output_dir}"
    terminal_execute(cmd_str, suppress_shell_output=False)


def plot_tree_in_browser(
    input_tree: Path, output_dir: Path = None, feature_metadata: Path = None
) -> None:
    """Runs empress tree-plot and opens generated html in browser
    empress: https://github.com/biocore/empress

    Args:
        input_tree (Path): path to tree in newick format
        output_dir (Path, optional): path to output directory. Defaults to None.
        feature_metadata (Path, optional): path to fieature metadata
            table as output by make_feature_metadata_table. Defaults to None.
    """
    make_tree_html(input_tree, output_dir, feature_metadata)
    webbrowser.open_new_tab((output_dir / "empress.html").as_posix())
