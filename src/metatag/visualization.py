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
    """
    Construct feature metadata tsv classifiying reference and
    query sequences for empress
    https://github.com/biocore/empress/issues/548
    """
    feature_dict = {
        seq_name if original_labels else seq_id: "ref" if "ref_" in seq_id else "query"
        for seq_id, seq_name in label_dict.items()
    }
    df = pd.DataFrame.from_dict(feature_dict, orient="index", columns=["Sequence type"])
    df.index.name = "Feature ID"
    df.to_csv(output_tsv, sep="\t")


def plot_tree_in_browser(
    input_tree: Path, output_dir: Path = None, feature_metadata: Path = None
) -> None:
    """
    Runs empress tree-plot
    input_tree: tree in newick format
    feature_metadata: path to tsv file containing feature metadata
    empress: https://github.com/biocore/empress
    """
    input_tree = Path(input_tree)
    if output_dir is None:
        output_dir = input_tree.parent / "empress-plot"
    else:
        output_dir = Path(output_dir)
    if feature_metadata is not None:
        meta_str = f"-fm {feature_metadata}"
    else:
        meta_str = ""
    cmd_str = f"empress tree-plot -t {input_tree} -o {output_dir} {meta_str}"
    terminal_execute(cmd_str, suppress_shell_output=True)
    webbrowser.open_new_tab(output_dir / "empress.html")
