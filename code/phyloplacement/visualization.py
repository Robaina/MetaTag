#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to visualize phylogenetic trees with and
without short read placement.
"""

import os
import webbrowser
import pandas as pd

from phyloplacement.utils import terminal_execute


def make_feature_metadata_table(
    label_dict: dict, output_tsv: str, original_labels: bool = True
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
    input_tree: str, output_dir: str = None, feature_metadata: str = None
) -> None:
    """
    Runs empress tree-plot
    input_tree: tree in newick format
    feature_metadata: path to tsv file containing feature metadata
    empress: https://github.com/biocore/empress
    """
    if output_dir is None:
        output_dir = os.path.join(os.path.dirname(input_tree), "empress-plot")
    if feature_metadata is not None:
        meta_str = f"-fm {feature_metadata}"
    else:
        meta_str = ""
    cmd_str = f"empress tree-plot -t {input_tree} -o {output_dir} {meta_str}"
    terminal_execute(cmd_str, suppress_shell_output=False)
    webbrowser.open_new_tab(os.path.join(output_dir, "empress.html"))
