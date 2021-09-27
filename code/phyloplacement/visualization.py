"""
Functions to help visualize phylogenetic trees with and
without short read placement.
"""
import os
import webbrowser
from .utils import terminalExecute


def reformatIQTreeForiTOL(input_tree: str):
    """
    NOTES:

    script: itoliqtree.py
    """
    pass

def reformatFastTreeForiTOL(input_tree: str):
    """
    NOTES:

    script: itolfasttree.py
    """
    pass

def reformatTreeForiTOL(input_tree: str, tree_algorithm: str) -> None:
    """
    Reformat tree data to suit iTOL input requirements
    """
    if 'fast' in tree_algorithm.lower():
        reformatFastTreeForiTOL(input_tree)
    elif 'iq' in tree_algorithm.lower():
        reformatIQTreeForiTOL(input_tree)
    else:
        raise ValueError(
            'Tree algorithm not found. Valid algorithms are: iqtree and fasttree'
            )

def plotTreeInBrowser(input_newick: str, output_dir: str = None,
                      feature_metadata: str = None) -> None:
    """
    Runs empress tree-plot
    feature_metadata: path to tsv file containing feature metadata
    """
    if output_dir is None:
        output_dir = os.path.join(
            os.path.dirname(input_newick), 'empress-viz'
            )
    if feature_metadata is not None:
        meta_str = f'-fm {feature_metadata}'
    else:
        meta_str = ''
    cmd_str = f'empress tree-plot -t {input_newick} -o {output_dir} {meta_str}'
    terminalExecute(cmd_str, suppress_output=True)
    webbrowser.open_new_tab(os.path.join(
        output_dir, 'empress.html'
    ))