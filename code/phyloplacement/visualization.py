"""
Functions to help visualize phylogenetic trees with and
without short read placement.
"""

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
            'Tree algorithm not found. Valid algorithm are: iqtree and fasttree'
            )