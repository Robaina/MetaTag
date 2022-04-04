#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to perform phylogenetic tree reconstructions and
query sequence placements onto trees
"""


from __future__ import annotations
import os

from Bio import Phylo
import pyfastx



class PhyloTree:
    """
    Methods to help defining clusters in phylo trees
    """
    def __init__(self, tree_path: str, tree_format: str = 'newick', bootstrap_threshold: float = None):
        self._tree = list(Phylo.parse(tree_path, tree_format))[0]
        if bootstrap_threshold is not None:
            self.collapsePoorQualityNodes(bootstrap_threshold)
        self.nameInternalNodes()
        self._tree_path = tree_path
        self._tree_format = tree_format

    def _scaleBootstrapValues(self):
        """
        Scale bootstrap values to be percentages if necessary.
        This is to deal with discrepancies between how fasttree and
        iqtree report bootstrap values (fractions vs percentages)
        """
        conf_values = [
            c.confidence for c in self._tree.find_clades() 
            if c.confidence is not None
            ]
        if conf_values:
            max_value = max(conf_values)
            if max_value < 2:
                for clade in self._tree.find_clades():
                    if clade.confidence is not None:
                        clade.confidence *= 100
        else:
            raise ValueError('Tree does not contain confidence values. Change tree or set bootstrap_threshold to None')

    def collapsePoorQualityNodes(self, bootstrap_threshold: float = 95) -> None:
        """
        Collapse all nodes with a bootstrap value smaller than threshold
        """
        self._scaleBootstrapValues()
        self._tree.collapse_all(
            lambda c: c.confidence is not None and c.confidence < bootstrap_threshold
            )

    def nameInternalNodes(self) -> None:
        """
        Give unique identifier to internal nodes
        including bootstrap value in identifier
        as: IN_n_b, where n is a unique identifying
        number and b the bootstrap value (or None
        if not present)
        """
        for n, clade in enumerate(self._tree.get_nonterminals()):
            if clade.name is None:
                if clade.confidence is not None:
                    clade.name = f'IN_{n}_{clade.confidence}'
                else:
                    clade.name = f'IN_{n}'
                clade.confidence = None

    def getTreeObject(self):
        return self._tree

    def exportTree(self, outfile: str, tree_format: str = 'newick') -> None:
        """
        Export tree object to file
        """
        Phylo.write(self._tree, outfile, tree_format)

    def getAllDescendantsOfTargetNode(self, target_name: str) -> list:
        """
        Get all leaf names from given target (internal) node name
        """
        target = next(self._tree.find_clades(target=target_name))
        return [n.name for n in target.get_terminals()]

    def getClosestCommonAncestor(self, target_names: list[str]) -> str:
        """
        Get name of closest common ancestor given list of leaf names
        """
        clade = self._tree.common_ancestor(target_names)
        return clade.name

    def extractClustersFromInternalNodes(self, filter_by_pattern: str = None) -> dict:
        """
        Extract all terminal nodes which are descendants of
        each internal node in the tree

        @params
        filter_by_pattern: filter clusters by leaf name pattern, such that
        all leaves within a cluster contain the same pattern
        """
        cluster_dict = {}
        for clade in self._tree.get_nonterminals():
            terminal_nodes = clade.get_terminals()
            cluster_dict[clade.name] = [n.name for n in terminal_nodes]
        if filter_by_pattern is None:
            return cluster_dict
        else:
            return self.filterClustersByCommonLeafNamePattern(cluster_dict, filter_by_pattern)
    
    @staticmethod
    def filterClustersByCommonLeafNamePattern(clusters: dict, pattern: str) -> dict:
        """
        Obtain (if any) tree clusters only containing leaves with leaf
        names containing specified pattern (sorry)
        """
        return {
            cluster_id: cluster
            for cluster_id, cluster in clusters.items()
            if all([pattern in leaf_name for leaf_name in cluster])
        }   


def setDefaultOutputPath(input_path: str, tag: str = None,
                         extension: str = None,
                         only_filename: bool = False,
                         only_dirname: bool = False) -> str:
    """
    Get default path to output file or directory
    """
    basename = os.path.basename(input_path)
    dirname = os.path.abspath(os.path.dirname(input_path))
    fname, ext = os.path.splitext(basename)
    if extension is None:
        extension = ext
    if tag is None:
        tag = ''
    default_file = f'{fname}{tag}{extension}'
    if only_filename:
        return default_file
    if only_dirname:
        return os.path.abspath(dirname)
    else:
        return os.path.abspath(os.path.join(dirname, default_file))

def filterFASTAbyIDs(input_fasta: str, record_ids: list[str],
                     output_fasta: str = None) -> None:
    """
    Filter records in fasta file matching provided IDs
    """
    if output_fasta is None:
       output_fasta = setDefaultOutputPath(input_fasta, '_fitered')
    record_ids = set(record_ids)
    fa = pyfastx.Fasta(input_fasta)
    with open(output_fasta, 'w') as fp:
        for record_id in record_ids:
            try:
                record_obj = fa[record_id]
                fp.write(record_obj.raw)
            except:
                pass
    os.remove(input_fasta + ".fxi")

def exportTreeClustersToFile(clusters: dict, outfile: str) -> None:
    """
    Write tsv file containing the definition of tree clusters
    """
    def getNodeCluster(node_name: str, clusters: dict):
        for cluster_name, cluster in clusters.items():
             if node_name in cluster:
                 return cluster_name

    with open(outfile, "w") as file:
        lines = ["id\tcluster\n"]
        node_names = [nname for cluster in clusters.values() for nname in cluster]
        for nname in node_names:
            cluster_name = getNodeCluster(nname, clusters)
            line = f"{nname}\t{cluster_name}\n"
            lines.append(line)
        file.writelines(lines)