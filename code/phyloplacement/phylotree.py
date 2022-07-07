#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to perform phylogenetic tree reconstructions and
query sequence placements onto trees
"""

from __future__ import annotations
import os
import shutil
import re

from Bio import Phylo 

from phyloplacement.utils import setDefaultOutputPath, easyPatternMatching
import phyloplacement.wrappers as wrappers



class PhyloTree:

    def __init__(self, tree: str, tree_format: str = 'newick',
                 bootstrap_threshold: float = None, name_internal_nodes: bool = True):
        """
        Methods to help defining clusters in phylo trees
        @params:
        tree: str containing either the path to the tree file or directly the tree
        """
        self._tree = next(Phylo.parse(tree, tree_format))
        if bootstrap_threshold is not None:
            self.collapsePoorQualityNodes(bootstrap_threshold)
        if name_internal_nodes:
            self.nameInternalNodes()
        self._tree_path = tree
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
                clade.name = f'IN_{n}_{clade.confidence}'
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

    def getAllLeafNames(self) -> list[str]:
        """
        Get list of all leaves (terminal nodes names)
        """
        return [leaf.name for leaf in self._tree.get_terminals()]

    def extractClustersFromInternalNodes(self) -> dict:
        """
        Extract all terminal nodes which are descendants of
        each internal node in the tree
        """
        cluster_dict = {}
        for clade in self._tree.get_nonterminals():
            terminal_nodes = clade.get_terminals()
            cluster_dict[clade.name] = [n.name for n in terminal_nodes]
        return cluster_dict

    def computeTreeDiameter(self) -> float:
        """
        Find maximum (pairwise) distance between two tips
        (leaves) in the tree
        """
        root = self._tree.root
        max_distance = 0.0
        tips = self._tree.get_terminals()
        for tip in tips:
            self._tree.root_with_outgroup(tip)
            new_max = max(self._tree.depths().values())
            if new_max > max_distance:
                max_distance = new_max
        self._tree.root_with_outgroup(root)
        return max_distance

        

def inferTree(ref_aln: str,
              method: str = 'iqtree',
              substitution_model: str = 'modeltest',
              output_dir: str = None,
              additional_args: str = None) -> None:
    """
    Infer tree from reference msa. Best substitution model
    selected by default.
    """
    if method.lower() in 'iqtree':
        if 'modeltest' in substitution_model.lower():
            print('Selecting best subsitution model per modeltest-ng...')
            wrappers.runModelTest(
                input_algns=ref_aln,
                n_processes=None,
                output_dir=output_dir
            )
            best_model = getTreeModelFromModeltestLog(
                modeltest_log=os.path.join(output_dir, 'modeltest_result.log'),
                criterion='BIC'
            )
            modeltest_starting_tree = os.path.join(
                output_dir, 'modeltest_result.tree'
            )
        elif 'iqtest' in substitution_model.lower():
            best_model = 'TEST'
            modeltest_starting_tree = None
        else:
            best_model = substitution_model
            modeltest_starting_tree = None

        wrappers.runIqTree(
        input_algns=ref_aln,
        output_dir=output_dir,
        output_prefix='ref_database',
        keep_recovery_files=True,
        substitution_model=best_model,
        starting_tree=modeltest_starting_tree,
        additional_args=additional_args
        )
        shutil.move(
            os.path.join(output_dir, 'ref_database.contree'),
            os.path.join(output_dir, 'ref_database.newick')
            )

    elif method.lower() in 'fasttree':
        wrappers.runFastTree(
            input_algns=ref_aln,
            output_file=os.path.join(output_dir, 'ref_database.newick'),
            additional_args=additional_args
        )
    else:
        raise ValueError('Wrong method, enter iqtree or fasttree')

def sanityCheckForiTOL(label: str) -> str:
    """
    Reformat label to comply with iTOL requirements, remove:
    1. white spaces
    2. double underscores
    3. symbols outside english letters and numbers 
    """
    legal_chars = re.compile('[^a-zA-Z0-9]')
    itol_label = legal_chars.sub('_', label).replace('__', '_')
    return itol_label

def relabelTree(input_newick: str,
                label_dict: dict,
                output_file: str = None,
                iTOL=True) -> None: 
    """
    Relabel tree leaves with labels from 
    provided dictionary. If iTOL is set, then
    labels are checked for iTOL compatibility
    """
    if output_file is None:
        output_file = setDefaultOutputPath(
            input_newick, tag='_relabel'
        )
    if iTOL:
        sanity_check = sanityCheckForiTOL
    else:
        sanity_check = lambda x: x
    tree = next(Phylo.parse(input_newick, 'newick'))
    leaves = tree.get_terminals()
    for leaf in leaves:
        if leaf.name in label_dict.keys():
            leaf.name = sanity_check(label_dict[leaf.name])
    Phylo.write(tree, output_file, 'newick')

def newRelabelTree(input_newick: str,
                   label_dict: dict,
                   output_file: str = None,
                   iTOL=True) -> None: 
    """
    Relabel tree leaves with labels from 
    provided dictionary. If iTOL is set, then
    labels are checked for iTOL compatibility
    """
    if output_file is None:
        output_file = setDefaultOutputPath(
            input_newick, tag='_relabel'
        )
    if iTOL:
        sanity_check = sanityCheckForiTOL
    else:
        sanity_check = lambda x: x
    with open(input_newick, "r") as file:
        data = file.read()
    with open(output_file, "w") as file:
        for k, v in label_dict.items():
            data = data.replace(k, sanity_check(v))
        file.write(data)

def getIqTreeModelFromLogFile(iqtree_log: str) -> str:
    """
    Parse iqtree log file and return best fit model

    If model supplied, search model in Command: iqtree ... -m 'model' 
    If not, then -m TEST or -m MFP
    If one of those, continue to line:
    Best-fit model: 'model' chosen according to BIC
    """
    with open(iqtree_log, 'r') as log:
        text = log.read()
        subtext = re.search('(?<=Command: iqtree)(.*)(?=\\n)', text).group(1)
        model = re.search('(?<=-m )(.*)(?= -bb)', subtext).group(1)
        if model.lower() in ['mfp', 'test']:
            model = re.search('(?<=Best-fit model: )(.*)(?= chosen)', text).group(1)
    return model

def getTreeModelFromModeltestLog(modeltest_log: str, criterion: str = 'BIC') -> str:
    """
    Parse modeltest-ng log file and return best fit model
    according to selected criterion: BIC, AIC or AICc
    """
    with open(modeltest_log, 'r') as log:
        text = log.read()
        model = easyPatternMatching(
            easyPatternMatching(
                text, f'Best model according to {criterion}\n', '\nlnL'),
                left_pattern='Model:'
                ).strip()
        return model

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
