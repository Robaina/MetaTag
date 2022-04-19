#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to quantify and assign labels to placed sequences
"""
from __future__ import annotations
import os
import re
import json
from io import StringIO

import pandas as pd
from Bio import Phylo 

import phyloplacement.wrappers as wrappers
from phyloplacement.utils import (setDefaultOutputPath,
                                  TemporaryFilePath)
from phyloplacement.taxonomy import TaxonomyAssigner, TaxonomyCounter, Taxopath
from phyloplacement.alignment import alignShortReadsToReferenceMSA
from phyloplacement.phylotree import getIqTreeModelFromLogFile
from phyloplacement.database.manipulation import getFastaRecordIDs, splitReferenceFromQueryAlignments


class JplaceParser():
    """
    Methods to parse jplace files, as specified in 
    https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009
    """
    def __init__(self, path_to_jplace: str) -> None:
        with open(path_to_jplace, 'r') as JSON:
            jplace = json.load(JSON)
        self.jplace = jplace
    
    @property
    def meta(self):
        """
        Print metadata
        """
        return self.jplace['metadata']

    @property
    def fields(self): 
        """
        Print data fields
        """
        return self.jplace['fields']

    def tree(self, newick=False):
        """
        Return tree in original or newick format
        Original format contains branch labels
        in curly brackets. Newick format removes
        these labels.
        """
        if newick:
            return self.newickfyTree(self.jplace['tree'])
        else:
            return self.jplace['tree']
    
    @property
    def placements(self):
        """
        Return placement objects
        """
        return self.jplace['placements']
    
    @staticmethod
    def newickfyTree(tree_str: str) -> str:
        """
        Remove branch IDs from jplace tree string
        """
        subs_tree = re.sub("\{(\d+)\}", '', tree_str)
        return next(Phylo.parse(StringIO(subs_tree), 'newick'))

    def getReferenceSequences(self) -> list:
        """
        Get list of reference sequences in the placement tree
        """
        tree = self.tree(newick=True)
        return [c.name for c in tree.get_terminals()]
    
    def buildBranchDict(self) -> dict:
        """
        Build dictionary with edge/branch numbers as keys and 
        reference tree leaves as values
        """
        def get_id(s):
            return int(re.search("\{(\d+)\}", s).group(1))
        
        original_tree = self.jplace['tree']
        tree = self.newickfyTree(original_tree)
        leaves = tree.get_terminals()

        branches = {
            get_id(original_tree[original_tree.find(leaf.name):]): leaf.name
            for leaf in leaves
        }
        return branches

    def extractPlacementFields(self, pfielddata: list) -> dict:
        """
        Get dict with placement field values from list of values
        """
        fields = self.jplace['fields']
        return {field: pfielddata[i] for i, field in enumerate(fields)}

    def selectBestPlacement(self, placement_object: dict) -> dict:
        """
        Select placement with lowest likelihood
        """
        pdata = [
            self.extractPlacementFields(pfielddata)
            for pfielddata in placement_object['p']
            ]
        lowest_like_placement = sorted(pdata, key=lambda x: x['likelihood'])[0]
        return {'p': lowest_like_placement, 'n': placement_object['n']}

    def selectBestPlacements(self):
        """
        Select placement with lowest likelihood for 
        all placement objects in placements
        """
        best_placements = [
            self.selectBestPlacement(placement)
            for placement in self.jplace['placements']
        ]
        return best_placements

    def assignClusterCommonLabelsToQueries(self) -> dict:
        """
        Assign cluster function and lowest common taxonomy
        to query sequences placed within cluster
        """
        pass

    def labelPlacedQueries(self):
        """
        Assign reference label to best query placements
        """
        ref_dict = self.buildBranchDict()
        best_placements = self.selectBestPlacements()
        labeled_placements = {
            place_object['n'][0]: ref_dict[place_object['p']['edge_num']]
            for place_object in best_placements
            if place_object['p']['edge_num'] in ref_dict.keys()
        }
        return labeled_placements


class TaxAssignParser():
    """
    Parse function and taxonomy placement assignments table
    """
    def __init__(self, tax_assign_path: str):
                 self._tax_assign = pd.read_csv(tax_assign_path, sep='\t')
                 self._tax_assign = self._tax_assign[self._tax_assign.cluster_id != 'DISTANT']
                 self._tax_assign.cluster_score = self._tax_assign.cluster_score.apply(lambda x: float(x))
                 self._tax_assign.cluster_taxopath = self._tax_assign.cluster_taxopath.apply(lambda x: "Unspecified" if pd.isna(x) else x)
                 self._tax_assign.taxopath = self._tax_assign.taxopath.apply(lambda x: "Unspecified" if pd.isna(x) else x)
    @property
    def taxlevels(self):
        return Taxopath().taxlevels

    def countHits(self, 
                  cluster_ids: list[str] = None, score_threshold: float = None,
                  taxopath_type: str = 'taxopath') -> TaxonomyCounter:
        """
        Count hits within given cluster ids and at specificied taxon level
        @Params:
        normalize=True, then results are reported as fractions of total
        cluster_ids: IDs of tree clusters to be included in the counting of placements
        score_threshold: global placement score threshold to filter
                         low-quality placements out
        taxopath_type: 'taxopath' to use gappa-assign taxopath or 'cluster_taxopath'
                        to use lowest common taxopath of the reference tree cluster
        """
        if cluster_ids is None:
            cluster_ids = self._tax_assign.cluster_id.unique()
        if score_threshold is None:
            score_threshold = 0.0

        taxopath_hits = self._tax_assign[
            (
                (self._tax_assign.cluster_id.isin(cluster_ids)) &
                (self._tax_assign.cluster_score >= score_threshold)
                )
            ][taxopath_type].values
        if len(taxopath_hits) == 0:
            raise ValueError('No placement hits returned for the provided filter parameters')

        taxlevel_counter = TaxonomyCounter(taxopath_list=taxopath_hits)
        return taxlevel_counter


def placeReadsOntoTree(input_tree: str, 
                       tree_model: str,
                       ref_aln: str,
                       query_seqs: str,
                       aln_method: str = 'papara',
                       output_dir: str = None) -> None:
    """
    Performs short read placement onto phylogenetic tree
    tree_model: str, either the model name or path to log output by iqtree
    workflow example: https://github.com/Pbdas/epa-ng/wiki/Full-Stack-Example
    """
    if output_dir is None:
        output_dir = setDefaultOutputPath(query_seqs, only_dirname=True)
    else:
        output_dir = os.path.abspath(output_dir)
    
    if os.path.isfile(tree_model):
        tree_model = getIqTreeModelFromLogFile(tree_model)
        print(f'Running EPA-ng with inferred substitution model: {tree_model}')
    
    ref_ids = getFastaRecordIDs(ref_aln)
    ref_query_msa = os.path.join(
        output_dir, setDefaultOutputPath(query_seqs, extension='.faln',
                                         only_filename=True)
        )
    aln_ref_frac = os.path.splitext(ref_query_msa)[0] + '_ref_fraction.faln'
    aln_query_frac = os.path.splitext(ref_query_msa)[0] + '_query_fraction.faln'

    alignShortReadsToReferenceMSA(
        ref_msa=ref_aln,
        query_seqs=query_seqs,
        method=aln_method,
        tree_nwk=input_tree,
        output_dir=output_dir
    )
    
    splitReferenceFromQueryAlignments(
        ref_query_msa=ref_query_msa,
        ref_ids=ref_ids,
        out_dir=output_dir
    )

    wrappers.runEPAng(
        input_tree=input_tree,
        input_aln_ref=aln_ref_frac,
        input_aln_query=aln_query_frac,
        model=tree_model,
        output_dir=output_dir,
        n_threads=None,
        additional_args=None)

def parseTreeClusters(clusters_tsv: str, cluster_as_key: bool = True, sep='\t') -> dict:
    """
    Parse clusters text file into dictionary
    @param
    clusters_as_key: if True then dict keys are cluster IDs and values
    are lists of reference IDs. If False, dict keys are reference IDs and 
    values the cluster to which they belong.
    """
    df = pd.read_csv(clusters_tsv, sep='\t')
    if cluster_as_key:
        cluster_ids = df.cluster.unique()
        return {
            cluster_id: df[df.cluster == cluster_id].id.values.tolist()
            for cluster_id in cluster_ids
        }
    else:
        return dict(df.values)

def parseTreeClusterQualityScores(cluster_scores_tsv: str, sep='\t') -> dict:
    """
    Parse cluster quality scores file into dictionary
    """
    df = pd.read_csv(cluster_scores_tsv, sep='\t')
    return dict(df.values)

def addClustersToTaxTable(in_taxtable: str, clusters: dict,
                          out_taxtable: str = None) -> None:
    """
    Add tree cluster info at the beginning of each taxopath
    according to clusters defined in dictionary 'clusters'
    """
    if out_taxtable is None:
        out_taxtable = setDefaultOutputPath(in_taxtable, tag='_clustered')
    taxtable = pd.read_csv(in_taxtable, sep='\t', header=None)
    for i, row in taxtable.iterrows():
        row[1] = clusters[row[0]] + ';' + row[1]
    taxtable.to_csv(out_taxtable, sep='\t', index=False, header=None)

def parseGappaAssignTable(input_table: str, has_cluster_id: bool = True,
                          cluster_scores: dict = None,
                          clusters_taxopath: dict = None,
                          output_file: str = None) -> None:
    """
    Parse gappa assign per query taxonomy assignment result tsv
    has_cluster_id: set to True if results table includes reference
                    cluster info in the first element of taxopath
    cluster_scores: dictionary with values set to cluster quality
                    scores. It is only used if has_cluster_id = True.
    """
    if output_file is None:
        output_file = setDefaultOutputPath(input_table, tag='_parsed')
    if (has_cluster_id) and (cluster_scores is not None):
        cluster_str = 'cluster_id' + '\t' + 'cluster_score' + '\t' + 'cluster_taxopath' + '\t'
    elif (has_cluster_id) and (cluster_scores is None):
        cluster_str = 'cluster_id' + '\t' + 'cluster_taxopath' + '\t'
    elif not has_cluster_id:
         cluster_str = ""
    table = pd.read_csv(input_table, sep='\t')
    with open(output_file, 'w') as file:
        lines = []
        header = 'query_id' + '\t' + cluster_str + 'taxopath' + '\n'
        lines.append(header)
        for i, row in table.iterrows():
            if not ';' in row.taxopath:
                row.taxopath = row.taxopath + ';Unspecified'
            if has_cluster_id:
                elems = row.taxopath.split(';')
                cluster_id = elems[0]
                if cluster_id in clusters_taxopath.keys():
                    cluster_taxopath = clusters_taxopath[cluster_id]
                else:
                    cluster_taxopath = ''
                if not cluster_taxopath:
                    cluster_taxopath = 'Unspecified'
                taxopath = cluster_taxopath + '\t' + ';'.join(elems[1:])
            else:
                cluster_id, taxopath = '', row.taxopath
            if cluster_scores is not None:
                if cluster_id in cluster_scores.keys():
                    cluster_score = str(cluster_scores[cluster_id])
                else:
                    cluster_score = str(None)
                line = row['name'] + '\t' + cluster_id + '\t' + cluster_score + '\t' + taxopath + '\n'
            else:
                line = row['name'] + '\t' + cluster_id + '\t' + taxopath + '\n'
    
            lines.append(line)
        file.writelines(lines)

def addQueryLabelsToAssignTable(input_table: str,
                                query_labels: dict,
                                output_table: str = None) -> None:
    """
    Add new column containing actual query labels to query taxonomy/cluster assignment
    table
    """
    def relabel_query(query_id: str, query_labels: dict) -> str:
        try:
            return query_labels[query_id]
        except:
            return query_id

    if output_table is None:
        output_table = setDefaultOutputPath(input_table, tag="_relabel")
    df = pd.read_csv(input_table, sep="\t")
    df.insert(1, "query_name", df.query_id.apply(lambda x: relabel_query(x, query_labels)))
    df = df.set_index("query_id")
    df.to_csv(output_table, sep="\t")

def assignLabelsToPlacements(jplace: str, id_dict: dict,
                             output_dir: str = None,
                             output_prefix: str = None,
                             only_best_hit: bool = True,
                             ref_clusters_file: str = None,
                             ref_cluster_scores_file: str = None,
                             gappa_additional_args: str = None) -> None:
    """
    Assign taxonomy and/or tree cluster IDs to placed query sequences based on
    taxonomy assigned to tree reference sequences using
    gappa examine assign.
    @parameter
    only_best_hit: only report taxonomy with largest LWR
                   per query
    ref_clusters: dict (optionally) add tree cluster to the
                  beginning of the taxopath so query sequences
                  can also be classified according to tree
                  cluster (e.g., assigned function)
    """
    if output_dir is None:
        output_dir = setDefaultOutputPath(jplace, only_dirname=True)
    if output_prefix is None:
        output_prefix = setDefaultOutputPath(jplace, only_filename=True)

    if ref_clusters_file is not None:
        has_cluster_id = True
        ref_clusters = parseTreeClusters(ref_clusters_file, cluster_as_key=False, sep='\t')
        ref_clusters_as_keys = parseTreeClusters(ref_clusters_file, cluster_as_key=True, sep='\t')
    else:
        has_cluster_id = False

    if ref_cluster_scores_file is not None:
        ref_cluster_scores = parseTreeClusterQualityScores(ref_cluster_scores_file)
    else:
        ref_cluster_scores = None

    gappa_assign_out = os.path.join(output_dir, output_prefix + 'per_query.tsv')
    query_taxo_out = os.path.join(output_dir, output_prefix + 'assignments.tsv')
    
    # Remove references that are not in placement tree
    ref_in_jplace_tree = JplaceParser(jplace).getReferenceSequences()
    id_dict = {k: v for k,v in id_dict.items() if k in ref_in_jplace_tree}

    with TemporaryFilePath() as temptax:
        taxonomy = TaxonomyAssigner(
            taxo_file='./data/taxonomy/merged_taxonomy.tsv'
        )
        taxonomy.buildGappaTaxonomyTable(id_dict, output_file=temptax)
        if has_cluster_id:
            addClustersToTaxTable(
                in_taxtable=temptax,
                clusters=ref_clusters,
                out_taxtable=temptax
            )
            clusters_taxopath = taxonomy.assignLowestCommonTaxonomyToClusters(
                clusters=ref_clusters_as_keys,
                label_dict=id_dict
            )
        else:
            clusters_taxopath = None
        wrappers.runGappaAssign(
            jplace=jplace,
            taxonomy_file=temptax,
            output_dir=output_dir,
            output_prefix=output_prefix,
            only_best_hit=only_best_hit,
            additional_args=gappa_additional_args)
    
    with TemporaryFilePath() as tempout:
        parseGappaAssignTable(
            input_table=gappa_assign_out,
            has_cluster_id=has_cluster_id,
            cluster_scores=ref_cluster_scores,
            output_file=tempout,
            clusters_taxopath=clusters_taxopath
        )
        addQueryLabelsToAssignTable(
            input_table=tempout,
            query_labels=id_dict,
            output_table=query_taxo_out
        )

    
