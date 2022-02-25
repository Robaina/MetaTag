#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to assign taxonomy to reference and query (placed) sequences
"""

from __future__ import annotations
import os
from typing import List

import pandas as pd

from phyloplacement.database.parsers.mardb import MARdbLabelParser



class Taxopath():
    """
    Object to contain taxopath
    """
    def __init__(self, taxopath_str: str = None, delimiter: str = ";"):
        self._taxopath = taxopath_str
        self._delimiter = delimiter
        self._tax_levels = [
            'domain', 'phylum', 'class',
            'order', 'family', 'genus', 'species'
            ]
        self.taxodict = self._dictFromTaxopath()

    def _dictFromTaxopath(self):
        if self._taxopath is None:
            taxolist = []
        else:
            taxolist = [elem.strip() for elem in self._taxopath.split(self._delimiter)]
        taxolist.extend([None for _ in range(len(self._tax_levels) - len(taxolist))])
        return {taxlevel: taxon for taxlevel, taxon in zip(self._tax_levels, taxolist)}

    @classmethod
    def from_dict(cls, taxodict: dict, delimiter: str = ";") -> Taxopath:
        """
        Instantiate Taxopath object from dict
        """
        taxa = []
        for taxon in taxodict.values():
            if taxon is None:
                break
            else:
                taxa.append(taxon)
        taxostr = delimiter.join(taxa)
        return cls(taxopath_str=taxostr, delimiter=delimiter)

    @classmethod
    def getLowestCommonTaxopath(cls, taxopaths: list[str]) -> Taxopath:
        """
        compute lowest common taxopath (ancestor) of a list 
        of taxopaths
        """
        taxopath_dicts = [cls(taxostr).taxodict for taxostr in taxopaths]
        common_taxodict = cls().taxodict
        for taxlevel in cls().taxlevels:
            taxa = set([taxdict[taxlevel] for taxdict in taxopath_dicts])
            if len(taxa) > 1:
                break
            else:
                common_taxodict[taxlevel] = list(taxa)[0]
        return cls.from_dict(common_taxodict)

    @property
    def taxostring(self):
        return self._taxopath

    @property
    def taxlevels(self):
        return self._tax_levels



class TaxonomyAssigner():
    """
    Methods to assign taxonomy to reference sequences
    """
    def __init__(self, taxo_file: str):
        self._taxodata = pd.read_csv(
            os.path.abspath(taxo_file), sep='\t'
            ).drop_duplicates(subset='genome').set_index('genome')
    
    @staticmethod
    def lowestCommonTaxonomy(taxopaths: List[str]) -> str:
        """
        Find lowest common taxonomy among set of taxopaths
        """
        data = pd.DataFrame([t.split(';')for t in taxopaths])
        taxlevels = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
        lowest_tax = []
        for n, tax_level in enumerate(taxlevels):
                taxa = data.iloc[:, n].drop_duplicates().values
                if len(taxa) == 1:
                    lowest_tax.append(taxa[0])
                else:
                    break
        return ';'.join(lowest_tax)

    def _extractGenomeIDfromLabel(self, label: str) -> str:
        labelParser = MARdbLabelParser()
        mmp_id = labelParser.extractMMPid(label)
        if mmp_id:
            genome_id = mmp_id
        else:
            genome_id = label.split('__')[0]
        return genome_id

    def assignTaxonomyToLabel(self, label: str) -> str:
        """
        Assign GTDB taxonomy to label based on genome ID
        """
        genome_id = self._extractGenomeIDfromLabel(label)
        if genome_id in self._taxodata.index:
            return self._taxodata.loc[genome_id].item()
        else:
            return 'No_taxonomy_found'

    def assignLowestCommonTaxonomyToLabels(self, labels: List[str]) -> str:
        """
        Assing taxonomy to set of labels and find lowest common taxonomy
        among them
        """
        taxopaths = [
            taxopath for taxopath in map(self.assignTaxonomyToLabel, labels)
            if taxopath != 'No_taxonomy_found'
        ]
        if taxopaths:
            return self.lowestCommonTaxonomy(taxopaths)
        else:
            return 'Unspecified'

    def assignLowestCommonTaxonomyToClusters(self, clusters: dict, label_dict: dict = None) -> dict:
        """
        Find lowest possible common taxonomy to reference labels in clusters
        If reference labels do not contain genome IDs, a dictionary, label_dict,
        of reference labels and genome ids (or labels with genome ids) must be passed
        """
        clusters_taxopath = {}
        for cluster_id, cluster in clusters.items():
            if label_dict is not None:
                cluster_labels = [label_dict[ref_id] for ref_id in cluster]
            else:
                cluster_labels = cluster
            taxopath = self.assignLowestCommonTaxonomyToLabels(cluster_labels)
            clusters_taxopath[cluster_id] = taxopath
        return clusters_taxopath

    def buildGappaTaxonomyTable(self, ref_id_dict: dict, output_file: str = None) -> None:
        """
        Build gappa-compatible taxonomy file as specified here:
        https://github.com/lczech/gappa/wiki/Subcommand:-assign
        """
        if output_file is None:
            output_file = os.path.join(os.getcwd(), 'gappa_taxonomy.tsv')

        with open(output_file, 'w') as outfile:
            lines = []
            for ref_id, label in ref_id_dict.items():
                taxon_str = self.assignTaxonomyToLabel(label)
                taxon_str = 'Unspecified' if ('No_taxonomy_found' in taxon_str) else taxon_str
                if taxon_str != 'Unspecified':
                    lines.append(f'{ref_id}\t{taxon_str}\n')
            outfile.writelines(lines)