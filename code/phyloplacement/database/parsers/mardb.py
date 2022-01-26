#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to process MARdb data
"""

import os
import re
import shutil
import warnings

import pyfastx
import pandas as pd

from phyloplacement.utils import (setDefaultOutputPath,
                                  terminalExecute,
                                  createTemporaryFilePath) 


class MARdbLabelParser():
    """
    Parse MARdb entry label to extract coded info
    """
    def __init__(self) -> None:
        pass

    @staticmethod
    def extractMMPid(label: str) -> str:
        """
        Extract mardb mmp id from reference label
        """
        db_entry = re.compile('_MMP\d+')#{7}|_MMP\d{8}|_MMP\d{9}')
        try:
            return re.search(db_entry, label).group(0).strip('_')
        except Exception:
            return ''
    
    def parse(self, label: str) -> dict:
        """
        Parse MarDB sequence labels to obtain contig and locus info
        """
        parsed_dict = {
            'full': label,
            'species': '',
            'mmp_id': '',
            'contig': '',
            'gene_pos': None,
            'locus_pos': None,
            'strand': ''
            } 
        try: 
            entry = label.split('__')[0]
            mmp_id = self.extractMMPid(label)
            species = entry.strip(mmp_id).strip('_')
            meta = label.split('__')[1]
            strand = meta.split('_')[-1]
            locus_pos = tuple([int(pos) for pos in meta.split('_')[-3:-1]])
            gene_pos = int(meta.split('_')[-4])
            contig = '_'.join(meta.split('_')[:-4])
            parsed_dict['species'] = species
            parsed_dict['mmp_id'] = mmp_id
            parsed_dict['contig'] = contig
            parsed_dict['gene_pos'] = gene_pos 
            parsed_dict['locus_pos'] = locus_pos
            parsed_dict['strand'] = strand
        except Exception:
            pass
        return parsed_dict

    def parse_from_list(self, labels=list) -> pd.DataFrame: 
        """
        Parse labels in list of labels and return DataFrame
        """
        return pd.DataFrame(
            [self.parse(label) for label in labels]
        )


class MMPtaxonomyAssigner():
    """
    Methods to assign taxonomy to mardb reference sequeces
    @Params:
    complete: path to complete gtdb database
    partial: path to partial gtdb database
    """
    def __init__(self, complete: str = None, partial: str = None):
        if complete is None:
            self.complete = os.path.abspath('data/taxonomy/CurrentComplete.tsv')
        else:
            self.complete = os.path.abspath(complete)
        if partial is None:
            self.partial = os.path.abspath('data/taxonomy/CurrentPartial.tsv')
        else:
            self.partial = os.path.abspath(partial)
        tax_cols = [
        'mmp_ID', 'taxon_lineage_names', 'full_scientific_name',
        'kingdom','phylum','class', 'order', 'family', 'genus', 'species'
        ]
        tax_complete = pd.read_csv(self.complete, sep='\t')
        tax_complete = tax_complete[tax_cols]
        tax_partial = pd.read_csv(self.partial, sep='\t')
        tax_partial = tax_partial[tax_cols]
        tax_mmp = pd.concat([tax_complete, tax_partial]).drop_duplicates()
        self.tax_mmp = tax_mmp.set_index('mmp_ID')

    def lowestAvailableCommonTaxonomy(self, mmp_ids: list,
                                      only_classified_taxa: bool = True) -> dict:
        """
        Find lowest common GTDB taxonomy of selected MarDB entry ids
        @Params:
        only_classified_taxa: report only taxa with assigned taxonomy
        """
        mmps = [mmp for mmp in mmp_ids if mmp in self.tax_mmp.index]
        if len(mmps) < len(mmp_ids):
            warnings.warn('Not all labels had taxopath. Computing common taxonomy from a subset')
        tax_levels = [
            'kingdom', 'phylum', 'class',
            'order', 'family', 'genus', 'species'
            ]
        lowest_tax, taxdict = {}, {}
        data = self.tax_mmp.loc[mmps, :]
        for tax_cat in tax_levels:
            taxa = data[tax_cat].drop_duplicates().values
            if len(taxa) == 1:
                lowest_tax[tax_cat] = taxa[0]
            else:
                lowest_tax[tax_cat] = 'Unclassified'
        if only_classified_taxa:
            for tax, value in lowest_tax.items():
                if value == 'Unclassified':
                    break
                taxdict[tax] = value
        else:
            taxdict = lowest_tax
        if not taxdict:
            warnings.warn('No common taxonomy found for given MMPs')
        return taxdict

    def assignTaxonomyToLabel(self, label: str,
                              root_level: str = 'kingdom',
                              full_label: bool = False,
                              only_taxonomy: bool = False,
                              separator: str = '_') -> str:
        """
        Assign GTDB taxonomy to labels containing MardDB entry codes.
        Returns original label when no taxonomical info found.
        @Arguments:
            root_level: selects the highest taxon to be included in the label
            (kingdom, phylum, class, order, family)
            full_label: whether to output a simplified label with MMP id and
                        assigned taxonomy or the full label (defaults to False)
            only_taxonomy: whether to ouput only taxopath (defaults to False)
        """
        tax_levels = ['kingdom', 'phylum', 'class', 'order', 'family']
        selected_levels = tax_levels[tax_levels.index(root_level):]
        labelParser = MARdbLabelParser()
        mmp_id = labelParser.extractMMPid(label)
        if mmp_id:
            tax_dict = self.lowestAvailableCommonTaxonomy([mmp_id])
            fil_tax_dict = {}
            for tax, value in tax_dict.items():
                if tax in selected_levels:
                    if value.lower() in 'unclassified':
                        break
                    fil_tax_dict[tax] = value
            if fil_tax_dict:
                tax_labels = separator.join(fil_tax_dict.values())
                if full_label:
                    return tax_labels if only_taxonomy else f'{label}_{tax_labels}'
                else:
                    return tax_labels if only_taxonomy else f'{mmp_id}_{tax_labels}'
            else:
                tax_labels = 'No_taxonomy_found'
                return tax_labels if only_taxonomy else f'{label}_{tax_labels}'
        else:
            tax_labels = 'No_taxonomy_found'
            return tax_labels if only_taxonomy else f'{label}_{tax_labels}'

    def assignLowestCommonTaxonomyToCluster(self, clusters: dict, label_dict: dict) -> dict:
        """
        Find lowest possible common taxonomy to reference labels in clusters
        """
        parser = MARdbLabelParser()
        clusters_taxopath = {}
        for cluster_id, cluster in clusters.items():
            cluster_labels = [label_dict[ref_id] for ref_id in cluster]
            cluster_mmp_ids = [parser.extractMMPid(label) for label in cluster_labels]
            taxo_dict = self.lowestAvailableCommonTaxonomy(cluster_mmp_ids)
            taxopath = ';'.join(taxo_dict.values())
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
                taxon_str = self.assignTaxonomyToLabel(
                    label=label, full_label=False,
                    only_taxonomy=True, separator=';'
                    )
                taxon_str = 'Undefined' if ('No_taxonomy_found' in taxon_str) else taxon_str
                if taxon_str != 'Undefined':
                    lines.append(f'{ref_id}\t{taxon_str}\n')
            outfile.writelines(lines)



# *** Tagged as possibly trash code ***

def getMarDBentryCode(label: str) -> str:
    db_entry = re.compile('\[mmp_id=(.*)\] ')
    return re.search(db_entry, label).group(1)

def filterMarDBrecordsbyEntryCodes(input_fasta: str, entry_codes: set,
                                   output_fasta: str = None) -> None:
    """
    Filter records in mardb fasta file matching provided entry codes
    """
    if output_fasta is None:
        output_fasta = setDefaultOutputPath(input_fasta, '_fitered')
    
    fasta = pyfastx.Fasta(input_fasta, build_index=False, full_name=True)
    with open(output_fasta, 'w') as outfile:
        for record_name, record_seq in fasta:
            entry_code = getMarDBentryCode(record_name)
            if entry_code in entry_codes:
                outfile.write(f'>{record_name}\n{record_seq}\n')

def getMARdbGenomeByEntryCode(entry_code: str, input_fasta: str,
                              output_fasta: str = None,
                              clean_seqs: bool = True) -> None:
    """
    Get full or partial genomes with given MARdb entry code.
    If clean = True, remove characters which are not letters
    """
    if output_fasta is None:
        output_fasta = setDefaultOutputPath(input_fasta,
                                            tag=f'_genome_{entry_code}',
                                            extension='.fa')

    def cleanOutputFasta(output_fasta: str) -> None:
        """
        Check if illegal symbols in sequences,
        then remove and tag file as cleaned
        """
        not_capital_letters = re.compile('[^A-Z]')
        fname, ext = os.path.splitext(output_fasta)
        was_cleaned = False
        cleaned_fasta = f'{fname}_cleaned{ext}'
        temp_file_path = createTemporaryFilePath()
        with open(output_fasta, 'r') as fasta, open(temp_file_path, 'a+') as tfasta:
            for line in fasta.readlines():
                if ('>' not in line) and (not_capital_letters.search(line)):
                    line = not_capital_letters.sub('', line)
                    was_cleaned = True
                tfasta.write(line)
        if was_cleaned:
            shutil.move(fasta.name, cleaned_fasta)
            shutil.move(temp_file_path, cleaned_fasta)
        else:
            os.remove(temp_file_path)

    cmd_str = (
        f'grep -A1 {entry_code} {input_fasta} > {output_fasta}'
    )
    terminalExecute(cmd_str, suppress_shell_output=False)
    if clean_seqs:
        cleanOutputFasta(output_fasta)

def relabelMarDB(label_dict: dict) -> dict:
    """
    Convert mardb long labels into short labels 
    displaying mardb id and species (if present)
    """
    db_code_pattern = re.compile('\[mmp_(.*)\]')
    species_pattern = re.compile('\[(.*?)\]')

    def editMarDBlabel(label: str) -> str:
        try:
            species = re.search(
                species_pattern,
                re.sub(db_code_pattern, '', label)
                ).group(1)
        except:
            species = 'Undetermined'
        mar_id = label.split(' ')[0]
        return f'{mar_id}_{species}'

    return {
        k: editMarDBlabel(v)
        for k, v in label_dict.items()
    }