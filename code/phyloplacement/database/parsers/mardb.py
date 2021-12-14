#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to process MARdb data
"""

import os
import re
import shutil

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
        db_entry = re.compile('_MMP\d{8}')
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
    """
    def __init__(self, current: str = None, partial: str = None):
        if current is None:
            self.current = os.path.abspath('data/taxonomy/CurrentComplete.tsv')
        else:
            self.current = os.path.abspath(current)
        if partial is None:
            self.partial = os.path.abspath('data/taxonomy/CurrentPartial.tsv')
        else:
            self.partial = os.path.abspath(partial)
        tax_cols = [
        'mmp_ID', 'taxon_lineage_names', 'full_scientific_name',
        'kingdom','phylum','class', 'order', 'family', 'genus', 'species'
        ]
        tax_complete = pd.read_csv(self.current, sep='\t')
        tax_complete = tax_complete[tax_cols]
        tax_partial = pd.read_csv(self.partial, sep='\t')
        tax_partial = tax_partial[tax_cols]
        tax_mmp = pd.concat([tax_complete, tax_partial]).drop_duplicates()
        self.tax_mmp = tax_mmp.set_index('mmp_ID')

    def lowestAvailableCommonTaxonomy(self, mmp_ids: list) -> dict:
        """
        Find lowest common GTDB taxonomy of selected MarDB entry ids
        """
        tax_levels = [
            'kingdom', 'phylum', 'class',
            'order', 'family', 'genus', 'species'
            ]
        lowest_tax = {}
        try:
            data = self.tax_mmp.loc[mmp_ids, :]
            for tax_cat in tax_levels:
                taxa = data[tax_cat].drop_duplicates().values
                if len(taxa) == 1:
                    lowest_tax[tax_cat] = taxa[0]
            return lowest_tax
        except Exception as e:
            raise ValueError(f'No taxonomy found. Exception: {e}')

    def assignTaxonomyToLabel(self, label: str,
                              root_level: str = 'kingdom',
                              full_label: bool = False) -> str:
        """
        Assign GTDB taxonomy to labels containing MardDB entry codes.
        Returns original label when no taxonomical info found.
        @Arguments:
            root_level: selects the highest taxon to be included in the label
            (kingdom, phylum, class, order, family)
            full_label: whether to output a simplified label with MMP id and
                        assigned taxonomy or the full label (defaults to False)
        """
        tax_levels = ['kingdom', 'phylum', 'class', 'order', 'family']
        selected_levels = tax_levels[tax_levels.index(root_level):]
        labelParser = MARdbLabelParser()
        mmp_id = labelParser.extractMMPid(label)

        try:
            tax_dict = self.lowestAvailableCommonTaxonomy([mmp_id])
            
            fil_tax_dict = {
                tax: value 
                for tax, value in tax_dict.items()
                if (
                    (value.lower() not in 'unclassified') and
                    (tax in selected_levels)
                    )
            }
            tax_labels = '_'.join(fil_tax_dict.values())
            if full_label:
                return f'{label}_{tax_labels}'
            else:
                return f'{mmp_id}_{tax_labels}'
        except Exception:
            tax_labels = 'No_taxomy_found'
            return f'{label}_{tax_labels}'



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