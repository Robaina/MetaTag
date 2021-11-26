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


class MMPtaxonomyAssigner():
    """
    Methods to assign taxonomy to mardb reference sequeces
    """
    def __init__(self, current: str = None, partial: str = None):
        self.current = (os.path.abspath(current)
                        or os.path.abspath('/data/taxonomy/CurrentComplete.tsv'))
        self.partial = (os.path.abspath(partial)
                        or os.path.abspath('/data/taxonomy/CurrentPartial.tsv'))
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
        tax_list = [
            'kingdom', 'phylum', 'class',
            'order', 'family', 'genus', 'species'
            ]
        lowest_tax = {}
        data = self.tax_mmp.loc[mmp_ids, :]
        for tax_cat in tax_list:
            taxa = data[tax_cat].drop_duplicates().values
            if len(taxa) == 1:
                lowest_tax[tax_cat] = taxa[0]
        return lowest_tax


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