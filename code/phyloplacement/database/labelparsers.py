#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to parse sequence labels from different databases
"""

import re
import warnings

import pandas as pd



class LabelParser():
    def __init__(self) -> None:
        """
        Parse labels to extract genome ID and metadata
        """
        return None
    
    @staticmethod
    def extractGenomeID(label: str) -> str:
        """
        Extract genome ID from sequence label
        """
        mmp_id = MARdbLabelParser.extractMMPid(label)
        taxid = UniprotLabelParser.extractNcbiTaxId(label)
        if mmp_id and taxid:
            warnings.warn("Label contains conflicting genome IDs")
            return ''
        elif mmp_id:
            genome_id = mmp_id
        elif taxid:
            genome_id = taxid
        else:
            genome_id_split = label.split('__')
            if len(genome_id_split) > 1:
                genome_id = genome_id_split[0]
            else:
                genome_id = ''
        return genome_id

    @staticmethod
    def parseMetaInfo(label: str) -> dict:
        """
        Parse meta information from label in canonical format
        (i.e, contig, gene position, locus, strand)
        """
        try:
            meta = label.split('__')[1]
            strand = meta.split('_')[-1]
            locus_pos = tuple([int(pos) for pos in meta.split('_')[-3:-1]])
            gene_pos = int(meta.split('_')[-4])
            contig = '_'.join(meta.split('_')[:-4])
        except:
            contig = None
            gene_pos = None
            locus_pos = None 
            strand = None

        return {
            "contig": contig,
            "gene_pos": gene_pos,
            "locus_pos": locus_pos,
            "strand": strand
        }
  

class UniprotLabelParser():

    def __init__(self) -> None:
        """
        Parse UniProt entry label to extract coded info
        """
        return None

    @staticmethod
    def extractNcbiTaxId(label: str) -> str:
        """
        Extract NCBI taxon id from reference label
        """
        db_entry = re.compile('(OX=)(\d+)')   
        try:
            return f'taxid_{re.search(db_entry, label).group(2)}'
        except Exception:
            return ''


class MARdbLabelParser():

    def __init__(self) -> None:
        """
        Parse MARdb entry label to extract coded info
        """
        return None

    @staticmethod
    def extractMMPid(label: str) -> str:
        """
        Extract mardb mmp id from reference label
        """
        db_entry = re.compile('_MMP\d+')
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