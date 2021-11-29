#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Obtain UreC (alpha subunit) peptide sequences from list of GenBank entries
"""

import os 
import pandas as pd
from Bio import SeqIO
from phyloplacement.utils import terminalExecute, fullPathListDir


gbk_dir = 'genes/ureA/data/koper_2004_gbks/'
output_fasta = 'genes/ureA/data/koper_2004_seqs.fasta'


def downloadGBKfromNCBI(entry_ids: list, output_dir: str = None) -> None: 
    """
    Download genbank files from NCBI from given list of entry IDs
    """
    for n, entry_id in enumerate(entry_ids): 
        print(f'Downloading entry: {entry_id} ({n + 1} / {len(entry_ids)})')
        outfasta = os.path.join(output_dir, f'{entry_id}.gbk')
        cmd_str = f'ncbi-acc-download -o {outfasta} {entry_id}'
        terminalExecute(cmd_str)

def getProteinSequenceFromGBK(gbk: str, keywords: list) -> dict:
    """
    Extract cds record matching keywords from gbk file
    """
    def contains_keywords(text: str, keywords: list) -> bool: 
        return all([key in text for key in keywords])

    gbk = list(SeqIO.parse(gbk, 'genbank'))[0]
    cds_records = [f.qualifiers for f in gbk.features[1:] if 'CDS' in f.type]
    return [cds for cds in cds_records if ('product' in cds.keys() and contains_keywords(cds['product'][0], keywords))]

def writeFASTAfromCDSqualifiers(records: dict, output_fasta: str = None) -> None: 
    """
    Write FASTA file from dict of gbk record qualifiers (Ordered dict)
    """
    with open(output_fasta, 'w') as file:
        for record_id, record in records.items():
            if record:
                record = record[0]
                ref_id = f'{record_id}_{record["protein_id"][0]}_{"_".join(record["product"][0].split())}'
                file.write(f'>{ref_id}\n')
                file.write(f'{record["translation"][0]}\n')


if __name__ == '__main__':
    
    # entry_list = pd.read_csv(
    #     'genes/ureA/data/Koper_2004_seq.txt', sep='\t', header=None,
    #     names=['species', 'entry']
    #     )
    # downloadGBKfromNCBI(entry_list.entry.values, output_dir=gbk_dir)

    records_dict = {
        gbk_file.split('.gbk')[0]: getProteinSequenceFromGBK(
            os.path.join(gbk_dir, gbk_file), keywords=['urease', 'alpha']
            ) for gbk_file in os.listdir(gbk_dir)
    }

    writeFASTAfromCDSqualifiers(records_dict, output_fasta=output_fasta)