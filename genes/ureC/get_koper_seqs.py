#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Obtain UreC (alpha subunit) peptide sequences from list of GenBank entries
"""

import pandas as pd
from phyloplacement.database.parsers.ncbi import downloadGBKfromNCBI, getFastaForGene


genome_entries = 'genes/ureA/data/Koper_2004_seq.txt'
gbk_dir = 'genes/ureA/data/koper_2004_gbks/'
output_fasta = 'genes/ureA/data/koper_2004_seqs.fasta'


if __name__ == '__main__':
    
    print('* Downloading GenBank files...')
    entry_list = pd.read_csv(
        genome_entries, sep='\t', header=None,
        names=['species', 'entry']
        )
    downloadGBKfromNCBI(entry_list.entry.values, output_dir=gbk_dir)
    
    print('* Filtering sequences from database...')
    getFastaForGene(
    gbk_dir=gbk_dir,
    gene_keywords={
        "gene": ["ureC"],
        "product": ["urease", "alpha"]
    },
    output_fasta=output_fasta,
    case_insensitive=False
)