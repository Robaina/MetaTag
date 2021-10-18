#!/usr/bin/env python
# conda activate traits

"""
Preprocessing:
0) Merge fasta files if more than one
1) Remove duplicates
2) Assert correct sequence format for downstream analysis
"""

import os
import argparse

from phyloplacement.utils import setDefaultOutputPath, createTemporaryFilePath
from phyloplacement.database.preprocessing import (assertCorrectSequenceFormat,
                                                   removeDuplicatesFromFasta,
                                                   mergeFASTAs)


parser = argparse.ArgumentParser(description='Database preprocessing')
parser.add_argument('--dna', dest='dna', type=bool, default=False,
                    help='Declare if sequences are nucleotides. Default to peptide sequences.')
parser.add_argument('--in', dest='data', type=str,
                    help='Path to fasta file or directory containing fasta files')
parser.add_argument('--out', dest='out', type=str, default=None,
                    help='Path to output directory')

args = parser.parse_args()
is_peptide = not args.dna
if args.out is None:
    outfasta = setDefaultOutputPath(args.data, tag='_cleaned')
else:
    outfasta = args.out
if os.path.isdir(args.data):
    data_path = os.path.join(args.data, 'merged_data.fa')
    mergeFASTAs(
        input_fastas_dir=args.data,
        output_fasta=data_path
    )
else:
    data_path = args.data

def main():
    
    tmp_file_path = createTemporaryFilePath()

    removeDuplicatesFromFasta(
        input_fasta=data_path,
        output_fasta=tmp_file_path
    )
    assertCorrectSequenceFormat(
        fasta_file=tmp_file_path,
        output_file=outfasta,
        is_peptide=is_peptide
    )

    os.remove(tmp_file_path)

if __name__ == '__main__':
    main()
