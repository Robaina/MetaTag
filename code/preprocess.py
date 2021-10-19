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

from phyloplacement.utils import (setDefaultOutputPath,
                                  createTemporaryFilePath)
from phyloplacement.database.preprocessing import (assertCorrectSequenceFormat,
                                                   removeDuplicatesFromFasta,
                                                   mergeFASTAs)


parser = argparse.ArgumentParser(description='Database preprocessing')
parser.add_argument('--dna', dest='dna', type=bool, default=False,
                    help='Declare if sequences are nucleotides. Default to peptide sequences.')
parser.add_argument('--in', dest='data', type=str,
                    help='Path to fasta file or directory containing fasta files')
parser.add_argument('--outfile', dest='outfile', type=str, default=None,
                    help='Path to output fasta file')

args = parser.parse_args()
is_peptide = not args.dna
if args.outfile is None:
    outfasta = setDefaultOutputPath(args.data, tag='_cleaned')
else:
    outfasta = args.outfile
if os.path.isdir(args.data):
    print('Merging input files...')
    output_dir = os.path.dirname(args.outfile)
    _, file_ext = os.path.splitext(os.listdir(args.data)[0])
    data_path = os.path.abspath(os.path.join(output_dir, f'merged_data{file_ext}'))
    mergeFASTAs(
        input_fastas_dir=args.data,
        output_fasta=data_path
    )
else:
    data_path = args.data

def main():
    
    tmp_file_path = createTemporaryFilePath()
    
    print('Removing duplicates...')
    removeDuplicatesFromFasta(
        input_fasta=data_path,
        output_fasta=tmp_file_path
    )
    print('Asserting correct sequence format...')
    assertCorrectSequenceFormat(
        fasta_file=tmp_file_path,
        output_file=outfasta,
        is_peptide=is_peptide
    )

    os.remove(tmp_file_path)
    print('Finished!')

if __name__ == '__main__':
    main()
