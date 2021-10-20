#!/usr/bin/env python
# conda activate traits

"""
Preprocessing:
0) Merge fasta files if more than one
1) Remove duplicates
2) Assert correct sequence format for downstream analysis
"""

import os
import shutil
import argparse

from phyloplacement.utils import (setDefaultOutputPath,
                                  createTemporaryFilePath)
from phyloplacement.database.preprocessing import (assertCorrectSequenceFormat,
                                                   removeDuplicatesFromFasta,
                                                   mergeFASTAs,
                                                   relabelRecordsInFASTA)


parser = argparse.ArgumentParser(description='Database preprocessing')
parser.add_argument('--dna', dest='dna', action='store_true', default=False,
                    help='Declare if sequences are nucleotides. Default to peptide sequences.')
parser.add_argument('--in', dest='data', type=str,
                    help='Path to fasta file or directory containing fasta files')
parser.add_argument('--outfile', dest='outfile', type=str, default=None,
                    help='Path to output fasta file')
parser.add_argument('--is_query', dest='isquery',
                    default=False, action='store_true',
                    help='Indicate if fasta file contains query sequences. Defaults to reference database.')

args = parser.parse_args()
is_peptide = not args.dna

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
    data_path = os.path.abspath(args.data)

if args.outfile is None:
    outfasta = setDefaultOutputPath(data_path, tag='_cleaned')
else:
    outfasta = os.path.abspath(args.outfile)

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

    if args.isquery:
        print('Relabelling records...')
        outfasta_short = setDefaultOutputPath(outfasta, tag='_short_ids')
        relabelRecordsInFASTA(
            input_fasta=outfasta,
            output_dir=os.path.dirname(args.outfile),
            prefix='query_'
            )
        shutil.move(outfasta_short, outfasta)
    print('Finished!')

if __name__ == '__main__':
    main()
