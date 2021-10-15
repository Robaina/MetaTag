#!/usr/bin/env python
# conda activate traits

"""
Preprocessing:
1) Remove duplicates
2) Assert correct sequence format for downstream analysis
"""

import tempfile
import argparse
from phyloplacement.utils import setDefaultOutputPath
from phyloplacement.database.preprocessing import (assertCorrectSequenceFormat,
                                                   removeDuplicatesFromFasta)


parser = argparse.ArgumentParser(description='Database preprocessing')
parser.add_argument('--dna', dest='dna', type=bool, default=False,
                    help='Declare if sequences are nucleotides. Default to peptide sequences.')
parser.add_argument('--in', dest='data', type=str,
                    help='Path to fasta file')
parser.add_argument('--out', dest='out', type=str, default=None,
                    help='Path to output directory')

args = parser.parse_args()
is_peptide = not args.dna
if args.out is None:
    outfasta = setDefaultOutputPath(args.data, tag='_cleaned')
else:
    outfasta = args.out

def main():
    
    with tempfile.NamedTemporaryFile(mode='w') as tfile:
        removeDuplicatesFromFasta(
            input_fasta=args.data,
            output_fasta=tfile.name
        )
        assertCorrectSequenceFormat(
            fasta_file=tfile.name,
            output_file=outfasta,
            is_peptide=is_peptide
        )

if __name__ == '__main__':
    main()
