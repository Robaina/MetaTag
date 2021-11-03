#!/usr/bin/env python
# conda activate traits

import os
import shutil
import argparse

from phyloplacement.utils import TemporaryFilePath
from phyloplacement.database.preprocessing import relabelRecordsInFASTA
from phyloplacement.database.manipulation import filterFASTAByHMM, filterFastaBySequenceLength
from phyloplacement.database.reduction import reduceDatabaseRedundancy

"""
Reference database:
1) Run hmmer to extract peptides of interest
2) Reduce redundancy: cd-hit and/or repset
3) Relabel entries with temporary ids to avoid donwstream conflicts
"""

parser = argparse.ArgumentParser(description='Build peptide reference database')
parser.add_argument('--hmm', dest='hmm', type=str,
                    help='Path to tigrfam or pfam model')
parser.add_argument('--in', dest='data', type=str,
                    help='Path to peptide database')
parser.add_argument('--outdir', dest='outdir', type=str,
                    help='Path to output directory')
parser.add_argument('--prefix', dest='prefix', type=str,
                    required=False,
                    default='',
                    help='Prefix to be added to output files')
parser.add_argument('--max_size', dest='maxsize',
                    required=False,
                    default=None, type=int,
                    help=(
                        'Maximum size of representative set of sequences. '
                        'Defaults to full set.'
                        )
                    )
parser.add_argument('--min_seq_length', dest='minseqlength',
                    default=None, type=int,
                    required=False,
                    help=(
                        'Minimum sequence length in reference database. '
                        'Defaults to zero'
                        )
                    )
parser.add_argument('--max_seq_length', dest='maxseqlength',
                    default=None, type=int,
                    required=False,
                    help=(
                        'Maximum sequence length in reference database. '
                        'Defaults to inf'
                        )
                    )
parser.add_argument('--relabel', dest='relabel', action='store_true',
                    required=False,
                    default=False,
                    help='Relabel record IDs with numeral ids')


args = parser.parse_args()
output_fasta = os.path.join(args.outdir, f'{args.prefix}ref_database.faa')
output_fasta_short = os.path.join(args.outdir, f'{args.prefix}ref_database_short_ids.faa')

def main():
    
    print('Making peptide-specific reference database...')
    with TemporaryFilePath() as tempfasta, TemporaryFilePath() as tempfasta2:
        filterFASTAByHMM(
            hmm_model=args.hmm,
            input_fasta=args.data,
            output_fasta=tempfasta,
            additional_args=f'--cut_nc'
        )
        
        if (args.minseqlength is not None) or (args.maxseqlength is not None):
            print("Filtering sequences by established length bounds")
            filterFastaBySequenceLength(
                input_fasta=tempfasta,
                minLength=args.minseqlength,
                maxLength=args.maxseqlength,
                output_fasta=tempfasta2
            )
            shutil.move(tempfasta2, tempfasta)
    
        reduceDatabaseRedundancy(
            input_fasta=tempfasta,
            output_fasta=output_fasta,
            cdhit=True,
            cdhit_args=None,
            maxsize=args.maxsize
        )

    if args.relabel:
        print('Relabelling records in reference database...')
        relabelRecordsInFASTA(
            input_fasta=output_fasta,
            output_dir=args.outdir,
            prefix=f'ref_{args.prefix}'
            )
        shutil.move(output_fasta_short, output_fasta)
    print('Finished!')

if __name__ == '__main__':
    main()
