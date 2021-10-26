#!/usr/bin/env python
# conda activate traits

import os
import shutil
import argparse

import phyloplacement.wrappers as wrappers
from phyloplacement.utils import TemporaryFilePath
from phyloplacement.database.preprocessing import relabelRecordsInFASTA, getRepresentativeSet
from phyloplacement.database.manipulation import filterFASTAByHMM, countRecordsInFasta

"""
Reference database:
1) Run hmmer to extract peptides of interest
2) Reduce redundancy: cd-hit
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
                    default='',
                    help='Prefix to be added to output files')
parser.add_argument('--max_size', dest='maxsize',
                    default=None, type=int,
                    help=(
                        'Maximum size of representative set of sequences. '
                        'Defaults to full set.'
                        )
                    )
parser.add_argument('--relabel', dest='relabel', action='store_true',
                    default=False,
                    help='Relabel record IDs with numeral ids')


args = parser.parse_args()
output_fasta = os.path.join(args.outdir, f'{args.prefix}ref_database.faa')
output_fasta_short = os.path.join(args.outdir, f'{args.prefix}ref_database_short_ids.faa')
output_fasta_PI = os.path.join(args.outdir, f'{args.prefix}PI.txt')
reduced_fasta = os.path.join(args.outdir, f'{args.prefix}ref_database_reduced.faa')

def main():
    
    print('Making peptide-specific reference database...')
    with TemporaryFilePath() as tempaln:
        filterFASTAByHMM(
            hmm_model=args.hmm,
            input_fasta=args.data,
            output_fasta=output_fasta,
            additional_args=f'--cut_nc -A {tempaln}'
        )
        wrappers.getPercentIdentityFromMSA(
            input_msa=tempaln,
            output_file=output_fasta_PI
        )
    
    print('Finding representative sequences for reference database...')
    getRepresentativeSet(
        input_seqs=output_fasta,
        input_PI=output_fasta_PI,
        max_size=args.maxsize,
        outfile=reduced_fasta
    )
    shutil.move(reduced_fasta, output_fasta)
    
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
