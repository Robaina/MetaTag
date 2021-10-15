#!/usr/bin/env python
# conda activate traits

import os
import shutil
import argparse
import phyloplacement.wrappers as wrappers
from phyloplacement.database.preprocessing import relabelRecordsInFASTA
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
parser.add_argument('--out', dest='outdir', type=str,
                    help='Path to output directory')
parser.add_argument('--reduce', dest='reduce', type=bool,
                    default=False,
                    help='Run cd-hit to reduce database redundancy')
args = parser.parse_args()

output_fasta = os.path.join(args.outdir, 'ref_database.faa')
reduced_fasta = os.path.join(args.outdir, 'ref_database_reduced.faa')

def main():
    
    # Make peptide-specific database
    filterFASTAByHMM(
        hmm_model=args.hmm,
        input_fasta=args.data,
        output_fasta=output_fasta
    )
    
    if args.reduce:
        # Reduce redundancy of reference database
        wrappers.runCDHIT(
            input_fasta=output_fasta,
            output_fasta=reduced_fasta,
            additional_args=None
            )
        n_records = countRecordsInFasta(output_fasta)
        n_reduced_records = countRecordsInFasta(reduced_fasta)
        shutil.move(reduced_fasta, output_fasta)
        print(f'Original database size: {n_records}. Reduced database size: {n_reduced_records}')

    # Assign numbers to reference sequence labels for data processing
    relabelRecordsInFASTA(
        input_fasta=output_fasta,
        output_dir=args.outdir,
        prefix='ref_'
        )

if __name__ == '__main__':
    main()
