#!/usr/bin/env python
# conda activate traits

"""
Preprocessing:
1) Assert correct sequence format for downstream analysis
2) Reduce redundancy: remove duplicates, get representatives with cd-hit
3) Relabel entries with temporary ids to avoid donwstream conflicts
"""

import argparse
from pathlib import Path

import phyloplacement.wrappers as wrappers
from phyloplacement.database.preprocessing import relabelRecordsInFASTA, assertCorrectSequenceFormat
from phyloplacement.database.manipulation import filterFASTAByHMM, countRecords


parser = argparse.ArgumentParser(description='Database preprocessing')
parser.add_argument('--hmm', dest='hmm', type=str,
                    help='Path to tigrfam or pfam model')
parser.add_argument('--data', dest='data', type=str,
                    help='Path to peptide database')
parser.add_argument('--out', dest='outdir', type=str,
                    help='Path to output directory')
args = parser.parse_args()



def main():
    # Make peptide-specific database
    filterFASTAByHMM(
        hmm_model=args.hmm,
        input_fasta=args.data,
        output_fasta=str(args.outdir / 'mardb_TIGR01580.1.fasta')
    )

    # 1) Reduce redundancy of reference database
    wrappers.runCDHIT(
        input_fasta=str(args.outdir / 'mardb_TIGR01580.1.fasta'),
        output_fasta=str(args.outdir / 'ref_reduced.fasta'),
        additional_args=None
        )
    n_records = countRecords(str(args.outdir / 'mardb_TIGR01580.1.fasta'))
    n_reduced_records = countRecords(str(args.outdir / 'ref_reduced.fasta'))

    # 2) Assert  correct format
    assertCorrectSequenceFormat(
        fasta_file=str(args.outdir / 'ref_reduced.fasta'),
        output_file=str(args.outdir / 'ref_reduced_clean.fasta'),
    )

    # 3) Assign numbers to reference sequence labels for data processing
    relabelRecordsInFASTA(
        input_fasta=str(args.outdir / 'ref_reduced_clean.fasta'),
        output_dir=str(args.outdir),
        prefix='ref_'
        )

if __name__ == '__main__':
    main()
