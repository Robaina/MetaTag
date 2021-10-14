#!/usr/bin/env python
# conda activate traits

import argparse
from pathlib import Path

import phyloplacement.wrappers as wrappers
from phyloplacement.database.preprocessing import relabelRecordsInFASTA, assertCorrectSequenceFormat
from phyloplacement.database.manipulation import filterFASTAByHMM, countRecords

"""
Preprocessing:
1) Assert correct sequence format for downstream analysis
2) Reduce redundancy: remove duplicates, get representatives with cd-hit
3) Relabel entries with temporary ids to avoid donwstream conflicts
"""

parser = argparse.ArgumentParser(description='Database preprocessing')
parser.add_argument('--hmm', dest='hmm', type=str,
                    help='Path to tigrfam or pfam model')
parser.add_argument('--in', dest='data', type=str,
                    help='Path to peptide database')
parser.add_argument('--out', dest='out', type=str,
                    help='Path to output directory')
args = parser.parse_args()

"""
TIGR00639.1.HMM (large peptide database after cd-hit)
TIGR01580.1.HMM (narG / nxr)
"""
mardb_data = Path('/home/robaina/Documents/MAR_database/')
tigr_data = Path('/home/robaina/Documents/tigrfams/hmm_PGAP/')
work_dir = Path('/home/robaina/Documents/TRAITS/tests/')

def main():
    # Make peptide-specific database
    filterFASTAByHMM(
        hmm_model=str(tigr_data / 'TIGR01580.1.HMM'),
        input_fasta=str(mardb_data / 'mardb_proteins_V6_no_duplicates.fasta'),
        output_fasta=str(work_dir / 'mardb_TIGR01580.1.fasta')
    )

    # 1) Reduce redundancy of reference database
    wrappers.runCDHIT(
        input_fasta=str(work_dir / 'mardb_TIGR01580.1.fasta'),
        output_fasta=str(work_dir / 'ref_reduced.fasta'),
        additional_args=None
        )
    n_records = countRecords(str(work_dir / 'mardb_TIGR01580.1.fasta'))
    n_reduced_records = countRecords(str(work_dir / 'ref_reduced.fasta'))

    # 2) Assert  correct format
    assertCorrectSequenceFormat(
        fasta_file=str(work_dir / 'ref_reduced.fasta'),
        output_file=str(work_dir / 'ref_reduced_clean.fasta'),
    )

    # 3) Assign numbers to reference sequence labels for data processing
    relabelRecordsInFASTA(
        input_fasta=str(work_dir / 'ref_reduced_clean.fasta'),
        output_dir=str(work_dir),
        prefix='ref_'
        )

if __name__ == '__main__':
    main()
