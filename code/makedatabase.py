#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Reference database:
1) Run hmmer to extract peptides of interest
2) Reduce redundancy: cd-hit and/or repset
3) Relabel entries with temporary ids to avoid donwstream conflicts
"""

import os
import shutil
import argparse

from phyloplacement.utils import setDefaultOutputPath, TemporaryFilePath
from phyloplacement.database.preprocessing import setTempRecordIDsInFASTA
from phyloplacement.database.manipulation import filterFASTAByHMM, filterFastaBySequenceLength
from phyloplacement.database.reduction import reduceDatabaseRedundancy


parser = argparse.ArgumentParser(
    description='Build peptide reference database',
    epilog='Semidán Robaina Estévez (srobaina@ull.edu.es), 2021'
    )

optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
parser._action_groups.append(optional)

required.add_argument('--hmm', dest='hmm', type=str, required=True,
                      help='path to tigrfam or pfam model')
required.add_argument('--in', dest='data', type=str, required=True,
                      help='path to peptide database')
optional.add_argument('--outdir', dest='outdir', type=str,
                      help='path to output directory')
optional.add_argument('--prefix', dest='prefix', type=str,
                      default='',
                      help='prefix to be added to output files')
optional.add_argument('--max_size', dest='maxsize',
                      default=None, type=int,
                      help=(
                          'maximum size of representative set of sequences. '
                          'Defaults to full set.'
                          )
                    )
optional.add_argument('--min_seq_length', dest='minseqlength',
                      default=None, type=int,
                      help=(
                        'minimum sequence length in reference database. '
                        'Defaults to zero'
                        )
                    )
parser.add_argument('--max_seq_length', dest='maxseqlength',
                    default=None, type=int,
                    required=False,
                    help=(
                        'maximum sequence length in reference database. '
                        'Defaults to inf'
                        )
                    )
parser.add_argument('--relabel', dest='relabel', action='store_true',
                    required=False,
                    default=False,
                    help=(
                        'relabel record IDs with numeral ids. '
                        'Unrequired to build database, but highly recommended '
                        'to avoid possible conflicts downstream the pipeline.'))


args = parser.parse_args()
if args.outdir is None:
    args.outdir = setDefaultOutputPath(args.data, only_dirname=True)
hmmer_output = os.path.join(args.outdir, 'hmmer_output.txt')
output_fasta = os.path.join(args.outdir, f'{args.prefix}ref_database.faa')
output_fasta_short = os.path.join(args.outdir, f'{args.prefix}ref_database_short_ids.faa')

def main():
    
    print('Making peptide-specific reference database...')
    with TemporaryFilePath() as tempfasta, TemporaryFilePath() as tempfasta2:
        filterFASTAByHMM(
            hmm_model=args.hmm,
            input_fasta=args.data,
            output_fasta=tempfasta,
            hmmer_output=hmmer_output,
            additional_args=f'--cut_nc'
        )
        
        if (args.minseqlength is not None) or (args.maxseqlength is not None):
            print("Filtering sequences by established length bounds...")
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
        setTempRecordIDsInFASTA(
            input_fasta=output_fasta,
            output_dir=args.outdir,
            prefix=f'ref_{args.prefix}'
            )
        shutil.move(output_fasta_short, output_fasta)
    print('Finished!')

if __name__ == '__main__':
    main()
