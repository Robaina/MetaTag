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

from phyloplacement.utils import TemporaryDirectoryPath, fullPathListDir, setDefaultOutputPath, TemporaryFilePath, DictMerger
from phyloplacement.database.preprocessing import setTempRecordIDsInFASTA, mergeFASTAs
from phyloplacement.database.manipulation import filterFASTAByHMM, filterFastaBySequenceLength
from phyloplacement.database.reduction import reduceDatabaseRedundancy


parser = argparse.ArgumentParser(
    description='Build peptide reference database',
    epilog='Semidán Robaina Estévez (srobaina@ull.edu.es), 2021'
    )

optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
parser._action_groups.append(optional)

required.add_argument('--hmms', dest='hmms', type=str, nargs='+',
                      required=True,
                      help='a single or a list of space-separated paths to tigrfam or pfam models')
required.add_argument('--in', dest='data', type=str, required=True,
                      help='path to peptide database')
optional.add_argument('--outdir', dest='outdir', type=str,
                      help='path to output directory')
optional.add_argument('--prefix', dest='prefix', type=str,
                      default='',
                      help='prefix to be added to output files')
optional.add_argument('--relabel_prefixes', dest='relabel_prefixes', type=str,
                      nargs='+',
                      help=(
                          'List of space-separated prefixes to be added to each hmm-derived set of '
                          'sequences after relabelling. Only used if "--relabel" is set'
                      )
                      )
optional.add_argument('--max_sizes', dest='maxsizes',
                      default=None, type=str, nargs='+',
                      help=(
                          'maximum size of representative set of sequences for each hmm model. '
                          'Each (space-separated) integer corresponds to a hmm model inputed in "--hmms", '
                          'thus, sorted in the same order. A value of "None" may be given to a hmm model '
                          'in the list, in which case the maximum number of sequences is unlimited '
                          'for that hmm.'
                          'Defaults to full set of sequences for all hmm modells inputed.'
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
if args.maxsizes is None:
    args.maxsizes = [None for _ in args.hmms]
else:
    args.maxsizes = [int(arg) if arg.isdigit() else None for arg in args.maxsizes]
if args.relabel_prefixes is None:
    args.relabel_prefixes = [None for _ in args.hmms]
if args.outdir is None:
    args.outdir = setDefaultOutputPath(args.data, only_dirname=True)
args.outdir = os.path.abspath(args.outdir)
if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)
output_fasta = os.path.join(args.outdir, f'{args.prefix}ref_database.faa')
output_pickle_short_ids = os.path.join(args.outdir, f'{args.prefix}ref_database_id_dict.pickle')


def main():
    
    print('* Making peptide-specific reference database...')
    with TemporaryDirectoryPath() as tempdir1, TemporaryDirectoryPath() as tempdir2:
        for hmm, maxsize, prefix in zip(args.hmms, args.maxsizes, args.relabel_prefixes):
            hmm_name = os.path.basename(hmm)
            print(f" * Processing hmm {hmm_name}")
            hmmer_output = os.path.join(args.outdir, f"hmmer_output_{hmm_name}.txt")

            with TemporaryFilePath() as tempfasta, TemporaryFilePath() as tempfasta2, TemporaryFilePath() as tempfasta3:
                filterFASTAByHMM(
                    hmm_model=hmm,
                    input_fasta=args.data,
                    output_fasta=tempfasta,
                    hmmer_output=hmmer_output,
                    additional_args='--cut_nc'
                )
                
                if (args.minseqlength is not None) or (args.maxseqlength is not None):
                    print("* Filtering sequences by established length bounds...")
                    filterFastaBySequenceLength(
                        input_fasta=tempfasta,
                        minLength=args.minseqlength,
                        maxLength=args.maxseqlength,
                        output_fasta=tempfasta2
                    )
                    shutil.move(tempfasta2, tempfasta)
            
                reduceDatabaseRedundancy(
                    input_fasta=tempfasta,
                    output_fasta=tempfasta3,
                    cdhit=True,
                    cdhit_args=None,
                    maxsize=maxsize
                )

                if args.relabel:
                    print('* Relabelling records in reference database...')
                    output_fasta_short = os.path.join(tempdir2, f"{tempfasta3}_short_ids")
                    setTempRecordIDsInFASTA(
                        input_fasta=tempfasta3,
                        output_dir=tempdir2,
                        prefix=prefix 
                        )
                    shutil.move(output_fasta_short, tempdir1)

        mergeFASTAs(
            input_fastas_dir=tempdir1,
            output_fasta=output_fasta
        )
         
        pickle_dict_paths = [file for file in fullPathListDir(tempdir2) if file.endswith('.pickle')]
        if pickle_dict_paths:
            DictMerger.fromPicklePaths(pickle_dict_paths).merge(save_pickle_path=output_pickle_short_ids)

    print('Finished!')

if __name__ == '__main__':
    main()
