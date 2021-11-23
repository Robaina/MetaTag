#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Preprocessing:
0) Merge fasta files if more than one
1) Remove duplicates
2) Assert correct sequence format for downstream analysis
3) Translate if nucleotide
4) Relabel if asked for
"""

import os
import shutil
import argparse

from phyloplacement.utils import (setDefaultOutputPath,
                                  TemporaryFilePath,
                                  TemporaryDirectoryPath)
from phyloplacement.wrappers import runProdigal
from phyloplacement.database.preprocessing import (assertCorrectSequenceFormat,
                                                   removeDuplicatesFromFasta,
                                                   mergeFASTAs,
                                                   setTempRecordIDsInFASTA)


parser = argparse.ArgumentParser(
    description='Database preprocessing: removal of duplicated sequences and of sequences with illegal symbols',
    epilog='Semidán Robaina Estévez (srobaina@ull.edu.es), 2021'
    )

optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
parser._action_groups.append(optional)

required.add_argument('--in', dest='data', type=str, required=True,
                      help='path to fasta file or directory containing fasta files')
optional.add_argument('--dna', dest='dna', action='store_true', default=False,
                      help='declare if sequences are nucleotides. Defaults to peptide sequences.')
optional.add_argument('--translate', dest='translate', action='store_true', default=False,
                      help='choose whether nucleotide sequences are translated with prodigal')
optional.add_argument('--outfile', dest='outfile', type=str,
                      help='path to output fasta file')
optional.add_argument('--relabel', dest='relabel',
                      default=False, action='store_true',
                      help='relabel record IDs with numeral ids')
optional.add_argument('--idprefix', dest='idprefix', type=str,
                      default='',
                      help='prefix to be added to sequence IDs')

args = parser.parse_args()
is_peptide = not args.dna
output_dir = os.path.abspath(os.path.dirname(args.outfile))

if os.path.isdir(args.data):
    print('Merging input files...')
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
    
    with TemporaryFilePath() as tmp_file_path:
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
    
    if args.translate:
        outprefix = setDefaultOutputPath(outfasta, only_filename=True)
        with TemporaryDirectoryPath() as tempdir:
            runProdigal(
                input_file=outfasta,
                output_prefix=outprefix,
                output_dir=tempdir, 
                metagenome=False, 
                additional_args=None
            )
            outfaa = os.path.join(tempdir, outprefix + '.faa')
            outgbk = os.path.join(tempdir, outprefix + '.gbk')
            shutil.move(outfaa, outfasta)
            shutil.move(outgbk, output_dir)

    if args.relabel:
        print('Relabelling records...')
        outfasta_short = setDefaultOutputPath(outfasta, tag='_short_ids')
        setTempRecordIDsInFASTA(
            input_fasta=outfasta,
            output_dir=os.path.dirname(args.outfile),
            prefix=args.idprefix
            )
        shutil.move(outfasta_short, outfasta)
    print('Finished!')

if __name__ == '__main__':
    main()
