#!/usr/bin/env python
# conda activate traits

"""
'/home/robaina/Documents/TRAITS/tests/ref_reduced_clean_id_dict.pickle'
'/home/robaina/Documents/MAR_database/mardb_assembly_V6.fa'
'/home/robaina/Documents/TRAITS/data/nxr/nxr_genomes'
"""

import os
import argparse

from phyloplacement.utils import parallelizeOverInputFiles, readFromPickleFile
from phyloplacement.database.mardb import getMARdbGenomeByEntryCode, getMarDBentryCode


parser = argparse.ArgumentParser(description='Get MarDB genomes given set of sequences')
parser.add_argument('--data', dest='data', type=str,
                    help='Path mardb assembly fasta file')
parser.add_argument('--ids', dest='ids', type=str,
                    help='Path to pickle with original query sequence ids')
parser.add_argument('--out', dest='outdir', type=str,
                    help='Path to output directory')
parser.add_argument('--proc', dest='proc', type=int,
                    help='Number of processes to use')                  
args = parser.parse_args()

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

def parallel_genome(input_id, input_fasta: str, output_dir: str):
    output_fasta = os.path.join(output_dir, f'{input_id}.fa')
    getMARdbGenomeByEntryCode(input_id,
                              input_fasta=input_fasta,
                              output_fasta=output_fasta,
                              clean=True)

# Retrieve mardb genomes
label_dict = readFromPickleFile(args.ids)
nxr_entry_codes = {getMarDBentryCode(v) for v in label_dict.values()}

parallelizeOverInputFiles(
    parallel_genome,
    input_list=nxr_entry_codes,
    n_processes=args.proc,
    input_fasta=args.data,
    output_dir=args.outdir
)
