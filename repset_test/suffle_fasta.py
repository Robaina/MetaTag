#!/bin/env python

import os
import argparse
import random
from Bio import SeqIO 


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Random shuffle records in fasta")
    parser.add_argument("--infasta", type=str, required=True, help="Input fasta directory")
    parser.add_argument("--outfasta", type=str, required=False, default=None, help="Output fasta directory")
    parser.add_argument("--size", type=int, required=False, default=None, help="Size of randomized record set")
    args = parser.parse_args()

    if args.outfasta is None:
        fname, ext = os.path.splitext(args.infasta)
        args.outfasta = fname + '_shuffled' + ext

    records = list(SeqIO.parse(args.infasta, 'fasta'))
    random.shuffle(records)

    if args.size is not None:
        if args.size > len(records):
            raise ValueError("Output size can't be larger than input size")
        records = records[:args.size]
        
    SeqIO.write(records, args.outfasta, 'fasta')