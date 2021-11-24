#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Classify nifH sequences in clusters based on CART model,
CART model from: https://sfamjournals.onlinelibrary.wiley.com/doi/10.1111/1758-2229.12455
"""

import os
import argparse

from Bio import SeqIO

from phyloplacement.utils import readFromPickleFile, saveToPickleFile


def getnifHclusterID(seq: list, cart_model=dict) -> str:
    """
    Assign cluster to nifH sequence based on CART model in
    https://sfamjournals.onlinelibrary.wiley.com/doi/10.1111/1758-2229.12455

    """
    aa_pos = list(cart_model.keys())
    if len(seq) < max(aa_pos):
        return '0'
    if seq[aa_pos[0]] in cart_model[aa_pos[0]]:
        return 'I'
    elif seq[aa_pos[1]] in cart_model[aa_pos[1]]:
        return 'II'
    elif seq[aa_pos[2]] in cart_model[aa_pos[2]]:
        return 'III'
    else:
        return 'IV'

def getRecordAlignments(input_alignment: str) -> dict:
    """
    Get dict of lists containing sequence alignments
    """
    with open(input_alignment) as inalign:
        return {
            record.id: list(record.seq)
            for record in SeqIO.parse(inalign, 'fasta')
        }

def adjustCARTmodel(input_fasta: str, input_alignment: str):
    """
    Adjust CART aminoacid positions to current alignment.
    CART based on nifH sequence from Azotobacter
    """
    CART = {
        109: ['F', 'W', 'Y'], 
        49: ['A', 'D', 'I'],  
        53: ['L', 'M', 'W'] 
    }
    adjusted_CART = {}
    seq_dict = getRecordAlignments(input_fasta)
    first_id = list(seq_dict.keys())[0]
    aln_dict = getRecordAlignments(input_alignment)
    azo_nifH = seq_dict[first_id]
    azo_nifh_aln = ''.join(aln_dict[first_id])

    for aa_pos, aas in CART.items():
        pos_pattern = ''.join(azo_nifH)[aa_pos - 1: aa_pos + 9]
        adj_aa_pos = azo_nifh_aln.find(pos_pattern)
        adjusted_CART[adj_aa_pos] = aas

    return adjusted_CART

def addClusterToNifHfasta(input_fasta: str, input_alignment: str,
                     output_fasta: str = None) -> None:
    """
    Add assigned cluster to nifH sequences in fasta file
    """
    input_fasta = os.path.abspath(input_fasta)
    input_alignment = os.path.abspath(input_alignment)
    recordAligns = getRecordAlignments(input_alignment)
    cart_model = adjustCARTmodel(input_fasta, input_alignment)

    if output_fasta is None:
        base, ext = os.path.splitext(input_fasta)
        output_fasta = f'{base}_clustered{ext}'
    else:
        output_fasta = os.path.abspath(output_fasta)

    with open(input_fasta) as infasta, open(output_fasta, 'w') as outfasta:
        for record in SeqIO.parse(infasta, 'fasta'):
            record_align_seq = recordAligns[record.id.split('_')[0]]
            cluster_id = getnifHclusterID(record_align_seq, cart_model)
            record.id = f'{record.id}_cluster_{cluster_id}'
            record.name = ''
            record.description = ''
            SeqIO.write(record, outfasta, 'fasta')

def addClusterToNifHdict(input_fasta: str,
                         input_alignment: str,
                         input_dict: str,
                         output_dict: str = None) -> None:
    """
    Add assigned cluster to nifH sequence alignments
    in reference ID dictionary
    """
    input_fasta = os.path.abspath(input_fasta)
    input_alignment = os.path.abspath(input_alignment)
    input_dict = os.path.abspath(input_dict)
    cart_model = adjustCARTmodel(input_fasta, input_alignment)
    ref_dict = readFromPickleFile(input_dict)

    if output_dict is None:
        output_dict = input_dict
    else:
        output_dict = os.path.abspath(output_dict)

    for record in SeqIO.parse(input_alignment, 'fasta'):
        cluster_id = getnifHclusterID(list(record.seq), cart_model)
        ref_dict[record.id] += f'_cluster_{cluster_id}'
    
    saveToPickleFile(ref_dict, output_dict)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
    description='Classify nifH peptide sequences in clusters based on CART model',
    epilog='Semidán Robaina Estévez (srobaina@ull.edu.es), 2021'
    )

    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    parser._action_groups.append(optional)

    required.add_argument('--seqs', dest='seqs', type=str, required=True,
                        help='path to fasta file')
    required.add_argument('--aln', dest='aln', type=str, required=True,
                        help='path to fasta file')
    required.add_argument('--indict', dest='indict', type=str, required=True,
                        help='path to input dictionary with reference IDs and labels')
    optional.add_argument('--outdict', dest='outdict', type=str,
                        help='path to output ID dictionary in pickle format')

    args = parser.parse_args()
    if args.outfile is None:
        base, ext = os.path.splitext(args.indict)
        outdict = os.path.abspath(base + '_clustered' + ext)
    else:
        outdict = os.path.abspath(args.outdict)


    addClusterToNifHdict(args.seqs, args.aln, args.indict, outdict)
    