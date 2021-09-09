"""
Classify nifH sequences in clusters based on CART model,
CART model from: https://sfamjournals.onlinelibrary.wiley.com/doi/10.1111/1758-2229.12455

Usage: python3 nifhcart.py input_fasta input_alignment [output_fasta]
"""

import os
import sys
from Bio import SeqIO


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
            record.id.split('_')[0]: list(record.seq)
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
    azo_nifH = getRecordAlignments(input_fasta)['001']
    azo_nifh_aln = ''.join(getRecordAlignments(input_alignment)['001'])

    for aa_pos, aas in CART.items():
        pos_pattern = ''.join(azo_nifH)[aa_pos - 1: aa_pos + 9]
        adj_aa_pos = azo_nifh_aln.find(pos_pattern)
        adjusted_CART[adj_aa_pos] = aas

    return adjusted_CART

def addClusterToNifH(input_fasta: str, input_alignment: str,
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


if __name__ == '__main__':

    input_fasta = sys.argv[1]
    input_alignment = sys.argv[2]
    if len(sys.argv) > 3:
        output_fasta = sys.argv[3]
    else:
        output_fasta = None

    addClusterToNifH(input_fasta, input_alignment, output_fasta)
    