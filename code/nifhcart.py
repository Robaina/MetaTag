"""
Classify nifH sequences in clusters based on CART model,
CART model from: https://sfamjournals.onlinelibrary.wiley.com/doi/10.1111/1758-2229.12455

Usage: python3 nifhcart.py input_fasta [output_fasta]
"""

import os
import sys
from Bio import SeqIO


def getnifHclusterID(seq: list) -> str:
    """
    Assign cluster to nifH sequence based on CART model in
    https://sfamjournals.onlinelibrary.wiley.com/doi/10.1111/1758-2229.12455
    """
    CART = {
        109: ['F', 'W', 'Y'],
        49: ['A', 'D', 'I'],
        53: ['L', 'M', 'W']
    }
    aa_pos = list(CART.keys())
    if len(seq) < max(aa_pos):
        return '0'
    if seq[aa_pos[0] - 1] in CART[aa_pos[0]]:
        return 'I'
    elif seq[aa_pos[1] - 1] in CART[aa_pos[1]]:
        return 'II'
    elif seq[aa_pos[2] - 1] in CART[aa_pos[2]]:
        return 'III'
    else:
        return 'IV'

def addClusterToNifH(input_fasta: str, output_fasta: str = None) -> None:
    """
    Add assigned cluster to nifH sequences in fasta file
    """
    input_fasta = os.path.abspath(input_fasta)
    if output_fasta is None:
        base, ext = os.path.splitext(input_fasta)
        output_fasta = f'{base}_clustered{ext}'
    else:
        output_fasta = os.path.abspath(output_fasta)

    with open(input_fasta) as infasta, open(output_fasta, 'w') as outfasta:
        for record in SeqIO.parse(infasta, 'fasta'):
            cluster_id = getnifHclusterID(list(record.seq))
            record.id = f'{record.id}_cluster_{cluster_id}'
            record.name = ''
            record.description = ''
            SeqIO.write(record, outfasta, 'fasta')


if __name__ == '__main__':

    input_fasta = sys.argv[1]
    if len(sys.argv) > 2:
        output_fasta = sys.argv[2]
    else:
        output_fasta = None

    addClusterToNifH(input_fasta, output_fasta)
    