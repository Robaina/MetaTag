#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Classify nifH sequences in clusters based on CART model,
CART model from: https://sfamjournals.onlinelibrary.wiley.com/doi/10.1111/1758-2229.12455
"""

import argparse
import os

from Bio import SeqIO, pairwise2

from metatag.utils import read_from_pickle_file, save_to_pickle_file


def find_pattern_in_msa_record(msa_record: str, subsequence: str) -> int:
    """
    Find index of the first character of a match between a
    subsequence and a sequence alignment record
    """
    alns = pairwise2.align.globalxx(msa_record, subsequence)
    ind = list(set([a.seqB.index(subsequence[0]) for a in alns]))
    if len(ind) > 1:
        raise ValueError("More than one match found for subsequence")
    elif not ind:
        raise ValueError("No match for subsequence")
    elif len(ind) == 1:
        print(f"Only one for: {subsequence}")
        return ind[0]


def get_nifh_cluster_id(seq: list, cart_model: dict) -> str:
    """
    Assign cluster to nifH sequence based on CART model in
    https://sfamjournals.onlinelibrary.wiley.com/doi/10.1111/1758-2229.12455

    """
    aa_pos = list(cart_model.keys())
    if len(seq) < max(aa_pos):
        return "0"
    if seq[aa_pos[0]] in cart_model[aa_pos[0]]:
        return "I"
    elif seq[aa_pos[1]] in cart_model[aa_pos[1]]:
        return "II"
    elif seq[aa_pos[2]] in cart_model[aa_pos[2]]:
        return "III"
    else:
        return "IV"


def get_record_alignments(input_alignment: str) -> dict:
    """
    Get dict of lists containing sequence alignments
    """
    with open(input_alignment) as inalign:
        return {record.id: list(record.seq) for record in SeqIO.parse(inalign, "fasta")}


def adjust_cart_model(
    input_fasta: str,
    input_alignment: str,
    azo_id: str = "ref_azo",
    substring_length: int = 25,
):
    """
    Adjust CART aminoacid positions to current alignment.
    CART based on nifH sequence from Azotobacter.
    Using global alignment of subsequence to find position in aligment.
    """
    # CART = {
    #     109: ['F', 'W', 'Y'],
    #     49: ['A', 'D', 'I'],
    #     53: ['L', 'M', 'W']
    # }

    # adjusted_CART = {}
    # seq_dict = get_record_alignments(input_fasta)
    # aln_dict = get_record_alignments(input_alignment)
    # azo_nifH = seq_dict[azo_id]
    # azo_nifh_aln = ''.join(aln_dict[azo_id])

    # for aa_pos, aas in CART.items():
    #     pos_pattern = ''.join(azo_nifH)[aa_pos - 1: aa_pos - 1 + substring_length]
    #     aln_aa_pos = find_pattern_in_msa_record(azo_nifh_aln, pos_pattern)
    #     adjusted_CART[aln_aa_pos] = aas

    """
    NOTE: Adjusted manually based on msa containing Azotobacter's nifH sequence:
    001_WP_039801084.1 MULTISPECIES: nitrogenase iron protein [Azotobacter]
    and counting residue positions in the aligned Azotobacter sequence vs residues
    in the unaligned Azotobacter sequence.
    """
    # Adjusted model for original nifH tree built with --cut_nc
    # adjusted_CART = {
    #     (180 - 1): ['F', 'W', 'Y'],
    #     (95 - 1): ['A', 'D', 'I'],
    #     (99 - 1): ['L', 'M', 'W']
    # }
    # # Adjusted CART built for alternative tree built with -E 1e-10
    # adjusted_CART = {
    #     (561 - 1): ['F', 'W', 'Y'],
    #     (398 - 1): ['A', 'D', 'I'],
    #     (405 - 1): ['L', 'M', 'W']
    # }
    # Adjusted CART built for new tree built by Nuria
    adjusted_CART = {
        (275 - 1): ["F", "W", "Y"],
        (115 - 1): ["A", "D", "I"],
        (147 - 1): ["L", "M", "W"],
    }
    print(f"CART adjusted to {adjusted_CART}")
    return adjusted_CART


def add_cluster_to_nifh_fasta(
    input_fasta: str, input_alignment: str, output_fasta: str = None
) -> None:
    """
    Add assigned cluster to nifH sequences in fasta file
    """
    input_fasta = os.path.abspath(input_fasta)
    input_alignment = os.path.abspath(input_alignment)
    recordAligns = get_record_alignments(input_alignment)
    cart_model = adjust_cart_model(input_fasta, input_alignment)

    if output_fasta is None:
        base, ext = os.path.splitext(input_fasta)
        output_fasta = f"{base}_clustered{ext}"
    else:
        output_fasta = os.path.abspath(output_fasta)

    with open(input_fasta) as infasta, open(output_fasta, "w") as outfasta:
        for record in SeqIO.parse(infasta, "fasta"):
            record_align_seq = recordAligns[record.id.split("_")[0]]
            cluster_id = get_nifh_cluster_id(record_align_seq, cart_model)
            record.id = f"{record.id}_cluster_{cluster_id}"
            record.name = ""
            record.description = ""
            SeqIO.write(record, outfasta, "fasta")


def add_cluster_to_nifh_dict(
    input_fasta: str,
    input_alignment: str,
    input_dict: str,
    output_dict: str = None,
    out_clusters_file: str = None,
) -> None:
    """
    Add assigned cluster to nifH sequence alignments
    in reference ID dictionary
    """
    input_fasta = os.path.abspath(input_fasta)
    input_alignment = os.path.abspath(input_alignment)
    input_dict = os.path.abspath(input_dict)
    cart_model = adjust_cart_model(input_fasta, input_alignment)
    ref_dict = read_from_pickle_file(input_dict)

    if output_dict is None:
        output_dict = input_dict
    else:
        output_dict = os.path.abspath(output_dict)

    for record in SeqIO.parse(input_alignment, "fasta"):
        cluster_id = get_nifh_cluster_id(list(record.seq), cart_model)
        ref_dict[record.id] += f"_cluster_{cluster_id}"

    # Temporary fix to remove ref_dict entries not found in alignment
    ref_dict = {
        ref_id: label  # if 'cluster' in label[-15:] else label + "_cluster_IV"
        for ref_id, label in ref_dict.items()
        if "cluster" in label[-15:]
    }

    save_to_pickle_file(ref_dict, output_dict)
    if out_clusters_file is not None:
        lines = ["id\tcluster\n"]
        with open(out_clusters_file, "w") as file:
            for ref_id, label in ref_dict.items():
                cluster_id = "_".join(label.split("_")[-2:])
                cluster_id = "cluster_IV" if "cluster" not in cluster_id else cluster_id
                lines.append(f"{ref_id}\t{cluster_id}\n")
            file.writelines(lines)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Classify nifH peptide sequences in clusters based on CART model",
        epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2021",
    )

    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    parser._action_groups.append(optional)

    required.add_argument(
        "--seqs", dest="seqs", type=str, required=True, help="path to fasta file"
    )
    required.add_argument(
        "--aln", dest="aln", type=str, required=True, help="path to fasta file"
    )
    required.add_argument(
        "--indict",
        dest="indict",
        type=str,
        required=True,
        help="path to input dictionary with reference IDs and labels",
    )
    optional.add_argument(
        "--outdict",
        dest="outdict",
        type=str,
        help="path to output ID dictionary in pickle format",
    )
    optional.add_argument(
        "--out_clusters_file",
        dest="clusters_file",
        type=str,
        help="path to output tsv file containing defined clusters",
    )

    args = parser.parse_args()
    if args.outdict is None:
        base, ext = os.path.splitext(args.indict)
        outdict = os.path.abspath(base + "_clustered" + ext)
    else:
        outdict = os.path.abspath(args.outdict)

    print("* Classifying nifH sequences according to CART model")
    add_cluster_to_nifh_dict(
        args.seqs,
        args.aln,
        args.indict,
        output_dict=outdict,
        out_clusters_file=args.clusters_file,
    )
