#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to create peptide-specific sequence databases

1. Implement hmmr
2. Filter fasta files based on query sequences
"""

from __future__ import annotations
import os
import warnings
from collections import defaultdict

import pandas as pd
from Bio import SearchIO, SeqIO, AlignIO
import pyfastx

import phyloplacement.wrappers as wrappers
from phyloplacement.utils import set_default_output_path
from phyloplacement.database.labelparsers import MARdbLabelParser


def filter_fasta_by_sequence_length(
    input_fasta: str,
    min_length: int = None,
    max_length: int = None,
    output_fasta: str = None,
) -> None:
    """
    Filter sequences by length in fasta file
    """
    if (min_length is None) and (max_length is None):
        warnings.warn("Missing boundary values for sequence length")
        return
    input_fasta = os.path.abspath(input_fasta)
    fa = pyfastx.Fasta(input_fasta)
    record_ids = fa.keys()
    if min_length is None:
        min_length = 0
    if max_length is not None:
        max_tag = str(max_length)
        record_ids.filter(record_ids >= min_length, record_ids <= max_length)
    else:
        max_tag = ""
        record_ids.filter(record_ids >= min_length)
    if output_fasta is None:
        output_fasta = set_default_output_path(
            input_fasta, f"_length_{min_length}_{max_tag}"
        )
    if not record_ids:
        raise ValueError("No records found with given sequence length bounds")
    with open(output_fasta, "w") as fp:
        for record_id in record_ids:
            record_obj = fa[record_id]
            fp.write(record_obj.raw)
    os.remove(input_fasta + ".fxi")


def parse_hmmsearch_output(hmmer_output: str) -> pd.DataFrame:
    """
    Parse hmmsearch or hmmscan summary table output file
    """
    attribs = ["id", "bias", "bitscore", "description"]
    hits = defaultdict(list)
    with open(hmmer_output) as handle:
        for queryresult in SearchIO.parse(handle, "hmmer3-tab"):
            for hit in queryresult.hits:
                for attrib in attribs:
                    hits[attrib].append(getattr(hit, attrib))
    return pd.DataFrame.from_dict(hits)


def filter_fasta_by_ids(
    input_fasta: str, record_ids: list, output_fasta: str = None
) -> None:
    """
    Filter records in fasta file matching provided IDs
    """
    if output_fasta is None:
        output_fasta = set_default_output_path(input_fasta, "_fitered")
    if os.path.exists(input_fasta + ".fxi"):
        os.remove(input_fasta + ".fxi")
    record_ids = set(record_ids)
    fa = pyfastx.Fasta(input_fasta)
    with open(output_fasta, "w") as fp:
        for record_id in record_ids:
            try:
                record_obj = fa[record_id]
                fp.write(record_obj.raw)
            except:
                pass
    os.remove(input_fasta + ".fxi")


def filter_fasta_by_hmm(
    hmm_model: str,
    input_fasta: str,
    output_fasta: str = None,
    hmmer_output: str = None,
    method: str = "hmmsearch",
    additional_args: str = None,
) -> None:
    """
    Generate protein-specific database by filtering
    sequence database to only contain sequences
    corresponing to protein of interest

    @Arguments:
    additional_args: additional arguments to hmmsearch or hmmscan
    """
    hmm_name, _ = os.path.splitext(os.path.basename(hmm_model))
    if hmmer_output is None:
        hmmer_output = set_default_output_path(
            input_fasta, tag=f"_{hmm_name}", extension=".txt"
        )
    if output_fasta is None:
        output_fasta = set_default_output_path(input_fasta, tag=f"filtered_{hmm_name}")

    print("Running Hmmer...")
    wrappers.run_hmmsearch(
        hmm_model=hmm_model,
        input_fasta=input_fasta,
        output_file=hmmer_output,
        method=method,
        additional_args=additional_args,
    )
    print("Parsing Hmmer output file...")
    hmmer_hits = parse_hmmsearch_output(hmmer_output)
    if not hmmer_hits.id.values.tolist():
        raise ValueError("No records found in database matching provided hmm")
    print("Filtering Fasta...")
    filter_fasta_by_ids(
        input_fasta, record_ids=hmmer_hits.id.values, output_fasta=output_fasta
    )


def convert_fasta_aln_to_phylip(
    input_fasta_aln: str, output_phylip: str = None
) -> None:
    """
    Convert alignments in Fasta to Phylip.
    """
    if output_phylip is None:
        output_phylip = set_default_output_path(input_fasta_aln, extension=".phylip")
    with open(input_fasta_aln, "r") as input_handle, open(
        output_phylip, "w"
    ) as output_handle:
        alignments = AlignIO.parse(input_handle, "fasta")
        AlignIO.write(alignments, output_handle, "phylip-relaxed")


def convert_phylip_to_fasta_aln(input_phylip: str, output_file: str = None) -> None:
    """
    Convert alignments in Phylip to Fasta format
    """
    if output_file is None:
        output_file = set_default_output_path(input_phylip, extension=".faln")
    alignments = AlignIO.parse(input_phylip, "phylip-relaxed")
    AlignIO.write(alignments, output_file, "fasta")


def convert_stockholm_to_fasta_aln(
    input_stockholm: str, output_fasta: str = None
) -> None:
    """
    Convert alignment file in Stockholm format to fasta
    """
    if output_fasta is None:
        output_fasta = set_default_output_path(input_stockholm, extension=".faln")
    alignments = AlignIO.read(input_stockholm, "stockholm")
    AlignIO.write(alignments, output_fasta, "fasta")


def split_reference_from_query_alignments(
    ref_query_msa: str, ref_ids: set = None, ref_prefix: str = None, out_dir: str = None
) -> None:
    """
    Separate reference sequences from query sequences in msa fasta file
    """
    if out_dir is None:
        out_dir = os.path.dirname(ref_query_msa)
    if (ref_ids is not None) and (ref_prefix is None):

        def is_reference(record_name):
            return record_name in ref_ids

    elif (ref_ids is None) and (ref_prefix is not None):

        def is_reference(record_name):
            return record_name.lower().startswith(ref_prefix.lower())

    else:
        raise ValueError("Provide either set of ref ids or ref prefix")
    out_ref_msa = set_default_output_path(ref_query_msa, tag="_ref_fraction")
    out_query_msa = set_default_output_path(ref_query_msa, tag="_query_fraction")

    fasta = pyfastx.Fasta(ref_query_msa, build_index=False, full_name=True)
    with open(out_ref_msa, "w") as outref, open(out_query_msa, "w") as outquery:
        for record_name, record_seq in fasta:
            if is_reference(record_name):
                outref.write(f">{record_name}\n{record_seq}\n")
            else:
                outquery.write(f">{record_name}\n{record_seq}\n")


def get_fasta_record_ids(fasta_file: str) -> set:
    """
    Extract record ids from fasta
    """
    fasta = pyfastx.Fasta(fasta_file, full_name=True)
    return set(fasta.keys())


def remove_records_from_fasta(
    input_fasta: str, record_ids: list[str], output_fasta: str = None
) -> None:
    """
    Remove records from Fasta file based on provided list of record IDs.
    """
    if output_fasta is None:
        output_fasta = set_default_output_path(input_fasta, tag="_reduced")
    fasta_record_ids = get_fasta_record_ids(input_fasta)
    records_to_keep = fasta_record_ids - set(record_ids)
    filter_fasta_by_ids(input_fasta, records_to_keep, output_fasta)
    return None


def count_records_in_fasta(fasta_file: str) -> int:
    """
    Count records in fasta file
    """
    with open(fasta_file, "r") as file:
        n_records = sum([1 for line in file if line.startswith(">")])
    return n_records


def slice_fasta(input_file, output_file, N):
    n = 0
    records = SeqIO.parse(input_file, "fasta")
    sliced_records = []
    for record in records:
        if n < N:
            sliced_records.append(record)
        else:
            break
        n += 1
    SeqIO.write(sliced_records, output_file, "fasta")


def is_empty_fasta(fasta_file):
    with open(fasta_file, "r") as file:
        return ">" not in file.read()


class LinkedHMMfilter:
    """
    Tools to search for linkage structures among sets of hmm models
    """

    def __init__(self, hmm_hits: dict) -> None:
        """
        Search for contigs that satisfy the given gene linkage structure

        @param: hmm_hits, a dict of pandas DataFrames, as output by
                parse_hmmsearch_output with keys corresponding to hmm names
        """
        self._hmm_hits = hmm_hits

    @staticmethod
    def parse_linkage_string(link_str: str) -> list[list]:
        """
        Parse linkage structure string. A linkage structure
        is a string like the following:

        >hmm_a n_ab <hmm_b n_bc hmm_c

        where '>' indicates a hmm target located on the positive strand,
        '<' a target located on the negative strand, and n_ab cooresponds
        to the maximum number of genes separating matched gene a and b.
        Multiple hmms may be employed (limited by computational capabilities).
        No order symbol in a hmm indicates that results should be independent
        of strand location.
        """

        def split_strand_from_locus(locus_str: str) -> tuple:
            if locus_str[0] == "<" or locus_str[0] == ">":
                sense = locus_str[0]
                locus_str = locus_str[1:]
                strand = "pos" if sense == ">" else "neg"
            else:
                strand = None
            return (strand, locus_str)

        links = link_str.split()
        max_dists = [int(dist) for dist in links if dist.isdigit()]
        hmms = [h for h in links if not h.isdigit()]
        pairs = list(zip(hmms, hmms[1:]))
        parsed_struc = []
        for pair, dist in zip(pairs, max_dists):
            sense_a, locus_a = split_strand_from_locus(pair[0])
            sense_b, locus_b = split_strand_from_locus(pair[1])
            parsed_struc.append([locus_a, locus_b, dist, sense_a, sense_b])
        return parsed_struc

    @staticmethod
    def filter_hits_by_linkage_pair(
        hit_labels_a: pd.DataFrame,
        hit_labels_b: pd.DataFrame,
        dist_ab: int,
        strand_a: str = None,
        strand_b: str = None,
    ) -> pd.DataFrame:
        """
        Get labels compatible with linked pair
        """

        def in_right_strand(hit_strand: str, strand: str):
            if strand is not None:
                return hit_strand == strand
            else:
                return True

        def get_linked_hit_labels_with_strand():
            linked_labels = set()
            linked_hit_labels = pd.DataFrame(columns=hit_labels_a.columns)

            for i, hit_a in hit_labels_a.iterrows():
                if (
                    in_right_strand(hit_a.strand, strand_a)
                    and hit_a.full not in linked_labels
                ):
                    contig_a = hit_a.contig
                    pos_a = hit_a.gene_pos
                    linked_b_hits = hit_labels_b[
                        (
                            (hit_labels_b.contig == contig_a)
                            & (abs(hit_labels_b.gene_pos - pos_a) <= dist_ab)
                            & (hit_labels_b.strand == strand_b)
                        )
                    ]
                    linked_b_labels = linked_b_hits.full.values.tolist()

                    if linked_b_labels:
                        linked_labels.update([hit_a.full] + linked_b_labels)
                        linked_hit_labels = linked_hit_labels.append(
                            linked_b_hits
                        ).append(hit_a)

            return linked_hit_labels.drop_duplicates()

        def get_linked_hit_labels():
            linked_labels = set()
            linked_hit_labels = pd.DataFrame(columns=hit_labels_a.columns)

            for i, hit_a in hit_labels_a.iterrows():
                if hit_a.full not in linked_labels:
                    contig_a = hit_a.contig
                    pos_a = hit_a.gene_pos
                    linked_b_hits = hit_labels_b[
                        (
                            (hit_labels_b.contig == contig_a)
                            & (abs(hit_labels_b.gene_pos - pos_a) <= dist_ab)
                        )
                    ]
                    linked_b_labels = linked_b_hits.full.values.tolist()

                    if linked_b_labels:
                        linked_labels.update([hit_a.full] + linked_b_labels)
                        linked_hit_labels = linked_hit_labels.append(
                            linked_b_hits
                        ).append(hit_a)

            return linked_hit_labels.drop_duplicates()

        if strand_a is None and strand_b is None:
            return get_linked_hit_labels()
        else:
            return get_linked_hit_labels_with_strand()

    def filter_hits_by_linked_hmm_structure(self, link_structure: str) -> pd.DataFrame:
        """
        Search for contigs that satisfy the given gene linkage structure
        @param: link_structure, a str describing the desired linkage structure,
                structured as follows:

                'hmm_a N_ab hmm_b'

                where N_ab corresponds to the maximum number of genes separating
                gene found by hmm_a and gene found by hmm_b, and hmm_ corresponds
                to the name of the hmm as provided in the keys of hmm_hits.
                More than two hmms can be concatenated.
        """
        parsed_struc = self.parse_linkage_string(link_structure)

        hit_labels = {}
        mardblabel = MARdbLabelParser()
        for hmm, hits in self._hmm_hits.items():
            labels = hits.id.values.tolist()
            if not labels:
                raise ValueError(f"No records found in database matching HMM: {hmm}")
            hit_labels[hmm] = mardblabel.parse_from_list(labels)

        for n, linked_pair in enumerate(parsed_struc):

            if n < 1:
                locus_a, locus_b, dist_ab, strand_a, strand_b = linked_pair
                hit_labels_a, hit_labels_b = hit_labels.pop(locus_a), hit_labels.pop(
                    locus_b
                )
                linked_hit_labels = self.filter_hits_by_linkage_pair(
                    hit_labels_a, hit_labels_b, dist_ab, strand_a, strand_b
                )
            else:
                locus_a = f"{locus_a}_{locus_b}"
                new_locus_b, dist_ab, strand_a, strand_b = linked_pair[1:]
                hit_labels[locus_a] = linked_hit_labels
                hit_labels_a, hit_labels_b = hit_labels.pop(locus_a), hit_labels.pop(
                    new_locus_b
                )
                linked_hit_labels = self.filter_hits_by_linkage_pair(
                    hit_labels_a, hit_labels_b, dist_ab, strand_a, strand_b
                )

        return linked_hit_labels

    def partition_linked_labels_by_hmm(
        self, linked_hit_labels: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Partition linked labels dataframe into several dataframes containing
        hmm-specific hits.  This is a workaround to avoid extensive modification of
        LinkedHMMfilter.filter_hits_by_linked_hmm_structure.
        """
        return {
            hmm_name: linked_hit_labels[linked_hit_labels.full.isin(hits.id)]
            for hmm_name, hits in self._hmm_hits.items()
        }


def filter_fasta_by_hmm_structure(
    hmm_structure: str,
    target_hmm: str,
    input_fasta: str,
    input_hmms: list[str],
    output_fasta: str = None,
    output_dir: str = None,
    hmmer_output_dir: str = None,
    reuse_hmmer_results: bool = True,
    method: str = "hmmsearch",
    additional_args: list[str] = None,
) -> None:
    """
    Generate protein-specific database by filtering sequence database
    to only contain sequences which satisfy the provided (gene/hmm)
    structure

    @Arguments:
    additional_args: additional arguments to hmmsearch or hmmscan. Each
    element in the list is a string with additional arguments for each
    input hmm (arranged in the same order), an element can also take a
    value of None to avoid passing additional arguments for a specific
    input hmm. A single string may also be passed, in which case the
    same additional argument is passed to hmmsearch for all input hmms
    """
    if output_fasta is None:
        output_fasta = set_default_output_path(
            input_fasta, tag=f'filtered_{hmm_structure.replace(" ", "_")}'
        )
    if hmmer_output_dir is None:
        hmmer_output_dir = os.path.join(
            set_default_output_path(input_fasta, only_dirname=True), "hmmer_outputs"
        )

    if additional_args is None:
        additional_args = [None for _ in input_hmms]

    if type(additional_args) == str:
        additional_args = [additional_args for _ in input_hmms]

    elif type(additional_args) == list:
        if len(additional_args) == 1:
            additional_args = [additional_args[0] for _ in input_hmms]

        if (len(additional_args) > 1) and (len(additional_args) < len(input_hmms)):
            raise ValueError(
                "Provided additional argument strings are less than the number of input hmms."
            )
    else:
        raise ValueError(
            "Additional arguments must be: 1) a list[str], 2) a str, or 3) None"
        )
    if not os.path.isdir(hmmer_output_dir):
        os.mkdir(hmmer_output_dir)

    if output_dir is None:
        output_dir = set_default_output_path(input_fasta, only_dirname=True)
    else:
        output_dir = output_dir

    print("Running Hmmer...")
    hmm_hits = {}
    for hmm_model, add_args in zip(input_hmms, additional_args):
        hmm_name, _ = os.path.splitext(os.path.basename(hmm_model))
        hmmer_output = os.path.join(hmmer_output_dir, f"hmmer_output_{hmm_name}.txt")

        if not (reuse_hmmer_results and os.path.isfile(hmmer_output)):
            wrappers.run_hmmsearch(
                hmm_model=hmm_model,
                input_fasta=input_fasta,
                output_file=hmmer_output,
                method=method,
                additional_args=add_args,
            )

        hmm_hits[hmm_name] = parse_hmmsearch_output(hmmer_output)

    print("Filtering results by HMM structure...")
    linkedfilter = LinkedHMMfilter(hmm_hits)
    linked_hit_labels = linkedfilter.filter_hits_by_linked_hmm_structure(hmm_structure)

    if not linked_hit_labels.full.values.tolist():
        raise ValueError("No records found in database matching provided hmm structure")

    print("Filtering Fasta...")
    partitioned_hit_labels = linkedfilter.partition_linked_labels_by_hmm(
        linked_hit_labels
    )
    for hmm_name in hmm_hits:
        if hmm_name in target_hmm:
            outfasta = output_fasta
        else:
            outfasta = os.path.join(output_dir, f"{hmm_name}_hits.fasta")
        record_ids = partitioned_hit_labels[hmm_name].full.values
        filter_fasta_by_ids(input_fasta, record_ids=record_ids, output_fasta=outfasta)
