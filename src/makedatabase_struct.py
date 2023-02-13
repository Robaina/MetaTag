#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Reference database:
1) Run hmmer to extract peptides of interest following gene structure
2) Reduce redundancy: cd-hit and/or repset
3) Relabel entries with temporary ids to avoid donwstream conflicts
"""

import os
import shutil
import argparse

from phyloplacement.utils import set_default_output_path, TemporaryFilePath
from phyloplacement.database.preprocessing import set_temp_record_ids_in_fasta
from phyloplacement.database.manipulation import (
    filter_fasta_by_hmm_structure,
    filter_fasta_by_sequence_length,
)
from phyloplacement.database.reduction import reduce_database_redundancy


parser = argparse.ArgumentParser(
    description="Build peptide reference database from HMM structure",
    epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2021",
)

optional = parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
parser._action_groups.append(optional)

required.add_argument(
    "--hmms",
    dest="hmms",
    type=str,
    nargs="+",
    required=True,
    help="list of space-separated paths to tigrfam or pfam models",
)
required.add_argument(
    "--target_hmm",
    dest="target_hmm",
    type=str,
    required=True,
    help=(
        "name of target hmm model to be used to generate sequence database. "
        "Name must be equal to the name of one of the provided hmms"
    ),
)
required.add_argument(
    "--hmm_struc",
    dest="hmm_struc",
    type=str,
    required=True,
    help=(
        "string displaying hmm sctructure to search for, such as: \n"
        '">hmm_a n_ab <hmm_b n_bc hmm_c", \n'
        'where ">" indicates a hmm target located on the positive strand, '
        '"<" a target located on the negative strand, and n_ab cooresponds '
        "to the maximum number of genes separating matched gene a and b. \n"
        "Multiple hmms may be employed (limited by computational capabilities)."
        "No order symbol in a hmm indicates that results should be independent "
        "of strand location. "
        "The script outputs a main file containing sequences matching the provided "
        "hmm structure and corresponding to the main target indicated in the argument: "
        "target_hmm. These sequences are filtered by sequence length and by maximum "
        "number of total sequences. "
        "The script also outputs additional file containing the matched records for the "
        "other, non-target, hmms. However, these files are not processed any further."
    ),
)
required.add_argument(
    "--in", dest="data", type=str, required=True, help="path to peptide database"
)
optional.add_argument(
    "--outdir", dest="outdir", type=str, help="path to output directory"
)
optional.add_argument(
    "--prefix",
    dest="prefix",
    type=str,
    default="",
    help="prefix to be added to output files",
)
optional.add_argument(
    "--max_size",
    dest="maxsize",
    default=None,
    type=int,
    help=("maximum size of representative set of sequences. " "Defaults to full set."),
)
optional.add_argument(
    "--min_seq_length",
    dest="minseqlength",
    default=None,
    type=int,
    help=("minimum sequence length in reference database. " "Defaults to zero"),
)
optional.add_argument(
    "--hmmsearch_args",
    dest="hmmsearch_args",
    type=str,
    default=None,
    required=False,
    help=(
        "list of comma-separated additional arguments to hmmsearch for each input hmm. "
        "A single argument may be provided, in which case the same additional argument "
        "is employed in all hmms."
    ),
)
parser.add_argument(
    "--max_seq_length",
    dest="maxseqlength",
    default=None,
    type=int,
    required=False,
    help=("maximum sequence length in reference database. " "Defaults to inf"),
)
parser.add_argument(
    "--relabel",
    dest="relabel",
    action="store_true",
    required=False,
    default=False,
    help=(
        "relabel record IDs with numerical ids. "
        "Unrequired to build database, but highly recommended "
        "to avoid possible conflicts downstream the pipeline."
    ),
)


args = parser.parse_args()
if args.outdir is None:
    args.outdir = set_default_output_path(args.data, only_dirname=True)
hmmsearch_args = list(map(lambda x: x.strip(), args.hmmsearch_args.split(",")))
hmmsearch_args = list(map(lambda x: None if x == "None" else x, hmmsearch_args))
hmmer_output_dir = os.path.join(args.outdir, "hmmer_outputs/")
output_fasta = os.path.join(args.outdir, f"{args.prefix}ref_database.faa")
output_fasta_short = os.path.join(
    args.outdir, f"{args.prefix}ref_database_short_ids.faa"
)


def main():

    print("* Making peptide-specific reference database...")
    with TemporaryFilePath() as tempfasta, TemporaryFilePath() as tempfasta2:
        filter_fasta_by_hmm_structure(
            hmm_structure=args.hmm_struc,
            target_hmm=args.target_hmm,
            input_fasta=args.data,
            input_hmms=args.hmms,
            output_fasta=tempfasta,
            output_dir=args.outdir,
            hmmer_output_dir=hmmer_output_dir,
            reuse_hmmer_results=True,
            method="hmmsearch",
            additional_args=hmmsearch_args,  #'--cut_nc'
        )

        if (args.minseqlength is not None) or (args.maxseqlength is not None):
            print("* Filtering sequences by established length bounds...")
            filter_fasta_by_sequence_length(
                input_fasta=tempfasta,
                min_length=args.minseqlength,
                max_length=args.maxseqlength,
                output_fasta=tempfasta2,
            )
            shutil.move(tempfasta2, tempfasta)

        reduce_database_redundancy(
            input_fasta=tempfasta,
            output_fasta=output_fasta,
            cdhit=True,
            cdhit_args=None,
            maxsize=args.maxsize,
        )

    if args.relabel:
        print("* Relabelling records in reference database...")
        set_temp_record_ids_in_fasta(
            input_fasta=output_fasta,
            output_dir=args.outdir,
            prefix=f"ref_{args.prefix}",
        )
        shutil.move(output_fasta_short, output_fasta)
    print("Finished!")


if __name__ == "__main__":
    main()
