#!/usr/bin/env python
# conda activate traits

from pathlib import Path
import pyfastx

from phyloplacement.utils import readFromPickleFile
from phyloplacement.preprocessing import reIndexFASTA
from phyloplacement.database import runCDHIT, runHMMbuild
from phyloplacement.alignment import runMuscle, convertFastaAlnToPhylip, splitReferenceFromQueryAlignments, runTrimal
from phyloplacement.phylotree import runIqTree, runPapara, runHMMalign, runEPAng


# data_dir = Path('/home/robaina/Documents/MAR_database/')
# input_fasta = data_dir / 'mardb_proteins_V6_no_duplicates.fasta'
# Do Hmmer search and filtering

nxr_data = Path('/home/robaina/Documents/TRAITS/data/nxr/test_data/')
# nxr_fasta = nxr_data / 'mardb_proteins_V6_TIGR015180_1.fasta'
# nxr_fasta_reduced = nxr_data / "data_reduced.fasta"

# Preprocessing
"""
1) Assert correct sequence format for downstream analysis
2) Reduce redundancy: remove duplicates, get representatives with cd-hit
3) Relabel entries with temporary ids to avoid donwstream conflicts
"""

# 1) Assert  correct format


# 2) Reduce redundancy of database
# runCDHIT(
#     input_fasta=nxr_fasta,
#     output_fasta=str(nxr_data / 'data_reduced.fasta'),
#     additional_args=None
#     )

# print(f'Original database size: {len(pyfastx.Fasta(str(nxr_fasta)))}')
# print(f'Reduced database size: {len(pyfastx.Fasta(str(nxr_fasta_reduced)))}')

# Assign numbers to reference sequence labels for data processing
# reIndexFASTA(
#     input_fasta=str(nxr_fasta_reduced),
#     output_dir=str(nxr_data),
#     prefix='ref_'
#     )

# MSA on reduced database
# runMuscle(
#     input_fasta=str(nxr_data / 'data_reduced_short_ids_modified.fasta'),
#     output_file=None
# )

# Trimal
runTrimal(
    input_aln=str(nxr_data / 'data_reduced_short_ids_modified.fasta.aln'),
    output_aln=str(nxr_data / 'data_reduced_short_ids_modified.fasta.aln')
)

convertFastaAlnToPhylip(
    input_fasta_aln=str(nxr_data / 'data_reduced_short_ids_modified.fasta.aln'),
    output_file=str(nxr_data / 'data_reduced_short_ids_modified.phylip')
)

# Make tree
runIqTree(
    input_algns=str(nxr_data / 'data_reduced_short_ids_modified.fasta.aln'),
    output_dir=str(nxr_data),
    output_prefix=None,
    keep_recovery_files=False,
    substitution_model='TEST',
    additional_args=None
)

# Align query sequences with Papara 
# runPapara(
#     tree_nwk=str(nxr_data / 'data_reduced_short_ids.fasta.fasta.aln.contree'),
#     msa_phy=str(nxr_data / 'data_reduced_short_ids.fasta.fasta.aln.phylip'),
#     query_fasta=str(nxr_data / 'Nxr_kitzinger_2021_short_ids_modified.fasta'),
#     output_file=None,
#     additional_args=None
# )

# runPapara(
#     tree_nwk='/home/robaina/Documents/TRAITS/data/papara_test/alignment.phylip.contree',
#     msa_phy='/home/robaina/Documents/TRAITS/data/papara_test/alignment.phylip',
#     query_fasta='/home/robaina/Documents/TRAITS/data/papara_test/sequencesLongLabels.fasta',
#     output_file=None,
#     n_threads=None,
#     additional_args=None
# )

# Build HMM profile out of reference MSA
# runHMMbuild(
#     input_aln='/home/robaina/Documents/TRAITS/data/nxr/data_reduced_short_ids.fasta.fasta.aln',
#     output_hmm=None,
#     additional_args=None
#     )

# Align query sequences with hmmalign and split resulting fasta.aln
# runHMMalign(
#     input_hmm='/home/robaina/Documents/TRAITS/data/nxr/data_reduced_short_ids.fasta.fasta.aln.hmm',
#     input_aln='/home/robaina/Documents/TRAITS/data/nxr/data_reduced_short_ids.fasta.fasta.aln',
#     input_seqs='/home/robaina/Documents/TRAITS/data/nxr/kitzinger2021/Nxr_kitzinger_2021.fasta',
#     output_aln_seqs='/home/robaina/Documents/TRAITS/data/nxr/Nxr_kitzinger_2021_ref_aln.fasta.aln',
#     additional_args='--outformat afa'
#     )

# ref_dict = readFromPickleFile('/home/robaina/Documents/TRAITS/data/nxr/data_reduced_id_dict.pickle')
# ref_ids = set([str(k) for k in ref_dict.keys()])

# splitReferenceFromQueryAlignments(
#     ref_query_msa='/home/robaina/Documents/TRAITS/data/nxr/Nxr_kitzinger_2021_ref_aln.fasta.aln',
#     ref_ids=ref_ids,
#     out_dir=None
# )

# Run EPA-ng for placement of short reads
# runEPAng(
#     input_tree='/home/robaina/Documents/TRAITS/data/nxr/data_reduced_short_ids.fasta.fasta.aln.treefile',
#     input_aln_ref='/home/robaina/Documents/TRAITS/data/nxr/Nxr_kitzinger_2021_ref_aln.fasta.aln_ref_fraction.aln',
#     input_aln_query='/home/robaina/Documents/TRAITS/data/nxr/Nxr_kitzinger_2021_ref_aln.fasta.aln_query_fraction.aln',
#     model='LG+I+G4',
#     output_dir=None,
#     additional_args='--redo'
# )
