#!/usr/bin/env python
# conda activate traits
import os
from Bio import SeqIO
from pathlib import Path

from phyloplacement.utils import readFromPickleFile, terminalExecute
from phyloplacement.preprocessing import relabelRecordsInFASTA, reformatSequencesInFASTA
from phyloplacement.database import runCDHIT, runHMMbuild, filterFASTAByHMM
from phyloplacement.alignment import runMuscle, convertFastaAlnToPhylip, splitReferenceFromQueryAlignments, runTrimal
from phyloplacement.phylotree import runFastTree, runIqTree, runPapara, runHMMalign, runEPAng


def countRecords(fasta_file: str) -> None:
    cmd_str = f'grep -c ">" {fasta_file}'
    terminalExecute(cmd_str)

# data_dir = Path('/home/robaina/Documents/MAR_database/')
# input_fasta = data_dir / 'mardb_proteins_V6_no_duplicates.fasta'
# Do Hmmer search and filtering

test_data = Path('/home/robaina/Documents/TRAITS/tests/')
# nxr_fasta = test_data / 'mardb_proteins_V6_TIGR015180_1.fasta'
# nxr_fasta_reduced = test_data / "data_reduced.fasta"


# Make peptide-specific database
# filterFASTAByHMM(
#     hmm_model='/home/robaina/Documents/tigrfams/hmm_PGAP/TIGR00639.1.HMM',
#     input_fasta='/home/robaina/Documents/MAR_database/mardb_proteins_V6_no_duplicates.fasta',
#     output_fasta='/home/robaina/Documents/TRAITS/tests/mardb_TIGR00639.1.fasta'
# )


# Preprocessing
"""
1) Assert correct sequence format for downstream analysis
2) Reduce redundancy: remove duplicates, get representatives with cd-hit
3) Relabel entries with temporary ids to avoid donwstream conflicts
"""

# 1) Reduce redundancy of database
# runCDHIT(
#     input_fasta=str(test_data / 'mardb_TIGR00639.1.fasta'),
#     output_fasta=str(test_data / 'ref_reduced.fasta'),
#     additional_args='-c 0.6 -n 4'
#     )

# countRecords(str(test_data / "ref_reduced.fasta"))
# print(f'Original database size: {len(pyfastx.Fasta(str(test_data / "mardb_TIGR00639.1.fasta")))}')
# print(f'Reduced database size: {len(pyfastx.Fasta(str(test_data / "ref_reduced.fasta")))}')

# # 2) Assert  correct format
# reformatSequencesInFASTA(
#     fasta_file=str(test_data / 'ref_reduced.fasta'),
#     output_file=str(test_data / 'ref_reduced_clean.fasta'),
# )

# # Assign numbers to reference sequence labels for data processing
# relabelRecordsInFASTA(
#     input_fasta=str(test_data / 'ref_reduced_clean.fasta'),
#     output_dir=str(test_data),
#     prefix='ref_'
#     )

# MSA on reduced database
# runMuscle(
#     input_fasta=str(test_data / 'ref_reduced_clean_short_ids.fasta'),
#     output_file=str(test_data / 'ref_alignment.fasta.aln')
#     # additional_args='-physout alignment.phylip -phyiout alignment.iphylip'
# )

# # # Trimal
# # # runTrimal(
# # #     input_aln=str(test_data / 'ref_reduced_clean_short_ids.fasta.aln'),
# # #     output_aln=str(test_data / 'ref_alignment.fasta.aln')
# # # )

# convertFastaAlnToPhylip(
#     input_fasta_aln=str(test_data / 'ref_alignment.fasta.aln'),
#     output_file=str(test_data / 'ref_alignment.phylip')
# )

# # # # Make tree
# # # # runIqTree(
# # # #     input_algns=str(test_data / 'ref_alignment.phylip'),
# # # #     output_dir=str(test_data),
# # # #     output_prefix='ref_alignment',
# # # #     keep_recovery_files=True,
# # # #     substitution_model='TEST',
# # # #     additional_args=None
# # # # )

runFastTree(
    input_algns=str(test_data / 'ref_alignment.phylip'),
    output_file=str(test_data / 'ref_alignment.fasttree'),
    additional_args=None
)

# Align query sequences with Papara 
# runPapara(
#     tree_nwk=str(test_data / 'ref_alignment.fasttree'),
#     msa_phy='/home/robaina/Documents/TRAITS/alignment.phylip',#str(test_data / 'ref_alignment.phylip'),
#     query_fasta='/home/robaina/Documents/TRAITS/data/nxr/test_data/Nxr_kitzinger_2021_short_ids_modified.fasta',#'/home/robaina/Documents/TRAITS/data/papara_test/sequencesLongLabels.fasta',#str(test_data / 'Nxr_kitzinger_2021_short_ids_modified.fasta'),
#     output_file=None,
#     additional_args=None
# )
# os.remove('papara_log.aln')

# # Build HMM profile out of reference MSA
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
