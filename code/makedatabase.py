#!/usr/bin/env python
# conda activate traits
import os
from pathlib import Path

from phyloplacement.utils import readFromPickleFile, terminalExecute, countRecords
from phyloplacement.preprocessing import relabelRecordsInFASTA, reformatSequencesInFASTA
from phyloplacement.database import runCDHIT, runHMMbuild, filterFASTAByHMM
from phyloplacement.alignment import runMuscle, convertFastaAlnToPhylip, convertPhylipToFastaAln,splitReferenceFromQueryAlignments, runTrimal
from phyloplacement.phylotree import runFastTree, runIqTree, runPapara, runHMMalign, runEPAng, runTreeShrink


"""
TIGR00639.1.HMM (large peptide database after cd-hit)
TIGR01580.1.HMM (narG / nxr)
"""
mardb_data = Path('/home/robaina/Documents/MAR_database/')
tigr_data = Path('/home/robaina/Documents/tigrfams/hmm_PGAP/')
work_dir = Path('/home/robaina/Documents/TRAITS/tests/')

# # Make peptide-specific database
# filterFASTAByHMM(
#     hmm_model=str(tigr_data / 'TIGR01580.1.HMM'),
#     input_fasta=str(mardb_data / 'mardb_proteins_V6_no_duplicates.fasta'),
#     output_fasta=str(work_dir / 'mardb_TIGR01580.1.fasta')
# )

# # Preprocessing
# """
# 1) Assert correct sequence format for downstream analysis
# 2) Reduce redundancy: remove duplicates, get representatives with cd-hit
# 3) Relabel entries with temporary ids to avoid donwstream conflicts
# """

# # 1) Reduce redundancy of reference database
# runCDHIT(
#     input_fasta=str(work_dir / 'mardb_TIGR01580.1.fasta'),
#     output_fasta=str(work_dir / 'ref_reduced.fasta'),
#     additional_args=None
#     )
# countRecords(str(work_dir / "ref_reduced.fasta"))


# # 2) Assert  correct format
# reformatSequencesInFASTA(
#     fasta_file=str(work_dir / 'ref_reduced.fasta'),
#     output_file=str(work_dir / 'ref_reduced_clean.fasta'),
# )

# # 3) Assign numbers to reference sequence labels for data processing
# relabelRecordsInFASTA(
#     input_fasta=str(work_dir / 'ref_reduced_clean.fasta'),
#     output_dir=str(work_dir),
#     prefix='ref_'
#     )

# # MSA on reduced database
# runMuscle(
#     input_fasta=str(work_dir / 'ref_reduced_clean_short_ids.fasta'),
#     output_file=str(work_dir / 'ref_alignment.fasta.aln')
# )

# # # Trimal
# # # runTrimal(
# # #     input_aln=str(test_data / 'ref_reduced_clean_short_ids.fasta.aln'),
# # #     output_aln=str(test_data / 'ref_alignment.fasta.aln')
# # # )

# convertFastaAlnToPhylip(
#     input_fasta_aln=str(work_dir / 'ref_alignment.fasta.aln'),
#     output_file=str(work_dir / 'ref_alignment.phylip')
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

# runFastTree(
#     input_algns=str(work_dir / 'ref_alignment.phylip'),
#     output_file=str(work_dir / 'ref_alignment.fasttree'),
#     additional_args=None
# )

# Trim tree to remove outlier branches
runTreeShrink(
    input_tree='',
    input_aln='',
    output_dir='',
    output_deleted_nodes=True,
    additional_args=None
)

# # Preprocess query sequences
# reformatSequencesInFASTA(
#     fasta_file='/home/robaina/Documents/TRAITS/data/nxr/kitzinger2021/nxr_kitzinger_2021.fasta',
#     output_file=str(work_dir / 'query_sequences.fasta')
# )

# relabelRecordsInFASTA(
#     input_fasta=str(work_dir / 'query_sequences.fasta'),
#     output_dir=str(work_dir),
#     prefix='query_'
# )

# # Align query sequences with Papara 
# runPapara(
#     tree_nwk=str(work_dir / 'ref_alignment.fasttree'),
#     msa_phy=str(work_dir / 'ref_alignment.phylip'),
#     query_fasta=str(work_dir / 'query_sequences_short_ids.fasta'),
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

# # Split reference from query msa for EPA-ng 
# ref_ids = readFromPickleFile(
#     '/home/robaina/Documents/TRAITS/tests/ref_reduced_clean_id_dict.pickle'
#     )

# convertPhylipToFastaAln(
#     input_phylip=pfile,
#     output_file=None
# )

# splitReferenceFromQueryAlignments(
#     ref_query_msa='/home/robaina/Documents/TRAITS/tests/papara_alignment.fasta.aln',
#     ref_ids=ref_ids.keys(),
#     out_dir=None
# )

# Run EPA-ng for placement of short reads
runEPAng(
    input_tree=str(work_dir / 'ref_alignment.fasttree'),
    input_aln_ref=str(work_dir / 'papara_alignment.fasta_ref_fraction.aln'),
    input_aln_query=str(work_dir/ 'papara_alignment.fasta_query_fraction.aln'),
    model='JTT', #'LG+I+G4',
    output_dir=str(work_dir),
    additional_args='--redo'
)
