#!/usr/bin/env python
# conda activate traits

import re
import pyfastx
from phyloplacement.utils import readFromPickleFile, setDefaultOutputPath, parallelizeOverInputFiles
from phyloplacement.database.mardb import getMarDBentryCode, getMARdbGenomeByEntryCode


if __name__ == '__main__':

    # Retrieve mardb genomes
    label_dict = readFromPickleFile('/home/robaina/Documents/TRAITS/tests/ref_reduced_clean_id_dict.pickle')
    nxr_entry_codes = {getMarDBentryCode(v) for v in label_dict.values()}

    # for entry_code in nxr_entry_codes:
    #     print(f'Running entry code: {entry_code}')
    #     getMARdbGenomeByEntryCode(
    #         input_fasta='/home/robaina/Documents/MAR_database/mardb_assembly_V6.fa',
    #         entry_code=entry_code,
    #         output_fasta=f'/home/robaina/Documents/TRAITS/tests/nxr_genomes/{entry_code}.fa'
    #     )
    
    # parallelizeOverInputFiles(
    #     getMARdbGenomeByEntryCode,
    #     input_list=list(nxr_entry_codes)[:10],
    #     input_fasta='/home/robaina/Documents/MAR_database/mardb_assembly_V6.fa',
    #     n_processes=3
    # )
    print('Done!')


# # MSA on reduced database
# runMuscle(
#     input_fasta=str(work_dir / 'ref_reduced_clean_short_ids.fasta'),
#     output_file=str(work_dir / 'ref_alignment.fasta.aln')
# )

# Trimal
# runTrimal(
#     input_aln=str(work_dir / 'ref_alignment.fasta.aln'),
#     output_aln=str(work_dir / 'ref_alignment_trimal.fasta.aln')
# )

# convertFastaAlnToPhylip(
#     input_fasta_aln=str(work_dir / 'ref_alignment_trimal.fasta.aln'),
#     output_file=str(work_dir / 'ref_alignment_trimal.phylip')
# )

# # Make tree
# runIqTree(
#     input_algns=str(work_dir / 'ref_alignment_trimal.phylip'),
#     output_dir=str(work_dir),
#     output_prefix='ref_alignment',
#     keep_recovery_files=True,
#     substitution_model='TEST',
#     additional_args=None
# )

# runFastTree(
#     input_algns=str(work_dir / 'ref_alignment.phylip'),
#     output_file=str(work_dir / 'ref_alignment.fasttree'),
#     additional_args=None
# )

# Trim tree to remove outlier branches
# runTreeShrink(
#     input_tree=str(work_dir / 'ref_alignment.contree'),
#     input_aln=str(work_dir / 'ref_alignment_trimal.fasta.aln'),
#     output_dir=str(work_dir),
#     output_deleted_nodes=True,
#     additional_args=None
# )

# convertFastaAlnToPhylip(
#     input_fasta_aln=str(work_dir / 'ref_alignment_trimal.fasta_shrink.aln'),
#     output_file=str(work_dir / 'ref_alignment_trimal_shrink.phylip')
# )

# # Preprocess query sequences
# assertCorrectSequenceFormat(
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
# runEPAng(
#     input_tree=str(work_dir / 'ref_alignment.fasttree'),
#     input_aln_ref=str(work_dir / 'papara_alignment.fasta_ref_fraction.aln'),
#     input_aln_query=str(work_dir/ 'papara_alignment.fasta_query_fraction.aln'),
#     model='JTT', #'LG+I+G4',
#     output_dir=str(work_dir),
#     additional_args='--redo'
# )