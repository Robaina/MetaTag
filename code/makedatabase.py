#!/usr/bin/env python
# conda activate traits
import os
from pathlib import Path

from phyloplacement.utils import readFromPickleFile, terminalExecute, countRecords
from phyloplacement.database.preprocessing import relabelRecordsInFASTA, assertCorrectSequenceFormat
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
# assertCorrectSequenceFormat(
#     fasta_file=str(work_dir / 'ref_reduced.fasta'),
#     output_file=str(work_dir / 'ref_reduced_clean.fasta'),
# )

# # 3) Assign numbers to reference sequence labels for data processing
# relabelRecordsInFASTA(
#     input_fasta=str(work_dir / 'ref_reduced_clean.fasta'),
#     output_dir=str(work_dir),
#     prefix='ref_'
#     )
