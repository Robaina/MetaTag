"""
Tools to reduce the size of the peptide-specific
reference database

Currently based on:
1. CD-HIT
2. Repset: https://onlinelibrary.wiley.com/doi/10.1002/prot.25461
"""

import os
import shutil
import warnings

import phyloplacement.wrappers as wrappers
from phyloplacement.utils import setDefaultOutputPath
from phyloplacement.utils import TemporaryFilePath
from phyloplacement.database.preprocessing import getRepresentativeSet


def reduceDatabaseRedundancy(input_fasta: str,
                             output_fasta: str = None,
                             cdhit: bool = True, maxsize: int = None,
                             cdhit_args: str = None) -> None:
    """
    Reduce redundancy of peptide datatabase.
    Runs cd-hit, if selected, additional arguments to cdhit
    may be passed as a string (cdhit_args).
    Runs repset to obtain a final database size no larger 
    (number of sequences) than selected maxsize.
    If maxsize = None, repset is not run.
    """
    if (not cdhit) and (maxsize is None):
        warnings.warn("No reduction algorithm has been selected.")

    if output_fasta is None:
        output_fasta = setDefaultOutputPath(input_fasta,  tag="_reduced")

    with TemporaryFilePath() as tempaln,\
         TemporaryFilePath() as tempfasta,\
         TemporaryFilePath() as tempfasta2,\
         TemporaryFilePath() as tempident:

        if cdhit:
            wrappers.runCDHIT(
                input_fasta=input_fasta,
                output_fasta=tempfasta,
                additional_args=cdhit_args
                )
            os.remove(tempfasta + ".clstr")
        else:
            shutil.move(input_fasta, tempfasta)
    
        if maxsize is not None:
            wrappers.runMAFFT(
                input_fasta=tempfasta,
                output_file=tempaln,
                n_threads=-1,
                parallel=True,
                additional_args='--retree 1 --maxiterate 0'
            )
            
            wrappers.getPercentIdentityFromMSA(
                input_msa=tempaln,
                output_file=tempident
            )

            print('Finding representative sequences for reference database...')
            getRepresentativeSet(
                input_seqs=tempfasta,
                input_PI=tempident,
                max_size=maxsize,
                outfile=tempfasta2
            )
            shutil.move(tempfasta2, tempfasta)

        shutil.move(tempfasta, output_fasta)