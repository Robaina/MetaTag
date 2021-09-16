"""
Functions to preprocess sequence data
"""

import os
from posixpath import basename
import pandas as pd
from collections import defaultdict
from Bio import SearchIO, SeqIO

from .utils import terminalExecute


# def removeDuplicates():
#     """
#     Remove sequences in fasta files with repeated labels
#     or sequences.

#     NOTES:

#     script: joinseqs.py
#     """
#     pass

# def reformatFileName():
#     """
#     Checks if file name format is legal
#     and corrects otherwise


#     NOTES:

#     script: clean.py

#     1. File names, must contain only upper/lower case letters and digits and '_',
#     replace anything else (such as a space) by '_'

#     2. Sequences must be only composed of uppercase letters A-Z
#     """
#     pass

# def refomatSequencesAndLabels():
#     """
#     Checks for inconsistencies in sequence data and
#     reformat sequence labels

#     NOTES:

#     script: clean.py

#     1. File names, must contain only upper/lower case letters and digits and '_',
#     replace anything else (such as a space) by '_'

#     2. Sequences must be only composed of uppercase letters A-Z
#     """
#     pass

# def reformatProdigalLabels():
#     """
#     Reformat prodigal labels and file names

#     NOTES:

#     script: reformatabel.py
#     """
#     pass

# def translateDNA():
#     """
#     Translate DNA sequences with prodigal

#     NOTES:

#     script: loopparallel.py
#     """
#     pass

def runHMMER(hmm_model: str, input_fasta: str,
             output_file: str = None,
             method: str = 'hmmsearch') -> None:
    """
    Simple CLI wrapper to run hmmsearch or hmmscan
    Requires hmmer installed and accessible
    """
    if output_file is None:
        basename, ext = os.path.splitext(input_fasta)
        output_file = f'{basename}_hmmer_hits.txt'
    cmd_str = (f'{method} --cut_ga --tblout {output_file} '
               f'{hmm_model} {input_fasta}')
    terminalExecute(cmd_str, suppress_output=True)

def parseHMMERoutput(hmmer_output: str) -> pd.DataFrame:
    """
    Parse hmmsearch or hmmscan summary table output file
    """
    attribs = ['id', 'bias', 'bitscore', 'description']
    hits = defaultdict(list)
    with open(hmmer_output) as handle:
        for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
            for hit in queryresult.hits:
                for attrib in attribs:
                    hits[attrib].append(getattr(hit, attrib))
    return pd.DataFrame.from_dict(hits)

def filterFASTAbyIDs(input_fasta: str, record_ids: list,
                     output_fasta: str = None) -> None:
    """
    Filter records in fasta file matching provided IDs
    """
    if output_fasta is None:
        basename, ext = os.path.splitext(input_fasta)
        output_fasta = f'{basename}_filtered{ext}'
    hit_records = [record for record in SeqIO.parse(input_fasta, 'fasta')
                   if record.id in set(record_ids)]
    with open(output_fasta, 'w') as out_handle: 
         SeqIO.write(hit_records, out_handle, 'fasta')

def filterFASTAByHMM(hmm_model: str, input_fasta: str,
                     output_fasta: str = None,
                     method: str = 'hmmsearch') -> None:
    """
    Generate protein-specific database by filtering
    sequence database to only contain sequences 
    corresponing to protein of interest
    """
    
    basename, ext = os.path.splitext(input_fasta)
    hmmer_output = f'{basename}.txt'

    runHMMER(hmm_model=hmm_model,
             input_fasta=input_fasta,
             output_file=hmmer_output,
             method=method)

    hmmer_hits = parseHMMERoutput(hmmer_output)

    filterFASTAbyIDs(input_fasta, record_ids=hmmer_hits.id.values,
                     output_fasta=output_fasta)

