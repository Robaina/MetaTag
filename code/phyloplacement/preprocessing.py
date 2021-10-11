"""
Tools to preprocess sequence databases

1. Remove illegal characters from peptide sequences
2. Remove illegal symbols from file paths
3. Relabel fasta records and make dictionary with old labels
"""

import os
import re
import shutil
from Bio import SeqIO
import pyfastx

from .utils import saveToPickleFile, setDefaultOutputPath


def relabelRecordsInFASTA(input_fasta: str,
                 output_dir: str = None,
                 prefix: str = None):
    """
    Change record ids for numbers and store then in a dictionary
    Useful when converting fasta alignments to phylip format 
    (phylip ids are less than 10 char long)
    """
    if output_dir is None:
        output_dir = os.path.dirname(input_fasta)
    if prefix is not None:
        prefix_str = prefix 
    else:
        prefix_str = ''
    
    fasta_file = setDefaultOutputPath(input_fasta, 
                                      tag='_short_ids',
                                      only_filename=True)
    dict_file = setDefaultOutputPath(input_fasta,
                                     tag='_id_dict',
                                     extension='.pickle',
                                     only_filename=True)
    output_fasta = f'{os.path.join(output_dir, fasta_file)}'
    output_dict = f'{os.path.join(output_dir, dict_file)}'
    
    fasta = SeqIO.parse(input_fasta, 'fasta')
    # new_ids = map(lambda n: f'{prefix_str}{n}', range(len(fasta)))

    id_dict = dict()
    with open(output_fasta, 'w') as outfasta:
        for n, record in enumerate(fasta):
            new_id = f'{prefix_str}{n}'
            id_dict[new_id] = record.description
            outfasta.write(f'>{new_id}\n{record.seq}\n')
    saveToPickleFile(id_dict, output_dict)

def reformatFilePath(file_name: str) -> None:
    """
    Remove illegal symbols from file path
    """
    upper_lower_digits = re.compile('[^a-zA-Z0-9]')
    fdir = os.path.dirname(file_name)
    fname, ext = os.path.splitext(os.path.basename(file_name))
    clean_fname = upper_lower_digits.sub(
        '_', fname).replace('__', '_').strip('_')
    return os.path.join(fdir, f'{clean_fname}{ext}')

def isLegitPeptideSequence(record_seq: str) -> bool:
    """
    Assert that peptide sequence only contains valid symbols
    """
    aas = {
        'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
        'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
    }
    seq_symbols = {s for s in record_seq}
    return seq_symbols.issubset(aas)

def isLegitDNAsequence(record_seq: str) -> bool:
    """
    Assert that DNA sequence only contains valid symbols
    """
    nts = {'A', 'G', 'T', 'C'}
    seq_symbols = {s for s in record_seq}
    return seq_symbols.issubset(nts)

# def reformatPeptideSequence(record_seq) -> str:
#     """
#     Assert peptide sequence only contains upper case letters
#     """
#     upper_case_letters = re.compile('[^A-Z]')
#     return upper_case_letters.sub('', record_seq.upper())

def reformatSequencesInFASTA(fasta_file: str,
                             output_file: str = None,
                             is_peptide: bool = True) -> None:
    """
    Filter out (DNA or peptide) sequences containing illegal characters
    """
    dirname = os.path.dirname(fasta_file)
    basename = os.path.basename(fasta_file)
    fname, ext = os.path.splitext(basename)

    if output_file is None:
        output_file = os.path.join(dirname, f'{fname}_modified{ext}')
    else:
        output_file = os.path.abspath(output_file)
    if is_peptide:
        isLegitSequence = isLegitPeptideSequence
    else:
        isLegitSequence = isLegitDNAsequence

    fasta = pyfastx.Fasta(fasta_file, build_index=False, full_name=True)
    with open(output_file, 'w') as outfile:
        for record_name, record_seq in fasta:
            if isLegitSequence(record_seq):
                outfile.write(f'>{record_name}\n{record_seq}\n')

def pipe_line(fasta_path: str, id_type: int,
              output_file = None) -> None:
    """
    Pipeline to clean sequences and path to fasta file
    """
    def is_legit_path(fasta_path, legit_fasta_path):
        return fasta_path == legit_fasta_path
            
    clean_fasta_path = reformatFilePath(fasta_path)

    if not is_legit_path(fasta_path, clean_fasta_path):
        shutil.move(fasta_path, clean_fasta_path)
    
    reformatSequencesInFASTA(fasta_file=clean_fasta_path,
                             output_file=output_file)