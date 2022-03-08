#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to preprocess sequence databases

1. Remove illegal characters from peptide sequences
2. Remove illegal symbols from file paths
3. Relabel fasta records and make dictionary with old labels
"""

import os
import re

from Bio import SeqIO
import pyfastx


import phyloplacement.wrappers as wrappers
from phyloplacement.utils import (saveToPickleFile, setDefaultOutputPath,
                                  terminalExecute, handle_exceptions)


@handle_exceptions
def removeDuplicatesFromFasta(input_fasta: str,
                              output_fasta: str = None,
                              export_duplicates: bool = False,
                              method: str = 'seqkit') -> None:
    """
    Removes duplicate entries (either by sequence or ID) from fasta.
    """
    if output_fasta is None:
        output_fasta = setDefaultOutputPath(input_fasta, '_noduplicates')

    if 'bio' in method:
        seen_seqs, seen_ids = set(), set()
        def unique_records():
            for record in SeqIO.parse(input_fasta, 'fasta'):  
                if (record.seq not in seen_seqs) and (record.id not in seen_ids):
                    seen_seqs.add(record.seq)
                    seen_ids.add(record.id)
                    yield record

        SeqIO.write(unique_records(), output_fasta, 'fasta')

    else:
        wrappers.runSeqKitNoDup(input_fasta=input_fasta, output_fasta=output_fasta,
                                export_duplicates=export_duplicates)

def mergeFASTAs(input_fastas_dir: list, output_fasta: str = None) -> None:
    """
    Merge input fasta files into a single fasta
    """
    if output_fasta is None:
        output_fasta = os.path.join(input_fastas_dir, 'merged.fasta')
    cmd_str = f'awk 1 * > {output_fasta}'
    terminalExecute(
        cmd_str,
        work_dir=input_fastas_dir,
        )

def assertCorrectFilePath(file_name: str) -> None:
    """
    Remove illegal symbols from file path
    """
    upper_lower_digits = re.compile('[^a-zA-Z0-9]')
    fdir = os.path.dirname(file_name)
    fname, ext = os.path.splitext(os.path.basename(file_name))
    clean_fname = upper_lower_digits.sub(
        '_', fname).replace('__', '_').strip('_')
    return os.path.join(fdir, f'{clean_fname}{ext}')

def fastaContainsNucleotideSequences(fasta_file: str) -> bool:
    """
    Check whether fasta file contains nucleotide sequences
    """
    seq_1 = next(SeqIO.parse(fasta_file, 'fasta'))
    return set(seq_1.seq).issubset({'A', 'G', 'T', 'C'})

def isLegitPeptideSequence(record_seq: str) -> bool:
    """
    Assert that peptide sequence only contains valid symbols
    """
    aas = {
        'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
        'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
    }
    seq_symbols = {s.upper() for s in record_seq}
    return seq_symbols.issubset(aas)

def isLegitDNAsequence(record_seq: str) -> bool:
    """
    Assert that DNA sequence only contains valid symbols
    """
    nts = {'A', 'G', 'T', 'C'}
    seq_symbols = {s.upper() for s in record_seq}
    return seq_symbols.issubset(nts)

@handle_exceptions
def assertCorrectSequenceFormat(fasta_file: str,
                                output_file: str = None,
                                is_peptide: bool = True,
                                remove: bool = True) -> None:
    """
    Filter out (DNA or peptide) sequences containing illegal characters
    """
    dirname = os.path.dirname(fasta_file)
    basename = os.path.basename(fasta_file)
    fname, ext = os.path.splitext(basename)

    def removeStopCodonSignals(record_seq: str) -> str:
        return record_seq.replace('*', '')

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
            if is_peptide:
                record_seq = removeStopCodonSignals(record_seq)
            if isLegitSequence(record_seq):
                outfile.write(f'>{record_name}\n{record_seq}\n')

def setTempRecordIDsInFASTA(input_fasta: str,
                            output_dir: str = None,
                            prefix: str = None):
    """
    Change record ids for numbers and store then in a dictionary
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
    id_dict = dict()
    with open(output_fasta, 'w') as outfasta:
        for n, record in enumerate(fasta):
            new_id = f'{prefix_str}{n}'
            id_dict[new_id] = record.description
            outfasta.write(f'>{new_id}\n{record.seq}\n')
    saveToPickleFile(id_dict, output_dict)

def setOriginalRecordIDsInFASTA(input_fasta: str,
                                label_dict: dict = None,
                                output_fasta: str = None):
    """
    Relabel temporary record ID by original IDs
    """
    if output_fasta is None:
        output_fasta = setDefaultOutputPath(input_fasta, 
                                            tag='_original_ids')
    def relabel_records():
        for record in SeqIO.parse(input_fasta, 'fasta'):
            if record.name in label_dict.keys():  
                name = label_dict[record.name]
                record.name = name
                record.id = name
                record.description = name
            yield record

    SeqIO.write(relabel_records(), output_fasta, 'fasta')

def writeRecordNamesToFile(input_fasta: str, output_file: str = None):
    """
    Write a txt file containing a list of record IDs in fasta
    """
    if output_file is None:
        output_file = setDefaultOutputPath(input_fasta, extension='.txt')
    with open(output_file, 'w') as file:
        lines = []
        records = SeqIO.parse(input_fasta, 'fasta')
        for record in records:
            lines.append(record.name + '\n')
        file.writelines(lines)

def fastq2fasta(input_fastq: str, output_fasta: str = None,
                unzip: bool = True) -> None:
    """
    Convert Fastq to FASTA format via sed
    """
    if output_fasta is None:
        output_fasta = setDefaultOutputPath(input_fastq, extension='.fasta')
    
    if unzip:
        input_uncompressed = input_fastq.strip(".gz")
        cmd_str = f"gzip -d {input_fastq} > {input_uncompressed}"
        terminalExecute(cmd_str, suppress_shell_output=False)
        input_fastq = input_uncompressed

    cmd_str = f"sed -n '1~4s/^@/>/p;2~4p' {input_fastq} > {output_fasta}"
    terminalExecute(cmd_str, suppress_shell_output=False)

def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta) 