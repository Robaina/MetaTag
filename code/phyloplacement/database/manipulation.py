#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to create peptide-specific sequence databases

1. Implement hmmr
2. Filter fasta files based on query sequences
"""

import os
import warnings
from collections import defaultdict

import pandas as pd
from Bio import SearchIO, SeqIO, AlignIO
import pyfastx

import phyloplacement.wrappers as wrappers
from phyloplacement.utils import setDefaultOutputPath


def filterFastaBySequenceLength(input_fasta: str, minLength: int = None,
                                maxLength: int = None,
                                output_fasta: str = None) -> None:
    """
    Filter sequences by length in fasta file
    """  
    if (minLength is None) and (maxLength is None):
        warnings.warn("Missing boundary values for sequence length")
        return
    input_fasta = os.path.abspath(input_fasta)   
    fa = pyfastx.Fasta(input_fasta)
    record_ids = fa.keys()
    if minLength is None:
        minLength = 0
    if maxLength is not None:
        max_tag = str(maxLength)
        record_ids.filter(record_ids>=minLength, record_ids<=maxLength)
    else:
        max_tag = ''
        record_ids.filter(record_ids>=minLength)
    if output_fasta is None:
        output_fasta = setDefaultOutputPath(input_fasta, f'_length_{minLength}_{max_tag}')
    if not record_ids:
        raise ValueError("No records found with given sequence length bounds")
    with open(output_fasta, 'w') as fp:
        for record_id in record_ids:
            record_obj = fa[record_id]
            fp.write(record_obj.raw)
    os.remove(input_fasta + ".fxi")

def parseHMMsearchOutput(hmmer_output: str) -> pd.DataFrame:
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
       output_fasta = setDefaultOutputPath(input_fasta, '_fitered')
    record_ids = set(record_ids)
    fa = pyfastx.Fasta(input_fasta)
    with open(output_fasta, 'w') as fp:
        for record_id in record_ids:
            try:
                record_obj = fa[record_id]
                fp.write(record_obj.raw)
            except:
                pass
    os.remove(input_fasta + ".fxi")

def filterFASTAByHMM(hmm_model: str, input_fasta: str,
                     output_fasta: str = None,
                     hmmer_output: str = None,
                     method: str = 'hmmsearch',
                     additional_args: str = None) -> None:
    """
    Generate protein-specific database by filtering
    sequence database to only contain sequences 
    corresponing to protein of interest
    """
    hmm_name, _ = os.path.splitext(os.path.basename(hmm_model))
    if hmmer_output is None:
        hmmer_output = setDefaultOutputPath(input_fasta, tag=f'_{hmm_name}', extension='.txt')
    if output_fasta is None:
        output_fasta = setDefaultOutputPath(input_fasta, tag=f'filtered_{hmm_name}')
    
    print('Running Hmmer...')
    wrappers.runHMMsearch(
        hmm_model=hmm_model,
        input_fasta=input_fasta,
        output_file=hmmer_output,
        method=method,
        additional_args=additional_args
        )
    print('Parsing Hmmer output file...')
    hmmer_hits = parseHMMsearchOutput(hmmer_output)
    if not hmmer_hits.id.values:
        raise ValueError('No records found in database matching provided hmm')
    print('Filtering Fasta...')
    filterFASTAbyIDs(input_fasta, record_ids=hmmer_hits.id.values,
                     output_fasta=output_fasta)

def convertFastaAlnToPhylip(input_fasta_aln: str,
                            output_phylip: str = None) -> None:
    """
    Convert alignments in Fasta to Phylip.
    """
    if output_phylip is None:
        output_phylip = setDefaultOutputPath(input_fasta_aln, extension='.phylip')
    with open(input_fasta_aln, 'r') as input_handle, open(output_phylip, 'w') as output_handle:
        alignments = AlignIO.parse(input_handle, 'fasta')
        AlignIO.write(alignments, output_handle, 'phylip-relaxed')

def convertPhylipToFastaAln(input_phylip: str,
                            output_file: str = None) -> None:
    """
    Convert alignments in Phylip to Fasta format
    """
    if output_file is None:
        output_file = setDefaultOutputPath(input_phylip, extension='.faln')
    alignments = AlignIO.parse(input_phylip, 'phylip-relaxed')
    AlignIO.write(alignments, output_file, 'fasta')

def convertStockholmToFastaAln(input_stockholm: str,
                               output_fasta: str = None) -> None:
    """
    Convert alignment file in Stockholm format to fasta
    """
    if output_fasta is None:
        output_fasta = setDefaultOutputPath(input_stockholm, extension='.faln')
    alignments = AlignIO.read(input_stockholm, 'stockholm')
    AlignIO.write(alignments, output_fasta, 'fasta')

def splitReferenceFromQueryAlignments(ref_query_msa: str,
                                      ref_ids: set = None,
                                      ref_prefix: str = None,
                                      out_dir: str = None) -> None:
    """
    Separate reference sequences from query sequences in msa fasta file
    """
    if out_dir is None:
        out_dir = os.path.dirname(ref_query_msa)
    if (ref_ids is not None) and (ref_prefix is None):
        def is_reference(record_name): return record_name in ref_ids
    elif (ref_ids is None) and (ref_prefix is not None):
        def is_reference(record_name): return record_name.lower().startswith(ref_prefix.lower())
    else:
        raise ValueError('Provide either set of ref ids or ref prefix')
    out_ref_msa = setDefaultOutputPath(ref_query_msa, tag='_ref_fraction')
    out_query_msa = setDefaultOutputPath(ref_query_msa, tag='_query_fraction')
    
    fasta = pyfastx.Fasta(ref_query_msa, build_index=False, full_name=True)
    with open(out_ref_msa, 'w') as outref, open(out_query_msa, 'w') as outquery:
        for record_name, record_seq in fasta:
            if is_reference(record_name):
                outref.write(f'>{record_name}\n{record_seq}\n')
            else:
                outquery.write(f'>{record_name}\n{record_seq}\n')

def getFastaRecordIDs(fasta_file: str) -> set:
    """
    Extract record ids from fasta
    """
    fasta = pyfastx.Fasta(fasta_file, full_name=True)
    return set(fasta.keys())

def countRecordsInFasta(fasta_file: str) -> int:
    """
    Count records in fasta file
    """
    with open(fasta_file, 'r') as file:
        n_records = sum(
            [1 for line in file if line.startswith(">")]
            )
    return n_records

def sliceFasta(input_file, output_file, N):
    n = 0
    records = SeqIO.parse(input_file, 'fasta')
    sliced_records = []
    for record in records:
        if n < N:
            sliced_records.append(record)
        else:
            break
        n += 1
    SeqIO.write(sliced_records, output_file, 'fasta')

def is_empty_fasta(fasta_file):
    with open(fasta_file, 'r') as file:
        return '>' not in file.read()