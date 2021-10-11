"""
Tools to create peptide-specific sequence databases

1. Implement hmmr
2. Filter fasta files based on query sequences
"""

import os
import pandas as pd
from collections import defaultdict
from Bio import SearchIO, SeqIO
import pyfastx

import phyloplacement.wrappers as wrappers
from phyloplacement.utils import terminalExecute, setDefaultOutputPath


def removeDuplicatesFromFastaByID(input_fasta: str,
                                  output_fasta: str = None) -> None:
    """
    Remove entries with duplicated IDs from fasta.
    """
    if output_fasta is None:
        output_fasta = setDefaultOutputPath(input_fasta, '_noduplicates')
    seen_ids = set()
    records = []
    for record in SeqIO.parse(input_fasta, "fasta"):  
        if record.id not in seen_ids:
            seen_ids.add(record.id)
            records.append(record)
    with open(output_fasta, 'w') as out_handle: 
         SeqIO.write(records, out_handle, 'fasta')

def removeDuplicatesFromFasta(input_fasta: str,
                              output_fasta: str = None,
                              output_duplicates: bool = False) -> None:
    """
    Removes duplicate entries (either by sequence or ID) from fasta.
    TODO: implement output_duplicates
    """
    if output_fasta is None:
        output_fasta = setDefaultOutputPath(input_fasta, '_noduplicates')
    seen_seqs, seen_ids = set(), set()
    records = []
    for record in SeqIO.parse(input_fasta, "fasta"):  
        if (record.seq not in seen_seqs) and (record.id not in seen_ids):
            seen_seqs.add(record.seq)
            seen_ids.add(record.id)
            records.append(record)
    with open(output_fasta, 'w') as out_handle: 
         SeqIO.write(records, out_handle, 'fasta')

def filterFastaBySequenceLength(input_fasta: str, minLength: int = 0,
                                maxLength: int = None,
                                output_fasta: str = None) -> None:
    """
    Filter sequences by length in fasta file
    """     
    fa = pyfastx.Fasta(input_fasta)
    record_ids = fa.keys()
    if maxLength is not None:
        max_tag = str(maxLength)
        record_ids.filter(record_ids>=minLength, record_ids<=maxLength)
    else:
        max_tag = ''
        record_ids.filter(record_ids>=minLength)
    if output_fasta is None:
        output_fasta = setDefaultOutputPath(input_fasta, f'_length_{minLength}_{max_tag}')
    with open(output_fasta, 'w') as fp:
        for record_id in record_ids:
            record_obj = fa[record_id]
            fp.write(record_obj.raw)

def mergeFASTAs(input_fastas_dir: list, output_fasta: str = None) -> None:
    """
    Merge input fasta files into a single fast
    """
    if output_fasta is None:
        output_fasta = os.path.join(input_fastas_dir, 'merged.fasta')
    cmd_str = f'awk 1 *.fasta > {output_fasta}'
    terminalExecute(cmd_str, suppress_output=False)

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

def filterFASTAByHMM(hmm_model: str, input_fasta: str,
                     output_fasta: str = None,
                     method: str = 'hmmsearch',
                     remove_uninformative: bool = False) -> None:
    """
    Generate protein-specific database by filtering
    sequence database to only contain sequences 
    corresponing to protein of interest
    """
    basename, ext = os.path.splitext(input_fasta)
    hmm_name, _ = os.path.splitext(os.path.basename(hmm_model))
    hmmer_output = f'{basename}_{hmm_name}.txt'
    
    print('Running Hmmer...')
    wrappers.runHMMsearch(
        hmm_model=hmm_model,
        input_fasta=input_fasta,
        output_file=hmmer_output,
        method=method
        )
    print('Parsing Hmmer output file...')
    hmmer_hits = parseHMMsearchOutput(hmmer_output)
    print('Filtering Fasta...')
    filterFASTAbyIDs(input_fasta, record_ids=hmmer_hits.id.values,
                     output_fasta=output_fasta)
    if remove_uninformative:
        wrappers.runCDHIT(input_fasta=output_fasta)