"""
Tools to create peptide-specific sequence databases

1. Implement hmmr
2. Filter fasta files based on query sequences
"""

import os
import pandas as pd
from collections import defaultdict
from Bio import SearchIO, SeqIO, AlignIO
import pyfastx

import phyloplacement.wrappers as wrappers
from phyloplacement.utils import setDefaultOutputPath


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

def convertFastaAlnToPhylip(input_fasta_aln: str,
                            output_file: str = None) -> None:
    """
    Convert alignments in Fasta to Phylip.
    """
    if output_file is None:
        output_file = setDefaultOutputPath(input_fasta_aln, extension='.phylip')
    with open(input_fasta_aln, 'r') as input_handle, open(output_file, 'w') as output_handle:
        alignments = AlignIO.parse(input_handle, 'fasta')
        AlignIO.write(alignments, output_handle, 'phylip')

def convertPhylipToFastaAln(input_phylip: str,
                            output_file: str = None) -> None:
    """
    Convert alignments in Phylip to Fasta format
    """
    if output_file is None:
        output_file = setDefaultOutputPath(input_phylip, extension='.fasta.aln')
    with open(input_phylip, 'r') as input_handle, open(output_file, 'w') as output_handle:
        alignments = AlignIO.parse(input_handle, 'phylip')
        AlignIO.write(alignments, output_handle, 'fasta')

def convertStockholmToFastaAln(input_stockholm: str,
                            output_fasta: str = None) -> None:
    """
    Convert alignment file in Stockholm format to fasta
    """
    if output_fasta is None:
        output_fasta = setDefaultOutputPath(input_stockholm,
                                            extension='.fasta')
    with open(output_fasta, 'w') as fasta_file:
        align = AlignIO.read(input_stockholm, 'stockholm')
        print(align.format('fasta'), file=fasta_file)

def splitReferenceFromQueryAlignments(ref_query_msa: str,
                                      ref_ids: set = None,
                                      ref_prefix: str = None,
                                      out_dir: str = None) -> None:
    """
    Separate reference sequences from query sequences in msa fasta file
    """
    if out_dir is None:
        out_dir = os.path.dirname(ref_query_msa)
    if ref_ids is not None and ref_prefix is None:
        def is_reference(record_name): return record_name in ref_ids
    elif ref_ids is None and ref_prefix is not None:
        def is_reference(record_name): return record_name.startswith(ref_prefix)
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

def countRecords(fasta_file: str) -> None:
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