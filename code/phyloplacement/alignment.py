"""
Tools to perform multiple sequence alignments
"""

import os
import pyfastx
from Bio import AlignIO

import phyloplacement.wrappers as wrappers
from phyloplacement.utils import terminalExecute, setDefaultOutputPath


def alignShortReadsToReferenceMSA() -> None:
    """
    Align short read query sequences to reference MSA
    through the same hmm model employed to reconstruct
    reference database
    """
    # wrappers.runHMMalign()
    pass

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

def convertStockholmToFasta(input_stockholm: str,
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
                                      ref_ids: set,
                                      out_dir: str = None) -> None:
    """
    Separate reference sequences from query sequences in msa fasta file
    """
    if out_dir is None:
        out_dir = os.path.dirname(ref_query_msa)
    out_ref_msa = setDefaultOutputPath(ref_query_msa, tag='_ref_fraction')
    out_query_msa = setDefaultOutputPath(ref_query_msa, tag='_query_fraction')
    
    fasta = pyfastx.Fasta(ref_query_msa, build_index=False, full_name=True)
    with open(out_ref_msa, 'w') as outref, open(out_query_msa, 'w') as outquery:
        for record_name, record_seq in fasta:
            if record_name in ref_ids:
                outref.write(f'>{record_name}\n{record_seq}\n')
            else:
                outquery.write(f'>{record_name}\n{record_seq}\n')

    

 