"""
Tools to perform multiple sequence alignments
"""

import os
import tempfile
import pyfastx
from Bio import AlignIO

import phyloplacement.wrappers as wrappers
from phyloplacement.utils import setDefaultOutputPath
from phyloplacement.alignment import convertFastaAlnToPhylip, convertPhylipToFastaAln


def alignShortReadsToReferenceMSA(ref_msa: str, query_seqs: str,
                                  method: str = 'papara',
                                  tree_nwk: str = None,
                                  output_dir: str = None) -> None:
    """
    Align short read query sequences to reference MSA (fasta format).
    Outputs fasta msa alignment between query and reference sequences
    """
    if output_dir is None:
        output_dir = setDefaultOutputPath(ref_msa, only_dirname=True)
    output_hmm = os.path.join(
        output_dir, setDefaultOutputPath(ref_msa, extension='.hmm', only_filename=True)
        )
    output_aln_seqs = os.path.join(
        output_dir, setDefaultOutputPath(query_seqs, extension='.aln',
                                         only_filename=True)
    )
    
    if method.lower() in 'hmmalign':
        wrappers.runHMMbuild(input_aln=ref_msa,
                             output_hmm=output_hmm)
        
        with tempfile.TemporaryFile() as tempstock:
            wrappers.runHMMalign(input_aln=ref_msa,
                                 input_hmm=output_hmm,
                                 input_seqs=query_seqs,
                                 output_aln_seqs=tempstock)
            convertStockholmToFastaAln(input_stockholm=tempstock,
                                       output_fasta=output_aln_seqs)
    elif method.lower() in 'papara':
        with tempfile.TemporaryFile() as tempphy, \
             tempfile.TemporaryFile() as tempqueryaln:
             convertFastaAlnToPhylip(input_fasta_aln=output_aln_seqs,
                                     output_file=tempphy)
             wrappers.runPapara(tree_nwk=tree_nwk,
                                msa_phy=tempphy,
                                query_fasta=query_seqs,
                                output_aln=tempqueryaln)
             convertPhylipToFastaAln(input_phylip=tempqueryaln,
                                     output_file=output_aln_seqs)
    else:
        raise ValueError('Alignment method not implemented')
        
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

    

 