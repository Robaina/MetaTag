"""
Tools to perform multiple sequence alignments
"""

import os
import pyfastx
from Bio import AlignIO
from .utils import terminalExecute, setDefaultOutputPath


def runMAFFT(input_fasta: str, output_file: str = None,
             n_threads: int = -1, parallel: bool = True,
             additional_args: str = None) -> None:
    """
    Simple CLI wrapper to mafft (MSA)
    
    Manual: https://mafft.cbrc.jp/alignment/software/manual/manual.html
    
    CLI examples:
    mafft --globalpair --thread n in > out
    mafft --localpair --thread n in > out
    mafft --large --globalpair --thread n in > out
    """
    if output_file is None:
        output_file = setDefaultOutputPath(input_fasta, extension='.fasta.aln')
    if parallel:
        thread_str = f'--thread {n_threads}'
    else:
        thread_str = ''
    if additional_args is None:
        additional_args = ''
    cmd_str = f'mafft {thread_str} {additional_args} {input_fasta} > {output_file}'
    terminalExecute(cmd_str, suppress_output=False)

def runMuscle(input_fasta: str, output_file: str = None,
              maxiters: int = None) -> None:
    """
    Simple CLI wrapper to muscle (MSA)
    muscle: https://www.drive5.com/muscle/manual/output_formats.html

    output phylip and fasta.aln
    """
    if output_file is None:
        output_file = setDefaultOutputPath(input_fasta,
                                           extension='.fasta.aln',
                                           only_filename=True)
    if maxiters is None:
        maxiters = 2
    cmd_str = f'muscle -in {input_fasta} -out {output_file} -maxiters {maxiters}'
    terminalExecute(cmd_str, suppress_output=False)

def runTrimal(input_aln: str, output_aln: str = None) -> None:
    """
    Simple CLI wrapper to trimal
  
    I/O in phylip as well: https://vicfero.github.io/trimal/
    """
    if output_aln is None:
        output_aln = setDefaultOutputPath(input_aln, '_trimal')
    cmd_str = (f'trimal -in {input_aln} -out {output_aln} -fasta -automated1 '
               f'-resoverlap 0.55 -seqoverlap 60 -htmlout trimal.html')
    terminalExecute(cmd_str, suppress_output=False)

def runHMMalign(input_hmm: str, input_aln: str,
                input_seqs: str, 
                output_aln_seqs: str = None,
                additional_args: str = None) -> None:
    """
    Simple CLI wrapper to hmmalign
    Align short read query sequences to reference MSA
    """
    if output_aln_seqs is None:
        output_aln_seqs = setDefaultOutputPath(input_seqs, '_hmm', extension='.aln')
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ''
    cmd_str = (
        f'hmmalign -o {output_aln_seqs} --mapali {input_aln} --trim '
        f'--informat FASTA {args_str} {input_hmm} {input_seqs}'
        )
    terminalExecute(cmd_str, suppress_output=False)

def alignShortReadsToReferenceMSA() -> None:
    """
    Align short read query sequences to reference MSA
    through the same hmm model employed to reconstruct
    reference database
    """
    # runHMMalign()
    pass

def convertFastaAlnToPhylip(input_fasta_aln: str,
                            output_file: str = None) -> None:
    """
    Convert alignments in Fasta to Phylip.
    Note: id labels in phylip format are restricted to 10 characters.
          This restriction may caused repeated truncated labels.
    """
    if output_file is None:
        output_file = setDefaultOutputPath(input_fasta_aln, extension='.phylip')
    with open(input_fasta_aln, 'r') as input_handle, open(output_file, 'w') as output_handle:
        alignments = AlignIO.parse(input_handle, 'fasta')
        AlignIO.write(alignments, output_handle, 'phylip')
        # records = SeqIO.parse(input_handle, 'fasta')
        # SeqIO.write(records, output_handle, 'phylip')

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

    

 