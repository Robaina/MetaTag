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
    """
    if output_file is None:
        output_file = setDefaultOutputPath(input_fasta, extension='.fasta.aln')
    if maxiters is None:
        maxiters = 2
    cmd_str = f'muscle -in {input_fasta} -out {output_file} -maxiters {maxiters}'
    terminalExecute(cmd_str, suppress_output=False)

def runTrimal(input_aln: str, output_aln: str = None) -> None:
    """
    Simple CLI wrapper to trimal
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
    terminalExecute(cmd_str, suppress_output=True)

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
        output_file = setDefaultOutputPath(input_fasta_aln, extension='phylip')
    with open(input_fasta_aln, 'r') as input_handle, open(output_file, 'w') as output_handle:
        alignments = AlignIO.parse(input_handle, 'fasta')
        AlignIO.write(alignments, output_handle, 'phylip')

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