"""
Functions to perform multiple alignment of peptide sequences 
and phylogenetic tree reconstruction
"""

import os
import shutil
from Bio import AlignIO
from .utils import terminalExecute, setDefaultOutputPath


def alignPeptides():
    """
    Perform multiple alignment on a set of peptides

    NOTES:

    script: 0MusclePep.py, 0MusclePepTrimal.py
    """
    pass

def getPhyloTree():
    """
    Make phylogenetic tree out of peptide multiple aligment data

    NOTES:

    script: fasttree.sh, iqtree.sh
    """
    pass

def placeReadsOntoTree():
    """
    Performs short read placement onto phylogenetic tree

    NOTES:

    Runs papara
    """
    pass

def runMAFFT(input_fasta: str, output_file: str = None,
             n_threads: int = -1, parallel: bool = True,
             additional_args: str = None) -> None:
    """
    Simple CLI wrapper to mafft
    
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
    Simple CLI wrapper to muscle
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
    TODO: all
    """
    if output_aln is None:
        output_aln = setDefaultOutputPath(input_aln, '_trimal')
    cmd_str = (f'trimal -in {input_aln} -out {output_aln} -fasta -automated1 '
               f'-resoverlap 0.55 -seqoverlap 60 -htmlout trimal.html')
    terminalExecute(cmd_str, suppress_output=False)

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

def runFastTree(input_algns: str, output_file: str = None,
                nucleotides: bool = False,
                additional_args: str = None) -> None:
    """
    Simple CLI wrapper to fasttree.
    fasttree accepts multiple alignments in fasta or phylip formats

    additional_args: a string containing additional parameters and
                    parameter values to be passed to fasttree
    """
    if output_file is None:
        output_file = setDefaultOutputPath(input_algns, tag='_fasttree', extension='.newick')
    if nucleotides:
        nt_str = '-gtr -nt'
    else:
        nt_str = ''
    if additional_args is None:
        additional_args = ''
    cmd_str = f'fasttree {nt_str} {input_algns} {additional_args} > {output_file}'
    terminalExecute(cmd_str, suppress_output=False)

def runIqTree(input_algns: str, output_dir: str = None,
              output_prefix: str = None,
              keep_recovery_files: bool = False,
              nucleotides: bool = False, n_processes: int = None,
              substitution_model: str = 'TEST',
              bootstrap_replicates: int = 1000,
              additional_args: str = None) -> None:
    """
    Simple CLI wrapper to iqtree.
    iqtree accepts multiple alignments in fasta or phylip formats.

    additional_args: a string containing additional parameters and
                     parameter values to be passed to iqtree

    output: iqtree outputs several files
    """
    def removeAuxiliaryOutput(output_prefix):
        """
        Removes iqtree auxiliary output files
        """
        output_exts = [
            '.bionj', '.ckp.gz', '.contree', '.iqtree', '.log',
            '.model.gz', '.splits.nex', '.treefile', '.mldist'
        ]
        output_files = [output_prefix + ext for ext in output_exts]
        for file_path in output_files:
            if ('iqtree' not in file_path) and ('contree' not in file_path):
                os.remove(file_path)

    if output_dir is None:
        output_dir = os.path.dirname(input_algns)
    else:
        output_dir = os.path.abspath(output_dir)
    if output_prefix is None:
        input_file = os.path.basename(input_algns)
        output_prefix = os.path.join(output_dir, input_file)
        output_prefix_str = f'-pre {output_prefix}'
    else:
        output_prefix = os.path.join(output_dir, output_prefix)
        output_prefix_str = f'-pre {output_prefix}'
    if nucleotides:
        seq_type = 'DNA'
    else:
        seq_type = 'AA'
    if n_processes is None:
        n_processes = 'AUTO'
    if additional_args is None:
        additional_args = ''
    cmd_str = (f'iqtree -s {input_algns} -st {seq_type} -nt {n_processes} '
               f'-m {substitution_model} -bb {bootstrap_replicates} {output_prefix_str} {additional_args}')
    terminalExecute(cmd_str, suppress_output=False)
    if not keep_recovery_files:
        removeAuxiliaryOutput(output_prefix)
    