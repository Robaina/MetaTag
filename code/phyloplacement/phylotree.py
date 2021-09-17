"""
Functions to perform multiple alignment of peptide sequences 
and phylogenetic tree reconstruction
"""

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
    Convert alignments in Fasta to Phylip
    """
    if output_file is None:
        output_file = setDefaultOutputPath(input_fasta_aln, extension='phylip')
    with open(input_fasta_aln, 'r') as input_handle, open(output_file, 'w') as output_handle:
        alignments = AlignIO.parse(input_handle, 'fasta')
        AlignIO.write(alignments, output_handle, 'phylip')

def runFastTree(input_phylip: str, output_file: str = None,
                nucleotides: bool = False) -> None:
    """
    Simple CLI wrapper to fasttree
    """
    if output_file is None:
        output_file = setDefaultOutputPath(input_phylip, tag='_tree', extension='txt')
    if nucleotides:
        nt_str = '-gtr -nt'
    else:
        nt_str = ''
    cmd_str = f'fasttree {nt_str} {input_phylip} > {output_file}'
    terminalExecute(cmd_str, suppress_output=False)