"""
Functions to perform multiple alignment of peptide sequences 
and phylogenetic tree reconstruction
"""

import os
import shutil
import tempfile
from .utils import terminalExecute, setDefaultOutputPath
from .database import runHMMbuild
from .alignment import runHMMalign, convertStockholmToFasta


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

def placeReadsOntoTree(input_tree: str, input_aln: str,
                       input_query: str,
                       query_already_aligned: bool = False,
                       output_file: str = None) -> None:
    """
    Performs short read placement onto phylogenetic tree

    NOTES:

    Runs hmmbuild, hmmalign, epa-ng
    """
    if output_file is None:
        output_file = setDefaultOutputPath(input_query, tag='_placement',
                                           extension='.jplace')
    if not query_already_aligned:
         
        output_hmm = '/home/robaina/Documents/TRAITS/out.hmm'
        runHMMbuild(input_aln=input_aln,
                    output_hmm=output_hmm,
                    additional_args=None)
        
        output_aln_seqs = '/home/robaina/Documents/TRAITS/aln_query.stk'
        runHMMalign(input_aln=input_aln,
                    input_hmm=output_hmm,
                    input_seqs=input_query,
                    output_aln_seqs=output_aln_seqs,
                    additional_args=None)
        
        input_aln_query = '/home/robaina/Documents/TRAITS'
        convertStockholmToFasta(input_stockholm=output_aln_seqs,
                                output_fasta=input_aln_query)
    else:
        input_aln_query = input_query
    
    output_dir = '/home/robaina/Documents/TRAITS'
    runEPAng(input_tree=input_tree, input_aln=input_aln,
             input_query=input_aln_query, output_dir=output_dir,
             n_threads=None, additional_args=None)

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
        output_file = setDefaultOutputPath(input_algns, tag='_fasttree',
                                           extension='.newick')
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
        exts_to_remove = [
            '.bionj', '.ckp.gz', '.iqtree', '.log',
            '.model.gz', '.splits.nex', '.mldist'
        ]
        files_to_remove = [output_prefix + ext for ext in exts_to_remove]
        for file_path in files_to_remove:
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

def runTreeShrink(input_tree: str, input_aln: str,
                  output_dir: str = None,
                  output_deleted_nodes: bool = False,
                  additional_args: str = None) -> None: 
    """
    Run treeshrink to remove tree branch outliers. 
    Remove outliers from MSA file too.
    Tree file must be of newick format.
    see run_treeshrink.py  -h for help
    """
    if output_dir is None:
        output_dir = os.path.dirname(input_tree)
    else:
        output_dir = os.path.abspath(output_dir)
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ''
    
    out_tree = setDefaultOutputPath(input_tree, tag='_shrink', only_filename=True)
    out_aln = setDefaultOutputPath(input_aln, tag='_shrink', only_filename=True)

    # Handle treeshrink input/output requirements (temp/tree/input.tree)
    with tempfile.TemporaryDirectory() as temp_out_dir, \
         tempfile.TemporaryDirectory() as parent_in_temp, \
         tempfile.TemporaryDirectory(dir=parent_in_temp) as temp_in_dir:

        temp_tree_dir = os.path.basename(temp_in_dir)
        shutil.copy(input_tree, os.path.join(temp_in_dir, "input.tree"))
        shutil.copy(input_aln, os.path.join(temp_in_dir, "input.aln"))
        
        cmd_str = (
            f'run_treeshrink.py -i {parent_in_temp} -m per-gene '
            '-t input.tree -a input.aln '
            f'-o {temp_out_dir} -O output {args_str}'
                )
        terminalExecute(cmd_str, suppress_output=False)

        shutil.move(
            os.path.join(temp_out_dir, temp_tree_dir, "output.tree"),
            os.path.join(output_dir, out_tree)
        )
        shutil.move(
            os.path.join(temp_out_dir, temp_tree_dir, "output.aln"),
            os.path.join(output_dir, out_aln)
            )
        if output_deleted_nodes:
            out_txt = setDefaultOutputPath(input_aln, tag='_shrink_deleted',
                                           extension='.txt', only_filename=True)
            shutil.move(
                os.path.join(temp_out_dir, temp_tree_dir, "output.txt"),
                os.path.join(output_dir, out_txt)
            )

def runPapara() -> None:
    """
    Simple CLI wrapper to Papara
    papara -t tree.nwk -s alignment.phy -q query-seqs.fasta -r -n combined-aln (name of output alignment)
    
    -r 	Prevent PaPaRa from adding gaps in the reference alignment

    '-j <num threads>'

    Seems like EPA-NG is better, quicker and developed by same tema that made Papara..
    """
    pass

def runEPAng(input_tree: str, input_aln: str, input_aln_query: str,
             model: str = None, output_dir: str = None,
             n_threads: int = None,
             additional_args: str = None) -> None:
    """
    Simple CLI wrapper to EPA-ng
    See epa-ng -h for additional parameters

    input_tree: newick format
    input_aln: fasta format
    input_aln_query: fasta format (sequences must be alignned to reference 
    msa fasta and have the same length as the reference msa alignment)
    """
    if model is None:
        model = 'GTR+G'
    if output_dir is None:
        output_dir = os.path.dirname(input_tree)
    else:
        output_dir = os.path.abspath(output_dir)
    if n_threads is None:
        n_threads = os.cpu_count() - 1
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ''

    cmd_str = (
        f'epa-ng --ref-msa {input_aln} --tree {input_tree} --query {input_aln_query} '
        f'--model {model} --threads {n_threads} --outdir {output_dir} {args_str}'
        )
    terminalExecute(cmd_str, suppress_output=False)


        
    
