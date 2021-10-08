"""
Simple CLI wrappers to several tools
"""

import os
import shutil
import tempfile
from .utils import terminalExecute, setDefaultOutputPath
path_to_papara_exec = '/home/robaina/Software/papara'


def runHMMsearch(hmm_model: str, input_fasta: str,
             output_file: str = None,
             method: str = 'hmmsearch',
             n_processes: int = None) -> None:
    """
    Simple CLI wrapper to hmmsearch or hmmscan
    Requires hmmer installed and accessible
    """
    if n_processes is None:
        n_processes = os.cpu_count() - 1
    if output_file is None:
        output_file = setDefaultOutputPath(input_fasta, '_hmmer_hits', '.txt')
    cmd_str = (f'{method} --cut_ga --tblout {output_file} --cpu {n_processes} '
               f'{hmm_model} {input_fasta}')
    terminalExecute(cmd_str, suppress_output=False)

def runHMMbuild(input_aln: str, output_hmm: str = None,
                additional_args: str = None) -> None:
    """
    Simple CLI wrapper to hmmbuild (build HMM profile from MSA file)
    additional args: see hmmbuild -h
    """
    if output_hmm is None:
        output_hmm = setDefaultOutputPath(input_aln, extension='.hmm')
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ''
    cmd_str = f'hmmbuild {args_str} {output_hmm} {input_aln}'
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

def runCDHIT(input_fasta: str, output_fasta: str = None,
             additional_args: str = None) -> None:
    """
    Simple CLI wrapper to cd-hit to obain representative sequences
    CD-HIT may be used to remove duplicated sequences (keeps one representatie)
    with parameters -c 1 -t 1. However, it does require lots of RAM to store sequences,
    cannot run on Aquifex.

    """
    if output_fasta is None:
       output_fasta = setDefaultOutputPath(input_fasta, '_cdhit')
    if additional_args is None:
        additional_args = ''
    cmd_str = f'cd-hit -i {input_fasta} -o {output_fasta} {additional_args}'
    terminalExecute(cmd_str, suppress_output=False)

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
              maxiters: int = None,
              additional_args: str = None) -> None:
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
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ''
    cmd_str = (f'muscle -in {input_fasta} -out {output_file} '
               f'-maxiters {maxiters} {args_str}')
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

def runPapara(tree_nwk: str, msa_phy: str,
              query_fasta: str, n_threads: int = None,
              output_file: str = None,
              additional_args: str = None) -> None:
    """
    Simple CLI wrapper to Papara
    papara -t tree.nwk -s alignment.phy -q query-seqs.fasta -r -n combined-aln (name of output alignment)
    
    -r 	Prevent PaPaRa from adding gaps in the reference alignment

    '-j <num threads>'

    -a: sequences are protein data

    Run Papara to do query alignment to reference MSA and tree (required for EPA-ng)
    Alignment could be done with hmmalign or muscle as well, but these tools don't 
    consider the tree during alignment (would this be a justified improvement over hmmalign?)

    cd /home/robaina/Software/papara
    ./papara (to run papara)

    There seems to be a problem with enabling multithreading in papara when run as a static
    executable. It looks like it has to be enable during compilation (but compilation currently not working):
    https://stackoverflow.com/questions/19618926/thread-doesnt-work-with-an-error-enable-multithreading-to-use-stdthread-ope
    """
    if n_threads is not None:
        threads_str = f'-j {n_threads}'
    else:
        threads_str = ''
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ''

    cmd_str = (
        f'{os.path.join(path_to_papara_exec, "papara")} -t {tree_nwk} '
        f'-s {msa_phy} -q {query_fasta} {threads_str} -n phylip '
        f'-r {args_str} -a'
        )
    terminalExecute(cmd_str, suppress_output=False)

def runEPAng(input_tree: str, input_aln_ref: str, input_aln_query: str,
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

    epa-ng: https://github.com/Pbdas/epa-ng

    TODO: model parameters must be the same employed to make reference tree.
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
        f'epa-ng --ref-msa {input_aln_ref} --tree {input_tree} --query {input_aln_query} '
        f'--model {model} --threads {n_threads} --outdir {output_dir} {args_str}'
        )
    terminalExecute(cmd_str, suppress_output=False)

def runGappaHeatTree(input_jplace: str,
                          output_dir: str = None,
                          output_prefix: str = None, 
                          additional_args: str = None) -> None:
    """
    Run gappa examine heat-tree to obtain tree with short read placements 
    in newick format
    """
    if output_dir is None:
        outdir_str = ''
    else:
        outdir_str = f'--out-dir {os.path.abspath(output_dir)}'
    if output_prefix is None:
        output_prefix_str = f''
    else:
        output_prefix_str = f'--file-prefix {output_prefix}'
    if additional_args is None:
        args_str = ''
    else:
        args_str = additional_args

    cmd_str = (
        f'gappa examine heat-tree --jplace-path {os.path.abspath(input_jplace)} '
        f'--write-newick-tree --write-svg-tree '
        f'{outdir_str} {output_prefix_str} {args_str}'
        )
    terminalExecute(cmd_str, suppress_output=False)




