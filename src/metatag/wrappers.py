#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Simple CLI wrappers to several tools
"""

import os
import shutil
import tempfile
from pathlib import Path

from metatag.utils import set_default_output_path, terminal_execute

parent_dir = Path(__file__).parent
papara_exec = None
if papara_exec is None:
    papara_bin = parent_dir / "vendor" / "papara_static_x86_64"
    if not papara_bin.is_file():
        raise FileExistsError("Papara executable not found")
    else:
        papara_exec = papara_bin.as_posix()


def run_seqkit_nodup(
    input_fasta: str,
    output_fasta: str = None,
    export_duplicates: bool = False,
    duplicates_file: str = None,
):
    """
    Simpe CLI wrapper to seqkit rmdup
    """
    if output_fasta is None:
        output_fasta = set_default_output_path(input_fasta, tag="_no_duplicates")
    if export_duplicates:
        if duplicates_file is None:
            dup_file = set_default_output_path(
                input_fasta, tag="_duplicates", extension=".txt"
            )
        else:
            dup_file = duplicates_file
        dup_str = f"-D {dup_file}"
    else:
        dup_str = ""
    cmd_str = f"seqkit rmdup {input_fasta} --quiet -s {dup_str} -o {output_fasta}"
    terminal_execute(cmd_str, suppress_shell_output=True)


def run_fastp(
    input_fastq1: str,
    input_fastq2: str = None,
    merge: bool = False,
    merged_out: str = None,
    output_dir: str = None,
    n_threads: int = None,
    additional_args: str = None,
) -> None:
    """
    Simple CLI wrapper to Fastp
    https://github.com/OpenGene/fastp
    """
    if input_fastq2 is None:
        input_str = f"-i {input_fastq1}"
    else:
        input_str = f"-i {input_fastq1} -I {input_fastq2}"
    if n_threads is not None:
        threads_str = f"--thread {n_threads}"
    else:
        threads_str = ""
    if output_dir is None:
        output_dir = set_default_output_path(input_fastq1, only_dirname=True)
    if merge:
        merged_out = (
            merged_out
            if merged_out is not None
            else os.path.join(output_dir, "merged.fastq")
        )
        out_str = f"--merge --merged_out {merged_out}"
    else:
        out_str = f'-o {"fastp_" + input_fastq1} -O {"fastp_" + input_fastq2}'
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ""
    cmd_str = f"fastp {input_str} {out_str} {threads_str} {args_str}"
    terminal_execute(cmd_str, suppress_shell_output=True)


def run_prodigal(
    input_file: str,
    output_prefix: str = None,
    output_dir: str = None,
    metagenome: bool = False,
    additional_args: str = None,
):
    """
    Simple CLI wrapper to prodigal
    """
    if metagenome:
        procedure = "meta"
    else:
        procedure = "single"
    if output_dir is None:
        output_dir = set_default_output_path(input_file, only_dirname=True)
    if output_prefix is None:
        output_prefix = set_default_output_path(input_file, only_filename=True)
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ""
    output_gbk = os.path.join(output_dir, output_prefix + ".gbk")
    output_fasta = os.path.join(output_dir, output_prefix + ".faa")
    cmd_str = (
        f"prodigal -i {input_file} -o {output_gbk} -p {procedure} "
        f"-a {output_fasta} -q {args_str}"
    )
    terminal_execute(cmd_str, suppress_shell_output=True)


def run_hmmsearch(
    hmm_model: str,
    input_fasta: str,
    output_file: str = None,
    method: str = "hmmsearch",
    n_processes: int = None,
    additional_args: str = None,
) -> None:
    """
    Simple CLI wrapper to hmmsearch or hmmscan
    """
    if n_processes is None:
        n_processes = os.cpu_count() - 1
    if output_file is None:
        output_file = set_default_output_path(input_fasta, "_hmmer_hits", ".txt")
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ""
    cmd_str = (
        f"{method} --tblout {output_file} {args_str} --cpu {n_processes} "
        f"{hmm_model} {input_fasta}"
    )
    terminal_execute(cmd_str, suppress_shell_output=True)


def run_hmmbuild(
    input_aln: str, output_hmm: str = None, additional_args: str = None
) -> None:
    """
    Simple CLI wrapper to hmmbuild (build HMM profile from MSA file)
    additional args: see hmmbuild -h
    """
    if output_hmm is None:
        output_hmm = set_default_output_path(input_aln, extension=".hmm")
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ""
    cmd_str = f"hmmbuild {args_str} {output_hmm} {input_aln}"
    terminal_execute(cmd_str, suppress_shell_output=True)


def run_hmmalign(
    input_hmm: str,
    input_aln: str,
    input_seqs: str,
    output_aln_seqs: str = None,
    additional_args: str = None,
) -> None:
    """
    Simple CLI wrapper to hmmalign
    Align short read query sequences to reference MSA
    """
    if output_aln_seqs is None:
        output_aln_seqs = set_default_output_path(input_seqs, "_hmm", extension=".aln")
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ""
    input_hmm = os.path.abspath(input_hmm)
    input_aln = os.path.abspath(input_aln)
    input_seqs = os.path.abspath(input_seqs)
    cmd_str = (
        f"hmmalign -o {output_aln_seqs} --mapali {input_aln} --trim "
        f"--informat fasta {args_str} {input_hmm} {input_seqs}"
    )
    terminal_execute(cmd_str, suppress_shell_output=True)


def get_percent_identity_from_msa(input_msa: str, output_file: str = None) -> None:
    """
    Run esl-alipid to compute pairwise PI from a MSA.
    """
    input_msa = os.path.abspath(input_msa)
    if output_file is None:
        output_file = set_default_output_path(input_msa, tag="_PI", extension=".txt")
    cmd_str = f"esl-alipid {input_msa} > {output_file}"
    terminal_execute(cmd_str, suppress_shell_output=False)


def run_cdhit(
    input_fasta: str, output_fasta: str = None, additional_args: str = None
) -> None:
    """
    Simple CLI wrapper to cd-hit to obtain representative sequences
    CD-HIT may be used to remove duplicated sequences (keeps one representative)
    with parameters -c 1 -t 1.
    """
    if output_fasta is None:
        output_fasta = set_default_output_path(input_fasta, "_cdhit")
    if additional_args is None:
        additional_args = ""
    input_fasta = os.path.abspath(input_fasta)
    cmd_str = f"cd-hit -i {input_fasta} -o {output_fasta} {additional_args}"
    terminal_execute(cmd_str, suppress_shell_output=True)


def run_mafft(
    input_fasta: str,
    output_file: str = None,
    n_threads: int = -1,
    parallel: bool = True,
    additional_args: str = None,
) -> None:
    """
    Simple CLI wrapper to mafft (MSA)

    Manual: https://mafft.cbrc.jp/alignment/software/manual/manual.html

    CLI examples:
    mafft --globalpair --thread n in > out
    mafft --localpair --thread n in > out
    mafft --large --globalpair --thread n in > out
    """
    if output_file is None:
        output_file = set_default_output_path(
            input_fasta, extension=".fasta.aln", only_filename=True
        )
    if parallel:
        thread_str = f"--thread {n_threads}"
    else:
        thread_str = ""
    if additional_args is None:
        additional_args = ""
    input_fasta = os.path.abspath(input_fasta)
    cmd_str = f"mafft {thread_str} {additional_args} {input_fasta} > {output_file}"
    terminal_execute(cmd_str, suppress_shell_output=False)


def run_muscle(
    input_fasta: str,
    output_file: str = None,
    maxiters: int = None,
    additional_args: str = None,
) -> None:
    """
    Simple CLI wrapper to muscle (MSA)
    muscle: https://www.drive5.com/muscle/manual/output_formats.html

    output phylip and fasta.aln
    """
    if output_file is None:
        output_file = set_default_output_path(
            input_fasta, extension=".fasta.aln", only_filename=True
        )
    if maxiters is None:
        maxiters = 2
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ""
    input_fasta = os.path.abspath(input_fasta)
    cmd_str = (
        f"muscle -in {input_fasta} -out {output_file} "
        f"-maxiters {maxiters} {args_str}"
    )
    terminal_execute(cmd_str, suppress_shell_output=True)


def run_trimal(input_aln: str, output_aln: str = None) -> None:
    """
    Simple CLI wrapper to trimal

    I/O in phylip as well: https://vicfero.github.io/trimal/
    """
    if output_aln is None:
        output_aln = set_default_output_path(input_aln, "_trimal")
    input_aln = os.path.abspath(input_aln)
    cmd_str = (
        f"trimal -in {input_aln} -out {output_aln} -fasta -automated1 "
        f"-resoverlap 0.55 -seqoverlap 60 -htmlout trimal.html"
    )
    terminal_execute(cmd_str, suppress_shell_output=True)


def run_modeltest(
    input_algns: str, n_processes: int = None, output_dir: str = None
) -> None:
    """
    Simple CLI wrapper to modeltest-ng
    Repo: https://github.com/ddarriba/modeltest
    """
    if output_dir is None:
        output_dir = os.path.join(
            set_default_output_path(input_algns, only_dirname=True), "modeltest_result"
        )
    else:
        output_dir = os.path.abspath(os.path.join(output_dir, "modeltest_result"))
    if n_processes is None:
        n_processes = os.cpu_count() - 1
    cmd_str = (
        f"modeltest-ng -i {input_algns} -T raxml -d aa -p {n_processes} "
        f"-o {output_dir}"
    )
    terminal_execute(cmd_str, suppress_shell_output=True)


def run_fasttree(
    input_algns: str,
    output_file: str = None,
    nucleotides: bool = False,
    starting_tree: str = None,
    additional_args: str = None,
) -> None:
    """
    Simple CLI wrapper to fasttree.
    fasttree accepts multiple alignments in fasta or phylip formats
    It seems that fasttree does not allow inputing subsitution model.
    Default substitution model for protein seqs is JTT

    additional_args: a string containing additional parameters and
                     parameter values to be passed to fasttree
    """
    if output_file is None:
        output_file = set_default_output_path(
            input_algns, tag="_fasttree", extension=".newick"
        )
    if nucleotides:
        nt_str = "-gtr -nt"
    else:
        nt_str = ""
    if starting_tree is not None:
        start_t_str = f"-intree {starting_tree}"
    else:
        start_t_str = ""
    if additional_args is None:
        additional_args = ""
    input_algns = os.path.abspath(input_algns)
    cmd_str = f"fasttree {nt_str} {input_algns} {start_t_str} {additional_args} > {output_file}"
    terminal_execute(cmd_str, suppress_shell_output=False)


def run_iqtree(
    input_algns: str,
    output_dir: str = None,
    output_prefix: str = None,
    keep_recovery_files: bool = False,
    nucleotides: bool = False,
    n_processes: int = None,
    substitution_model: str = "TEST",
    starting_tree: str = None,
    bootstrap_replicates: int = 1000,
    max_bootstrap_iterations: int = 1000,
    overwrite_previous_results: bool = True,
    additional_args: str = None,
) -> None:
    """
    Simple CLI wrapper to iqtree.
    iqtree accepts multiple alignments in fasta or phylip formats.

    additional_args: a string containing additional parameters and
                     parameter values to be passed to iqtree

    output: iqtree outputs several files

    Reducing computational time via model selection:
    http://www.iqtree.org/doc/Command-Reference
    """

    def remove_auxiliary_output(output_prefix):
        """
        Removes iqtree auxiliary output files
        """
        exts_to_remove = [
            ".bionj",
            ".ckp.gz",
            ".iqtree",
            ".log",
            ".model.gz",
            ".splits.nex",
            ".mldist",
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
        output_prefix_str = f"-pre {output_prefix}"
    else:
        output_prefix = os.path.join(output_dir, output_prefix)
        output_prefix_str = f"-pre {output_prefix}"
    if nucleotides:
        seq_type = "DNA"
    else:
        seq_type = "AA"
    if n_processes is None:
        n_processes = "AUTO"
    if starting_tree is not None:
        start_t_str = f"-t {starting_tree}"
    else:
        start_t_str = ""
    if overwrite_previous_results:
        overwrite_str = "-redo"
    else:
        overwrite_str = ""
    if additional_args is None:
        additional_args = ""

    input_algns = os.path.abspath(input_algns)
    cmd_str = (
        f"iqtree -s {input_algns} -st {seq_type} -nt {n_processes} "
        f"-m {substitution_model} -bb {bootstrap_replicates} -mset raxml "
        f"-nm {max_bootstrap_iterations} "
        f"{output_prefix_str} {start_t_str} {overwrite_str} {additional_args}"
    )
    terminal_execute(cmd_str, suppress_shell_output=True)
    if not keep_recovery_files:
        remove_auxiliary_output(output_prefix)


def run_tree_shrink(
    input_tree: str,
    input_aln: str = None,
    output_dir: str = None,
    output_deleted_nodes: bool = False,
    additional_args: str = None,
) -> None:
    """
    Run treeshrink to remove tree branch outliers.
    Remove outliers from MSA file too.
    Tree file must be of newick format.
    See run_treeshrink.py  -h for help
    Currently using default parameter values.
    """
    if output_dir is None:
        output_dir = os.path.dirname(input_tree)
    else:
        output_dir = os.path.abspath(output_dir)
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ""
    input_tree = os.path.abspath(input_tree)
    out_tree = set_default_output_path(input_tree, tag="_shrink", only_filename=True)
    if input_aln is not None:
        input_aln = os.path.abspath(input_aln)
        out_aln = set_default_output_path(input_aln, tag="_shrink", only_filename=True)
        aln_str = "-a input.aln"
    else:
        aln_str = ""

    # Handle treeshrink input/output requirements (temp/tree/input.tree)
    with tempfile.TemporaryDirectory() as temp_out_dir, tempfile.TemporaryDirectory() as parent_in_temp, tempfile.TemporaryDirectory(
        dir=parent_in_temp
    ) as temp_in_dir:
        temp_tree_dir = os.path.basename(temp_in_dir)
        shutil.copy(input_tree, os.path.join(temp_in_dir, "input.tree"))
        if input_aln is not None:
            shutil.copy(input_aln, os.path.join(temp_in_dir, "input.aln"))

        cmd_str = (
            f"run_treeshrink.py -i {parent_in_temp} -m per-gene "
            f"-t input.tree {aln_str} --force "
            f"-o {temp_out_dir} -O output {args_str}"
        )
        terminal_execute(cmd_str, suppress_shell_output=True)

        shutil.move(
            os.path.join(temp_out_dir, temp_tree_dir, "output.tree"),
            os.path.join(output_dir, out_tree),
        )
        if input_aln is not None:
            shutil.move(
                os.path.join(temp_out_dir, temp_tree_dir, "output.aln"),
                os.path.join(output_dir, out_aln),
            )
        if output_deleted_nodes:
            out_txt = set_default_output_path(
                input_tree, tag="_shrink_deleted", extension=".txt", only_filename=True
            )
            shutil.move(
                os.path.join(temp_out_dir, temp_tree_dir, "output.txt"),
                os.path.join(output_dir, out_txt),
            )


def run_papara(
    tree_nwk: str,
    msa_phy: str,
    query_fasta: str,
    output_aln: str = None,
    additional_args: str = None,
) -> None:
    """
    Simple CLI wrapper to Papara. Output lignment in fasta format

    Run Papara to do query alignment to reference MSA and tree (required for EPA-ng)
    Alignment could be done with hmmalign or muscle as well, but these tools don't
    consider the tree during alignment (would this be a justified improvement over hmmalign?)

    There seems to be a problem with enabling multithreading in papara when run as a static
    executable. It looks like it has to be enabled during compilation (but compilation currently not working):
    https://stackoverflow.com/questions/19618926/thread-doesnt-work-with-an-error-enable-multithreading-to-use-stdthread-ope
    """
    tree_nwk = os.path.abspath(tree_nwk)
    msa_phy = os.path.abspath(msa_phy)
    query_fasta = os.path.abspath(query_fasta)

    if output_aln is None:
        output_aln = set_default_output_path(
            query_fasta, tag="_ref_aln", extension=".phylip"
        )
    else:
        output_aln = os.path.abspath(output_aln)
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ""
    cmd_str = (
        f"{papara_exec} -t {tree_nwk} "
        f"-s {msa_phy} -q {query_fasta} -n phylip "
        f"-r {args_str} -a"
    )
    with tempfile.TemporaryDirectory() as tempdir:
        terminal_execute(cmd_str, suppress_shell_output=True, work_dir=tempdir)
        shutil.move(os.path.join(tempdir, "papara_alignment.phylip"), output_aln)


def run_epang(
    input_tree: str,
    input_aln_ref: str,
    input_aln_query: str,
    model: str = None,
    output_dir: str = None,
    n_threads: int = None,
    overwrite_previous_results: bool = True,
    additional_args: str = None,
) -> None:
    """
    Simple CLI wrapper to EPA-ng
    See epa-ng -h for additional parameters
    input_tree: newick format
    input_aln: fasta format
    input_aln_query: fasta format (sequences must be alignned to reference
    msa fasta and have the same length as the reference msa alignment)
    epa-ng: https://github.com/Pbdas/epa-ng
    """
    if model is None:
        model = "GTR+G"
    if output_dir is None:
        output_dir = os.path.dirname(input_tree)
    else:
        output_dir = os.path.abspath(output_dir)
    if n_threads is None:
        n_threads = os.cpu_count() - 1
    if overwrite_previous_results:
        overwrite_str = "--redo"
    else:
        overwrite_str = ""
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ""

    cmd_str = (
        f"epa-ng --ref-msa {input_aln_ref} --tree {input_tree} --query {input_aln_query} "
        f"--model {model} --threads {n_threads} --outdir {output_dir} {overwrite_str} {args_str}"
    )
    terminal_execute(cmd_str, suppress_shell_output=True)


def run_gappa_heat_tree(
    input_jplace: str,
    output_dir: str = None,
    output_prefix: str = None,
    additional_args: str = None,
) -> None:
    """
    Run gappa examine heat-tree to obtain heat-map tree
    representing densities of placed short reads in svg
    or pdf formats
    """
    if output_dir is None:
        outdir_str = ""
    else:
        outdir_str = f"--out-dir {os.path.abspath(output_dir)}"
    if output_prefix is None:
        output_prefix_str = ""
    else:
        output_prefix_str = f"--file-prefix {output_prefix}"
    if additional_args is None:
        args_str = ""
    else:
        args_str = additional_args

    cmd_str = (
        f"gappa examine heat-tree --jplace-path {os.path.abspath(input_jplace)} "
        f"--write-svg-tree "  # --write-newick-tree
        f"{outdir_str} {output_prefix_str} {args_str}"
    )
    terminal_execute(cmd_str, suppress_shell_output=True)


def run_gappa_graft(
    input_jplace: str,
    output_dir: str = None,
    output_prefix: str = None,
    additional_args: str = None,
) -> None:
    """
    Run gappa examine graft to obtain tree with placements in
    newick format
    """
    if output_dir is None:
        outdir_str = ""
    else:
        outdir_str = f"--out-dir {os.path.abspath(output_dir)}"
    if output_prefix is None:
        output_prefix_str = ""
    else:
        output_prefix_str = f"--file-prefix {output_prefix}"
    if additional_args is None:
        args_str = ""
    else:
        args_str = additional_args

    cmd_str = (
        f"gappa examine graft --jplace-path {os.path.abspath(input_jplace)} "
        f"{outdir_str} {output_prefix_str} {args_str}"
    )
    terminal_execute(cmd_str, suppress_shell_output=True)


def run_gappa_assign(
    jplace: str,
    taxonomy_file: str,
    output_dir: str = None,
    output_prefix: str = None,
    only_best_hit: bool = True,
    additional_args: str = None,
    delete_output_tree: bool = True,
) -> None:
    """
    Use gappa examine assign to assign taxonomy to placed query sequences
    based on taxonomy assigned to tree reference sequences

    argument: --resolve-missing-paths alongside --root-outgroup can be
    added to find missing taxonomic info in labels.

    Info: https://github.com/lczech/gappa/wiki/Subcommand:-assign
    """
    if output_dir is None:
        output_dir = set_default_output_path(jplace, only_dirname=True)
    if output_prefix is None:
        output_prefix = set_default_output_path(jplace, only_filename=True)
    if only_best_hit:
        best_str = "--best-hit"
    else:
        best_str = ""
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ""

    cmd_str = (
        f"gappa examine assign --jplace-path {jplace} --taxon-file {taxonomy_file} "
        f"--out-dir {output_dir} --file-prefix {output_prefix} --allow-file-overwriting "
        f"--per-query-results {best_str} {args_str}"
    )
    terminal_execute(cmd_str, suppress_shell_output=True)
    outtree = os.path.join(output_dir, f"{output_prefix}labelled_tree.newick")
    if delete_output_tree and os.path.exists(outtree):
        os.remove(outtree)
