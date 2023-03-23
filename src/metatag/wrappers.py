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
    input_fasta: Path,
    output_fasta: Path = None,
    export_duplicates: bool = False,
    duplicates_file: Path = None,
):
    """Simpe CLI wrapper to seqkit rmdup

    Args:
        input_fasta (Path): _description_
        output_fasta (Path, optional): _description_. Defaults to None.
        export_duplicates (bool, optional): _description_. Defaults to False.
        duplicates_file (Path, optional): _description_. Defaults to None.
    """
    if output_fasta is None:
        output_fasta = set_default_output_path(input_fasta, tag="_no_duplicates")
    else:
        output_fasta = Path(output_fasta).resolve()
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


def run_prodigal(
    input_file: Path,
    output_prefix: str = None,
    output_dir: Path = None,
    metagenome: bool = False,
    additional_args: str = None,
):
    """Simple CLI wrapper to prodigal

    Args:
        input_file (Path): _description_
        output_prefix (str, optional): _description_. Defaults to None.
        output_dir (Path, optional): _description_. Defaults to None.
        metagenome (bool, optional): _description_. Defaults to False.
        additional_args (str, optional): _description_. Defaults to None.
    """
    input_file = Path(input_file).resolve()
    if output_dir is None:
        output_dir = set_default_output_path(input_file, only_dirname=True)
    else:
        output_dir = Path(output_dir).resolve()
    if output_prefix is None:
        output_prefix = set_default_output_path(input_file, only_filename=True)
    if metagenome:
        procedure = "meta"
    else:
        procedure = "single"
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ""
    output_gbk = output_dir / f"{output_prefix}.gbk"
    output_fasta = output_dir / f"{output_prefix}.faa"
    cmd_str = (
        f"prodigal -i {input_file} -o {output_gbk} -p {procedure} "
        f"-a {output_fasta} -q {args_str}"
    )
    terminal_execute(cmd_str, suppress_shell_output=True)


def run_hmmsearch(
    hmm_model: Path,
    input_fasta: Path,
    output_file: Path = None,
    method: str = "hmmsearch",
    n_processes: int = None,
    additional_args: str = None,
) -> None:
    """Simple CLI wrapper to hmmsearch or hmmscan

    Args:
        hmm_model (Path): _description_
        input_fasta (Path): _description_
        output_file (Path, optional): _description_. Defaults to None.
        method (str, optional): _description_. Defaults to "hmmsearch".
        n_processes (int, optional): _description_. Defaults to None.
        additional_args (str, optional): _description_. Defaults to None.
    """
    hmm_model = Path(hmm_model).resolve()
    input_fasta = Path(input_fasta).resolve()
    if output_file is None:
        output_file = set_default_output_path(input_fasta, "_hmmer_hits", ".txt")
    else:
        output_file = Path(output_file).resolve()
    if n_processes is None:
        n_processes = os.cpu_count() - 1
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
    input_aln: Path, output_hmm: Path = None, additional_args: str = None
) -> None:
    """Simple CLI wrapper to hmmbuild (build HMM profile from MSA file)
    additional args: see hmmbuild -h

    Args:
        input_aln (Path): _description_
        output_hmm (Path, optional): _description_. Defaults to None.
        additional_args (str, optional): _description_. Defaults to None.
    """
    input_aln = Path(input_aln).resolve()
    if output_hmm is None:
        output_hmm = set_default_output_path(input_aln, extension=".hmm")
    else:
        output_hmm = Path(output_hmm).resolve()
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ""
    cmd_str = f"hmmbuild {args_str} {output_hmm} {input_aln}"
    terminal_execute(cmd_str, suppress_shell_output=True)


def run_hmmalign(
    input_hmm: Path,
    input_aln: Path,
    input_seqs: Path,
    output_aln_seqs: Path = None,
    additional_args: str = None,
) -> None:
    """Simple CLI wrapper to hmmalign
    Align short read query sequences to reference MSA

    Args:
        input_hmm (Path): _description_
        input_aln (Path): _description_
        input_seqs (Path): _description_
        output_aln_seqs (Path, optional): _description_. Defaults to None.
        additional_args (str, optional): _description_. Defaults to None.
    """
    input_hmm = Path(input_hmm).resolve()
    input_aln = Path(input_aln).resolve()
    input_seqs = Path(input_seqs).resolve()
    if output_aln_seqs is None:
        output_aln_seqs = set_default_output_path(input_seqs, "_hmm", extension=".aln")
    else:
        output_aln_seqs = Path(output_aln_seqs).resolve()
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ""
    cmd_str = (
        f"hmmalign -o {output_aln_seqs} --mapali {input_aln} --trim "
        f"--informat fasta {args_str} {input_hmm} {input_seqs}"
    )
    terminal_execute(cmd_str, suppress_shell_output=True)


def get_percent_identity_from_msa(input_msa: Path, output_file: Path = None) -> None:
    """Run esl-alipid to compute pairwise PI from a MSA.

    Args:
        input_msa (Path): _description_
        output_file (Path, optional): _description_. Defaults to None.
    """
    input_msa = Path(input_msa).resolve()
    if output_file is None:
        output_file = set_default_output_path(input_msa, tag="_PI", extension=".txt")
    else:
        output_file = Path(output_file).resolve()
    cmd_str = f"esl-alipid {input_msa} > {output_file}"
    terminal_execute(cmd_str, suppress_shell_output=False)


def run_cdhit(
    input_fasta: Path, output_fasta: Path = None, additional_args: str = None
) -> None:
    """Simple CLI wrapper to cd-hit to obtain representative sequences
    CD-HIT may be used to remove duplicated sequences (keeps one representative)
    with parameters -c 1 -t 1.

    Args:
        input_fasta (Path): _description_
        output_fasta (Path, optional): _description_. Defaults to None.
        additional_args (str, optional): _description_. Defaults to None.
    """
    input_fasta = Path(input_fasta).resolve()
    if output_fasta is None:
        output_fasta = set_default_output_path(input_fasta, "_cdhit")
    else:
        output_fasta = Path(output_fasta).resolve()
    if additional_args is None:
        additional_args = ""
    cmd_str = f"cd-hit -i {input_fasta} -o {output_fasta} {additional_args}"
    terminal_execute(cmd_str, suppress_shell_output=True)


def run_mafft(
    input_fasta: Path,
    output_file: Path = None,
    n_threads: int = -1,
    parallel: bool = True,
    additional_args: str = None,
) -> None:
    """Simple CLI wrapper to mafft (MSA)

    Manual: https://mafft.cbrc.jp/alignment/software/manual/manual.html

    CLI examples:
    mafft --globalpair --thread n in > out
    mafft --localpair --thread n in > out
    mafft --large --globalpair --thread n in > out

    Args:
        input_fasta (Path): _description_
        output_file (Path, optional): _description_. Defaults to None.
        n_threads (int, optional): _description_. Defaults to -1.
        parallel (bool, optional): _description_. Defaults to True.
        additional_args (str, optional): _description_. Defaults to None.
    """
    input_fasta = Path(input_fasta).resolve()
    if output_file is None:
        output_file = set_default_output_path(
            input_fasta, extension=".fasta.aln", only_filename=True
        )
    else:
        output_file = Path(output_file).resolve()
    if parallel:
        thread_str = f"--thread {n_threads}"
    else:
        thread_str = ""
    if additional_args is None:
        additional_args = ""
    cmd_str = f"mafft {thread_str} {additional_args} {input_fasta} > {output_file}"
    terminal_execute(cmd_str, suppress_shell_output=False)


def run_muscle(
    input_fasta: Path,
    output_file: Path = None,
    maxiters: int = None,
    additional_args: str = None,
) -> None:
    """Simple CLI wrapper to muscle (MSA)
    muscle: https://www.drive5.com/muscle/manual/output_formats.html

    output phylip and fasta.aln

    Args:
        input_fasta (Path): _description_
        output_file (Path, optional): _description_. Defaults to None.
        maxiters (int, optional): _description_. Defaults to None.
        additional_args (str, optional): _description_. Defaults to None.
    """
    input_fasta = Path(input_fasta).resolve()
    if output_file is None:
        output_file = set_default_output_path(
            input_fasta, extension=".fasta.aln", only_filename=True
        )
    else:
        output_file = Path(output_file).resolve()
    if maxiters is None:
        maxiters = 2
    if additional_args is not None:
        args_str = additional_args
    else:
        args_str = ""
    cmd_str = (
        f"muscle -in {input_fasta} -out {output_file} "
        f"-maxiters {maxiters} {args_str}"
    )
    terminal_execute(cmd_str, suppress_shell_output=False)


def run_fasttree(
    input_algns: Path,
    output_file: Path = None,
    nucleotides: bool = False,
    starting_tree: str = None,
    additional_args: str = None,
) -> None:
    """Simple CLI wrapper to fasttree.
    fasttree accepts multiple alignments in fasta or phylip formats
    It seems that fasttree does not allow inputing subsitution model.
    Default substitution model for protein seqs is JTT

    additional_args: a string containing additional parameters and
                     parameter values to be passed to fasttree

    Args:
        input_algns (Path): _description_
        output_file (Path, optional): _description_. Defaults to None.
        nucleotides (bool, optional): _description_. Defaults to False.
        starting_tree (str, optional): _description_. Defaults to None.
        additional_args (str, optional): _description_. Defaults to None.
    """
    input_algns = Path(input_algns).resolve()
    if output_file is None:
        output_file = set_default_output_path(
            input_algns, tag="_fasttree", extension=".newick"
        )
    else:
        output_file = Path(output_file).resolve()
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
    cmd_str = f"fasttree {nt_str} {input_algns} {start_t_str} {additional_args} > {output_file}"
    terminal_execute(cmd_str, suppress_shell_output=False)


def remove_auxiliary_output(output_prefix: str):
    """Removes iqtree auxiliary output files

    Args:
        output_prefix (str): _description_
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
    files_to_remove = [Path(output_prefix + ext) for ext in exts_to_remove]
    for file_path in files_to_remove:
        file_path.unlink(missing_ok=True)


def run_iqtree(
    input_algns: Path,
    output_dir: Path = None,
    output_prefix: str = None,
    keep_recovery_files: bool = False,
    nucleotides: bool = False,
    n_processes: int = None,
    substitution_model: str = "TEST",
    starting_tree: Path = None,
    bootstrap_replicates: int = 1000,
    max_bootstrap_iterations: int = 1000,
    overwrite_previous_results: bool = True,
    additional_args: str = None,
) -> None:
    """Simple CLI wrapper to iqtree.
    iqtree accepts multiple alignments in fasta or phylip formats.

    additional_args: a string containing additional parameters and
                     parameter values to be passed to iqtree

    output: iqtree outputs several files

    Reducing computational time via model selection:
    http://www.iqtree.org/doc/Command-Reference

    Args:
        input_algns (Path): _description_
        output_dir (Path, optional): _description_. Defaults to None.
        output_prefix (str, optional): _description_. Defaults to None.
        keep_recovery_files (bool, optional): _description_. Defaults to False.
        nucleotides (bool, optional): _description_. Defaults to False.
        n_processes (int, optional): _description_. Defaults to None.
        substitution_model (str, optional): _description_. Defaults to "TEST".
        starting_tree (Path, optional): _description_. Defaults to None.
        bootstrap_replicates (int, optional): _description_. Defaults to 1000.
        max_bootstrap_iterations (int, optional): _description_. Defaults to 1000.
        overwrite_previous_results (bool, optional): _description_. Defaults to True.
        additional_args (str, optional): _description_. Defaults to None.
    """
    input_algns = Path(input_algns).resolve()
    if output_dir is None:
        output_dir = input_algns.parent
    else:
        output_dir = Path(output_dir).resolve()
    if output_prefix is None:
        input_file = input_algns.name
        output_prefix = (output_dir / input_file).as_posix()
        output_prefix_str = f"-pre {output_prefix}"
    else:
        output_prefix = (output_dir / output_prefix).as_posix()
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

    cmd_str = (
        f"iqtree -s {input_algns} -st {seq_type} -nt {n_processes} "
        f"-m {substitution_model} -bb {bootstrap_replicates} -mset raxml "
        f"-nm {max_bootstrap_iterations} "
        f"{output_prefix_str} {start_t_str} {overwrite_str} {additional_args}"
    )
    terminal_execute(cmd_str, suppress_shell_output=True)
    if not keep_recovery_files:
        remove_auxiliary_output(output_prefix)


def run_papara(
    tree_nwk: Path,
    msa_phy: Path,
    query_fasta: Path,
    output_aln: Path = None,
    additional_args: str = None,
) -> None:
    """Simple CLI wrapper to Papara. Output lignment in fasta format

    Run Papara to do query alignment to reference MSA and tree (required for EPA-ng)
    Alignment could be done with hmmalign or muscle as well, but these tools don't
    consider the tree during alignment (would this be a justified improvement over hmmalign?)

    There seems to be a problem with enabling multithreading in papara when run as a static
    executable. It looks like it has to be enabled during compilation (but compilation currently not working):
    https://stackoverflow.com/questions/19618926/thread-doesnt-work-with-an-error-enable-multithreading-to-use-stdthread-ope

    Args:
        tree_nwk (Path): _description_
        msa_phy (Path): _description_
        query_fasta (Path): _description_
        output_aln (Path, optional): _description_. Defaults to None.
        additional_args (str, optional): _description_. Defaults to None.
    """
    tree_nwk = Path(tree_nwk).resolve()
    msa_phy = Path(msa_phy).resolve()
    query_fasta = Path(query_fasta).resolve()

    if output_aln is None:
        output_aln = set_default_output_path(
            query_fasta, tag="_ref_aln", extension=".phylip"
        )
    else:
        output_aln = Path(output_aln).resolve()
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
    input_tree: Path,
    input_aln_ref: Path,
    input_aln_query: Path,
    model: str = None,
    output_dir: Path = None,
    n_threads: int = None,
    overwrite_previous_results: bool = True,
    additional_args: str = None,
) -> None:
    """Simple CLI wrapper to EPA-ng
    See epa-ng -h for additional parameters
    input_tree: newick format
    input_aln: fasta format
    input_aln_query: fasta format (sequences must be alignned to reference
    msa fasta and have the same length as the reference msa alignment)
    epa-ng: https://github.com/Pbdas/epa-ng

    Args:
        input_tree (Path): _description_
        input_aln_ref (Path): _description_
        input_aln_query (Path): _description_
        model (str, optional): _description_. Defaults to None.
        output_dir (Path, optional): _description_. Defaults to None.
        n_threads (int, optional): _description_. Defaults to None.
        overwrite_previous_results (bool, optional): _description_. Defaults to True.
        additional_args (str, optional): _description_. Defaults to None.
    """
    input_tree = Path(input_tree).resolve()
    input_aln_ref = Path(input_aln_ref).resolve()
    input_aln_query = Path(input_aln_query).resolve()
    if output_dir is None:
        output_dir = input_tree.parent
    else:
        output_dir = Path(output_dir).resolve()
    if model is None:
        model = "GTR+G"
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


def run_gappa_graft(
    input_jplace: Path,
    output_dir: Path = None,
    output_prefix: str = None,
    additional_args: str = None,
) -> None:
    """Run gappa examine graft to obtain tree with placements in
    newick format

    Args:
        input_jplace (Path): _description_
        output_dir (Path, optional): _description_. Defaults to None.
        output_prefix (str, optional): _description_. Defaults to None.
        additional_args (str, optional): _description_. Defaults to None.
    """
    input_jplace = Path(input_jplace).resolve()
    if output_dir is None:
        outdir_str = ""
    else:
        outdir_str = f"--out-dir {Path(output_dir)}"
    if output_prefix is None:
        output_prefix_str = ""
    else:
        output_prefix_str = f"--file-prefix {output_prefix}"
    if additional_args is None:
        args_str = ""
    else:
        args_str = additional_args

    cmd_str = (
        f"gappa examine graft --jplace-path {input_jplace} "
        f"{outdir_str} {output_prefix_str} {args_str}"
    )
    terminal_execute(cmd_str, suppress_shell_output=True)


def run_gappa_assign(
    jplace: Path,
    taxonomy_file: Path,
    output_dir: Path = None,
    output_prefix: str = None,
    only_best_hit: bool = True,
    additional_args: str = None,
    delete_output_tree: bool = True,
) -> None:
    """Use gappa examine assign to assign taxonomy to placed query sequences
    based on taxonomy assigned to tree reference sequences

    argument: --resolve-missing-paths alongside --root-outgroup can be
    added to find missing taxonomic info in labels.

    Info: https://github.com/lczech/gappa/wiki/Subcommand:-assign

    Args:
        jplace (Path): _description_
        taxonomy_file (Path): _description_
        output_dir (Path, optional): _description_. Defaults to None.
        output_prefix (str, optional): _description_. Defaults to None.
        only_best_hit (bool, optional): _description_. Defaults to True.
        additional_args (str, optional): _description_. Defaults to None.
        delete_output_tree (bool, optional): _description_. Defaults to True.
    """
    jplace = Path(jplace).resolve()
    taxonomy_file = Path(taxonomy_file).resolve()
    if output_dir is None:
        output_dir = set_default_output_path(jplace, only_dirname=True)
    else:
        output_dir = Path(output_dir).resolve()
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
    outtree = output_dir / f"{output_prefix}labelled_tree.newick"
    if delete_output_tree:
        outtree.unlink(missing_ok=True)
