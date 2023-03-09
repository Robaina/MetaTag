#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to quantify and assign labels to placed sequences
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
import re
import shutil
import sys
from io import StringIO
from pathlib import Path

import pandas as pd
from Bio import Phylo

import metatag.wrappers as wrappers
from metatag.alignment import align_short_reads_to_reference_msa
from metatag.database.manipulation import (
    get_fasta_record_ids,
    split_reference_from_query_alignments,
)
from metatag.phylotree import get_iq_tree_model_from_log_file
from metatag.taxonomy import TaxonomyAssigner, TaxonomyCounter, Taxopath
from metatag.utils import TemporaryFilePath, set_default_output_path

logger = logging.getLogger(__name__)

package_dir = Path(__file__).parent


class JplaceParser:
    """
    Methods to parse jplace files, as specified in
    https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009
    """

    def __init__(self, path_to_jplace: Path) -> None:
        self._path_to_jplace = Path(path_to_jplace)
        self.jplace = self.get_json_object()
        self._tree_obj = next(
            Phylo.parse(StringIO(self.newickfy_tree(self.jplace["tree"])), "newick")
        )

    def get_json_object(self) -> dict:
        with open(self._path_to_jplace, "r") as JSON:
            return json.load(JSON)

    @property
    def meta(self):
        """
        Print metadata
        """
        return self.jplace["metadata"]

    @property
    def fields(self):
        """
        Print data fields
        """
        return self.jplace["fields"]

    @property
    def placements(self):
        """
        Return placement objects
        """
        return self.jplace["placements"]

    def get_tree_str(self, newick=False) -> str:
        """
        Return tree string in original or newick format
        Original format contains branch labels
        in curly brackets. Newick format removes
        these labels.
        """
        if newick:
            return self.newickfy_tree(self.jplace["tree"])
        else:
            return self.jplace["tree"]

    @staticmethod
    def newickfy_tree(tree_str: str) -> str:
        """
        Remove branch IDs from jplace tree string
        """
        subs_tree = re.sub("\{(\d+)\}", "", tree_str)
        return subs_tree
        # return next(Phylo.parse(StringIO(subs_tree), 'newick'))

    def get_reference_sequences(self) -> list:
        """
        Get list of reference sequences in the placement tree
        """
        return [c.name for c in self._tree_obj.get_terminals()]

    def build_branch_dict(self) -> dict:
        """
        Build dictionary with edge/branch numbers as keys and
        reference tree leaves as values
        """

        def get_id(s):
            return int(re.search("\{(\d+)\}", s).group(1))

        original_tree = self.jplace["tree"]
        leaves = self._tree_obj.get_terminals()

        branches = {
            get_id(original_tree[original_tree.find(leaf.name) :]): leaf.name
            for leaf in leaves
        }
        return branches

    def extract_placement_fields(self, pfielddata: list) -> dict:
        """
        Get dict with placement field values from list of values
        """
        fields = self.jplace["fields"]
        return {field: pfielddata[i] for i, field in enumerate(fields)}

    def select_best_placement(self, placement_object: dict) -> dict:
        """
        Select placement with lowest likelihood
        """
        pdata = [
            self.extract_placement_fields(pfielddata)
            for pfielddata in placement_object["p"]
        ]
        lowest_like_placement = sorted(pdata, key=lambda x: x["likelihood"])[0]
        return {"p": lowest_like_placement, "n": placement_object["n"]}

    def select_best_placements(self):
        """
        Select placement with lowest likelihood for
        all placement objects in placements
        """
        best_placements = [
            self.select_best_placement(placement)
            for placement in self.jplace["placements"]
        ]
        return best_placements

    def compute_tree_diameter(self) -> float:
        """
        Find maximum (pairwise) distance between two tips
        (leaves) in the tree
        """
        root = self._tree_obj.root
        max_distance = 0.0
        tips = self._tree_obj.get_terminals()
        for tip in tips:
            self._tree_obj.root_with_outgroup(tip)
            new_max = max(self._tree_obj.depths().values())
            if new_max > max_distance:
                max_distance = new_max
        self._tree_obj.root_with_outgroup(root)
        return max_distance

    def filter_placements_by_minimum_lwr(
        self, minimum_lwr: float, output_file: Path = None
    ) -> None:
        """
        Filter placements by minimum LWR (from 0 to 1)
        """
        if output_file is None:
            output_file = set_default_output_path(self._path_to_jplace, tag=f"_min_lwr_{minimum_lwr}")
        else:
            output_file = Path(output_file)
        jplace = self.get_json_object()

        filtered_placement_objs = []
        for placement_object in jplace["placements"]:
            filtered_placements = []
            for placement in placement_object["p"]:
                edge_num, likelihood, lwr, distal_length, pendant_length = placement
                if lwr >= minimum_lwr:
                    filtered_placements.append(placement)
            if filtered_placements:
                placement_object["p"] = filtered_placements
                filtered_placement_objs.append(placement_object)
        jplace["placements"] = filtered_placement_objs
        with open(output_file, "w") as ofile:
            json.dump(jplace, ofile, indent=2)

    def filter_placements_by_max_pendant_to_tree_diameter_ratio(
        self, max_pendant_ratio: float, output_file: Path = None
    ) -> None:
        """
        Filter placements by maximum pendant length
        """
        if output_file is None:
            output_file = set_default_output_path(
                self._path_to_jplace, tag=f"_max_pendant_diameter_ratio_{max_pendant_ratio}"
                )
        else:
            output_file = Path(output_file)
        tree_diameter = self.compute_tree_diameter()
        logger.info(f"Filtering placements for tree diameter: {tree_diameter}")
        jplace = self.get_json_object()

        filtered_placement_objs = []
        for placement_object in jplace["placements"]:
            filtered_placements = []
            for placement in placement_object["p"]:
                edge_num, likelihood, lwr, distal_length, pendant_length = placement
                if pendant_length / tree_diameter <= max_pendant_ratio:
                    filtered_placements.append(placement)
            if filtered_placements:
                placement_object["p"] = filtered_placements
                filtered_placement_objs.append(placement_object)
        jplace["placements"] = filtered_placement_objs

        with open(output_file, "w") as ofile:
            json.dump(jplace, ofile, indent=2)

    def filter_placements_by_max_pendant_length(
        self, max_pendant_length: float, output_file: Path = None
    ) -> None:
        """
        Filter placements by maximum pendant length
        """
        if output_file is None:
            output_file = set_default_output_path(
                self._path_to_jplace, tag=f"_max_pendant_{max_pendant_length}"
                )
        else:
            output_file = Path(output_file)
        jplace = self.get_json_object()

        filtered_placement_objs = []
        for placement_object in jplace["placements"]:
            filtered_placements = []
            for placement in placement_object["p"]:
                edge_num, likelihood, lwr, distal_length, pendant_length = placement
                if pendant_length <= max_pendant_length:
                    filtered_placements.append(placement)
            if filtered_placements:
                placement_object["p"] = filtered_placements
                filtered_placement_objs.append(placement_object)
        jplace["placements"] = filtered_placement_objs

        with open(output_file, "w") as ofile:
            json.dump(jplace, ofile, indent=2)

    def filter_placements_by_max_pendant_to_distal_length_ratio(
        self, max_pendant_ratio: float, output_file: Path = None
    ) -> None:
        """
        Filter placements by maximum pendant length
        """
        if output_file is None:
            output_file = set_default_output_path(
                self._path_to_jplace, tag=f"_max_pendant_distal_ratio_{max_pendant_ratio}"
                )
        else:
            output_file = Path(output_file)
        jplace = self.get_json_object()

        filtered_placement_objs = []
        for placement_object in jplace["placements"]:
            filtered_placements = []
            for placement in placement_object["p"]:
                edge_num, likelihood, lwr, distal_length, pendant_length = placement
                if pendant_length / distal_length <= max_pendant_ratio:
                    filtered_placements.append(placement)
            if filtered_placements:
                placement_object["p"] = filtered_placements
                filtered_placement_objs.append(placement_object)
        jplace["placements"] = filtered_placement_objs

        with open(output_file, "w") as ofile:
            json.dump(jplace, ofile, indent=2)


class TaxAssignParser:
    """
    Parse function and taxonomy placement assignments table
    """

    def __init__(self, tax_assign_path: Path):
        self._tax_assign = pd.read_csv(tax_assign_path, sep="\t")
        self._tax_assign = self._tax_assign[self._tax_assign.cluster_id != "DISTANT"]
        self._tax_assign.cluster_score = self._tax_assign.cluster_score.apply(
            lambda x: float(x)
        )
        self._tax_assign.cluster_taxopath = self._tax_assign.cluster_taxopath.apply(
            lambda x: "Unspecified" if pd.isna(x) else x
        )
        self._tax_assign.taxopath = self._tax_assign.taxopath.apply(
            lambda x: "Unspecified" if pd.isna(x) else x
        )

    @property
    def taxlevels(self):
        return Taxopath().taxlevels

    def count_hits(
        self,
        cluster_ids: list[str] = None,
        score_threshold: float = None,
        taxopath_type: str = "taxopath",
        path_to_query_list: Path = None,
    ) -> TaxonomyCounter:
        """
        Count hits within given cluster ids and at specificied taxon level
        @Params:
        normalize=True, then results are reported as fractions of total
        cluster_ids: IDs of tree clusters to be included in the counting of placements
        score_threshold: global placement score threshold to filter
                         low-quality placements out
        taxopath_type: 'taxopath' to use gappa-assign taxopath or 'cluster_taxopath'
                        to use lowest common taxopath of the reference tree cluster
        path_to_query_list: Path, if not None, then a tsv is exported to defined location
                            containing those queries with correct cluster assignment (
                                according to defined 'valid' cluster ids or threshold
                            )
        """
        if cluster_ids is None:
            cluster_ids = self._tax_assign.cluster_id.unique()
        if score_threshold is None:
            score_threshold = 0.0

        query_hits = self._tax_assign[
            (
                (self._tax_assign.cluster_id.isin(cluster_ids))
                & (self._tax_assign.cluster_score >= score_threshold)
            )
        ]
        if path_to_query_list is not None:
            query_hits.to_csv(path_to_query_list, sep="\t", index=False)
        taxopath_hits = query_hits[taxopath_type].values
        if len(taxopath_hits) == 0:
            logger.error(
                "No placement hits returned for the provided filter parameters"
            )
            sys.exit(1)

        taxlevel_counter = TaxonomyCounter(taxopath_list=taxopath_hits)
        return taxlevel_counter


def place_reads_onto_tree(
    input_tree: Path,
    tree_model: str,
    ref_aln: Path,
    query_seqs: Path,
    aln_method: str = "papara",
    output_dir: Path = None,
) -> None:
    """
    Performs short read placement onto phylogenetic tree
    tree_model: str, either the model name or path to log output by iqtree
    workflow example: https://github.com/Pbdas/epa-ng/wiki/Full-Stack-Example
    """
    if output_dir is None:
        output_dir = set_default_output_path(query_seqs, only_dirname=True)
    else:
        output_dir = Path(output_dir)

    if Path(tree_model).is_file():
        tree_model = get_iq_tree_model_from_log_file(tree_model)
        logger.info(f"Running EPA-ng with inferred substitution model: {tree_model}")

    ref_ids = get_fasta_record_ids(ref_aln)
    ref_query_msa = output_dir / set_default_output_path(query_seqs, extension=".faln", only_filename=True)
    aln_ref_frac = set_default_output_path(ref_query_msa, tag="_ref_fraction", extension=".faln")
    aln_query_frac = set_default_output_path(ref_query_msa, tag="_query_fraction", extension=".faln")

    align_short_reads_to_reference_msa(
        ref_msa=ref_aln,
        query_seqs=query_seqs,
        method=aln_method,
        tree_nwk=input_tree,
        output_dir=output_dir,
    )

    split_reference_from_query_alignments(
        ref_query_msa=ref_query_msa, ref_ids=ref_ids, output_dir=output_dir
    )

    wrappers.run_epang(
        input_tree=input_tree,
        input_aln_ref=aln_ref_frac,
        input_aln_query=aln_query_frac,
        model=tree_model,
        output_dir=output_dir,
        n_threads=None,
        additional_args=None,
    )


def parse_tree_clusters(
    clusters_tsv: Path, cluster_as_key: bool = True
) -> dict:
    """
    Parse clusters text file into dictionary
    @param
    clusters_as_key: if True then dict keys are cluster IDs and values
    are lists of reference IDs. If False, dict keys are reference IDs and
    values the cluster to which they belong.
    """
    df = pd.read_csv(clusters_tsv, sep="\t", dtype=str)
    if cluster_as_key:
        cluster_ids = df.cluster.unique()
        return {
            cluster_id: df[df.cluster == cluster_id].id.values.tolist()
            for cluster_id in cluster_ids
        }
    else:
        return dict(df.values)


def parse_tree_cluster_quality_scores(cluster_scores_tsv: Path) -> dict:
    """
    Parse cluster quality scores file into dictionary
    """
    df = pd.read_csv(cluster_scores_tsv, sep="\t")
    df["cluster"] = df["cluster"].astype(str)
    return dict(df.values)


def add_clusters_to_tax_table(
    in_taxtable: Path, clusters: dict, out_taxtable: Path = None
) -> None:
    """
    Add tree cluster info at the beginning of each taxopath
    according to clusters defined in dictionary 'clusters'
    """
    if out_taxtable is None:
        out_taxtable = set_default_output_path(in_taxtable, tag="_clustered")
    else:
        out_taxtable = Path(out_taxtable)
    taxtable = pd.read_csv(in_taxtable, sep="\t", header=None, dtype=str)
    for i, row in taxtable.iterrows():
        row[1] = clusters[row[0]] + ";" + row[1]
    taxtable.to_csv(out_taxtable, sep="\t", index=False, header=None)


def parse_gappa_assign_table(
    input_table: Path,
    has_cluster_id: bool = True,
    cluster_scores: dict = None,
    clusters_taxopath: dict = None,
    output_file: Path = None,
) -> None:
    """
    Parse gappa assign per query taxonomy assignment result tsv
    has_cluster_id: set to True if results table includes reference
                    cluster info in the first element of taxopath
    cluster_scores: dictionary with values set to cluster quality
                    scores. It is only used if has_cluster_id = True.
    """
    if output_file is None:
        output_file = set_default_output_path(input_table, tag="_parsed")
    else:
        output_file = Path(output_file)
    if (has_cluster_id) and (cluster_scores is not None):
        cluster_str = (
            "cluster_id" + "\t" + "cluster_score" + "\t" + "cluster_taxopath" + "\t"
        )
    elif (has_cluster_id) and (cluster_scores is None):
        cluster_str = "cluster_id" + "\t" + "cluster_taxopath" + "\t"
    elif not has_cluster_id:
        cluster_str = ""
    table = pd.read_csv(input_table, sep="\t")
    with open(output_file, "w") as file:
        lines = []
        header = "query_id" + "\t" + "LWR" + "\t" + cluster_str + "taxopath" + "\n"
        lines.append(header)
        for i, row in table.iterrows():
            if ";" not in row.taxopath:
                row.taxopath = row.taxopath + ";Unspecified"
            if has_cluster_id:
                elems = row.taxopath.split(";")
                cluster_id = elems[0]
                if cluster_id in clusters_taxopath.keys():
                    cluster_taxopath = clusters_taxopath[cluster_id]
                else:
                    cluster_taxopath = ""
                if not cluster_taxopath:
                    cluster_taxopath = "Unspecified"
                taxopath = cluster_taxopath + "\t" + ";".join(elems[1:])
            else:
                cluster_id, taxopath = "", row.taxopath
            if cluster_scores is not None:
                if cluster_id in cluster_scores.keys():
                    cluster_score = str(cluster_scores[cluster_id])
                else:
                    cluster_score = str(None)
                line = (
                    row["name"]
                    + "\t"
                    + str(row["LWR"])
                    + "\t"
                    + cluster_id
                    + "\t"
                    + cluster_score
                    + "\t"
                    + taxopath
                    + "\n"
                )
            else:
                line = (
                    row["name"]
                    + "\t"
                    + str(row["LWR"])
                    + "\t"
                    + cluster_id
                    + "\t"
                    + taxopath
                    + "\n"
                )
            lines.append(line)
        file.writelines(lines)


def add_query_labels_to_assign_table(
    input_table: Path, query_labels: dict, output_table: Path = None
) -> None:
    """
    Add new column containing actual query labels to query taxonomy/cluster assignment
    table
    """

    def relabel_query(query_id: str, query_labels: dict) -> str:
        try:
            return query_labels[query_id]
        except Exception:
            return query_id

    if output_table is None:
        output_table = set_default_output_path(input_table, tag="_relabel")
    else:
        output_table = Path(output_table)
    df = pd.read_csv(input_table, sep="\t")
    df.insert(
        1, "query_name", df.query_id.apply(lambda x: relabel_query(x, query_labels))
    )
    df = df.set_index("query_id")
    df.to_csv(output_table, sep="\t")


def assign_labels_to_placements(
    jplace: Path,
    ref_labels: dict,
    query_labels: dict = None,
    output_dir: Path = None,
    output_prefix: str = None,
    only_best_hit: bool = True,
    ref_clusters_file: Path = None,
    ref_cluster_scores_file: Path = None,
    gappa_additional_args: str = None,
    only_unique_cluster: bool = True,
    taxo_file: Path = None,
) -> None:
    """
    Assign taxonomy and/or tree cluster IDs to placed query sequences based on
    taxonomy assigned to tree reference sequences using
    gappa examine assign.
    @parameter
    ref_labels: dictionary containing short IDs as keys and long labels as values
                for reference sequences
    query_labels: dictionary containing short IDs as keys and long labels as values
                  for query sequences
    only_best_hit: only report taxonomy with largest LWR
                   per query
    ref_clusters: dict (optionally) add tree cluster to the
                  beginning of the taxopath so query sequences
                  can also be classified according to tree
                  cluster (e.g., assigned function)
    only_unique_cluster: if True, keep only queries with multiple placement locations
                         if they were assigned to the same cluster.
    """
    if output_dir is None:
        output_dir = set_default_output_path(jplace, only_dirname=True)
    else:
        output_dir = Path(output_dir)
    if output_prefix is None:
        output_prefix = set_default_output_path(jplace, only_filename=True)

    if ref_clusters_file is not None:
        has_cluster_id = True
        ref_clusters = parse_tree_clusters(
            ref_clusters_file, cluster_as_key=False, sep="\t"
        )
        ref_clusters_as_keys = parse_tree_clusters(
            ref_clusters_file, cluster_as_key=True, sep="\t"
        )
    else:
        has_cluster_id = False

    if ref_cluster_scores_file is not None:
        ref_cluster_scores = parse_tree_cluster_quality_scores(ref_cluster_scores_file)
    else:
        ref_cluster_scores = None

    if taxo_file is None:
        taxo_file = package_dir / "data" / "merged_taxonomy.tsv"
    else:
        taxo_file = Path(taxo_file)

    gappa_assign_out = output_dir / f"{output_prefix}per_query.tsv"
    query_taxo_out = output_dir / f"{output_prefix}assignments.tsv"

    # Remove references that are not in placement tree
    ref_in_jplace_tree = JplaceParser(jplace).get_reference_sequences()
    ref_labels = {k: v for k, v in ref_labels.items() if k in ref_in_jplace_tree}

    with TemporaryFilePath() as temptax:
        taxonomy = TaxonomyAssigner(taxo_file=taxo_file)
        taxonomy.build_gappa_taxonomy_table(ref_labels, output_file=temptax)

        if has_cluster_id:
            add_clusters_to_tax_table(
                in_taxtable=temptax, clusters=ref_clusters, out_taxtable=temptax
            )
            clusters_taxopath = taxonomy.assign_lowest_common_taxonomy_to_clusters(
                clusters=ref_clusters_as_keys, label_dict=ref_labels
            )
        else:
            clusters_taxopath = None
        wrappers.run_gappa_assign(
            jplace=jplace,
            taxonomy_file=temptax,
            output_dir=output_dir,
            output_prefix=output_prefix,
            only_best_hit=only_best_hit,
            additional_args=gappa_additional_args,
        )

    with TemporaryFilePath() as tempout, TemporaryFilePath() as tempout2, TemporaryFilePath() as tempout3:
        parse_gappa_assign_table(
            input_table=gappa_assign_out,
            has_cluster_id=has_cluster_id,
            cluster_scores=ref_cluster_scores,
            output_file=tempout,
            clusters_taxopath=clusters_taxopath,
        )

        if only_unique_cluster:
            filter_non_unique_placement_assignments(
                placed_tax_assignments=tempout, output_file=tempout2
            )
            shutil.move(tempout2, tempout)

        pick_taxopath_with_highest_lwr(placed_tax_assignments=tempout, output_file=tempout3)

        if query_labels is not None:
            add_query_labels_to_assign_table(
                input_table=tempout3,
                query_labels=query_labels,
                output_table=query_taxo_out,
            )
        else:
            shutil.move(tempout3, query_taxo_out)


def parse_duplicates_from_seqkit(query_duplicates: Path) -> None:
    """
    Add a column with the query IDs of duplicated sequences to the taxonomy assignments file.
    @params:
    taxtable: tsv with query taxonomic and cluster assignments as output by labelplacements.py
    query_duplicates: text file containing query duplicate IDs as output by seqkit rmdup
    """
    df = pd.read_csv(query_duplicates, sep="\t", header=None)
    dup_labels = {
        dup_pair[1]
        .split(",")[0]
        .strip(): ",".join(
            [seq_name.strip() for seq_name in dup_pair[1].split(",")[1:]]
        )
        .strip()
        for dup_pair in df.values
    }
    return dup_labels


def add_duplicates_to_assignment_table(
    taxtable: Path, query_duplicates: Path, output_file: Path = None
) -> None:
    """
    Add duplicated query IDs to cluster and taxonomic assignment table
    """
    if output_file is None:
        output_file = set_default_output_path(taxtable, tag="_duplicates")
    else:
        output_file = Path(output_file)

    assigns = pd.read_csv(taxtable, sep="\t", index_col=False)
    dup_dict = parse_duplicates_from_seqkit(query_duplicates)
    for query_name, duplicate_string in dup_dict.items():
        duplicates = duplicate_string.split(",")

        for duplicate in duplicates:
            row = assigns[assigns.query_name == query_name].copy()
            if not row.empty:
                row = row.iloc[0, :]
                row.query_name = duplicate.strip()
                assigns = assigns.append(row, ignore_index=True)
    assigns.to_csv(output_file, sep="\t", index=False)


def find_queries_placed_in_several_clusters(
    placed_tax_assignments: Path,
) -> tuple[list, pd.DataFrame]:
    """
    Find queries placed in more than one clusterrunGappaAssign
    """
    df = pd.read_csv(placed_tax_assignments, sep="\t", dtype=str)
    dfu = df.groupby("query_id")["cluster_id"].agg(["unique"])

    queries_in_more_than_one_cluster = []
    for i, row in dfu.iterrows():
        if len(row.item()) > 1:
            queries_in_more_than_one_cluster.append(row.name)
    return queries_in_more_than_one_cluster, dfu


def filter_non_unique_placement_assignments(
    placed_tax_assignments: Path, output_file: Path = None
) -> None:
    """
    Remove queries that were assigned to more than one cluster from placements assignments table
    """
    if output_file is None:
        output_file = set_default_output_path(placed_tax_assignments, tag="_filtered")
    else:
        output_file = Path(output_file)

    df = pd.read_csv(placed_tax_assignments, sep="\t")
    queries_in_more_than_one_cluster, _ = find_queries_placed_in_several_clusters(
        placed_tax_assignments
    )
    fdf = df[~df.query_id.isin(queries_in_more_than_one_cluster)].set_index("query_id")
    fdf.to_csv(output_file, sep="\t")


def pick_taxopath_with_highest_lwr(
    placed_tax_assignments: Path, output_file: Path = None
) -> None:
    """
    Pick taxopath assigment with higuest LWR for each placed query
    """
    if output_file is None:
        output_file = set_default_output_path(
            placed_tax_assignments, tag="_unique_taxopath"
        )
    else:
        output_file = Path(output_file)

    df = pd.read_csv(placed_tax_assignments, sep="\t")
    df.groupby("query_id", group_keys=["LWR"]).aggregate("max").to_csv(
        output_file, sep="\t"
    )
