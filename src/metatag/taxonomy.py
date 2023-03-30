#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to assign taxonomy to reference and query (placed) sequences
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import pandas as pd

from metatag.database.labelparsers import LabelParser

matplotlib.use("Agg")  # do not display figs in Jupyter


class Taxopath:
    """
    Object to contain taxopath
    """

    def __init__(self, taxopath_str: str = None, delimiter: str = ";"):
        """_summary_

        Args:
            taxopath_str (str, optional): taxopath as a string. Defaults to None.
            delimiter (str, optional): taxa delimiter in taxopath. Defaults to ";".
        """
        self._taxopath = taxopath_str
        self._delimiter = delimiter
        self._tax_levels = [
            "domain",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ]
        self.taxodict = self._dict_from_taxopath()

    def _dict_from_taxopath(self):
        if self._taxopath is None:
            taxolist = []
        else:
            taxolist = [elem.strip() for elem in self._taxopath.split(self._delimiter)]
        taxolist.extend([None for _ in range(len(self._tax_levels) - len(taxolist))])
        return {taxlevel: taxon for taxlevel, taxon in zip(self._tax_levels, taxolist)}

    @classmethod
    def from_dict(cls, taxodict: dict, delimiter: str = ";") -> Taxopath:
        """Instantiate Taxopath object from dict

        Args:
            taxodict (dict): dict of taxonomic levels and taxa
            delimiter (str, optional): delimiter to separate taxon levels.
                Defaults to ";".

        Returns:
            Taxopath: Taxopath object
        """
        taxa = []
        for taxon in taxodict.values():
            if taxon is None:
                break
            else:
                taxa.append(taxon)
        taxostr = delimiter.join(taxa)
        return cls(taxopath_str=taxostr, delimiter=delimiter)


class TaxonomyAssigner:
    """
    Methods to assign taxonomy to reference sequences
    """

    def __init__(self, taxo_file: Path):
        self._taxodata = (
            pd.read_csv(taxo_file, sep="\t")
            .drop_duplicates(subset="genome")
            .set_index("genome")
        )
        """Methods to assign taxonomy to reference sequences
        """

    @staticmethod
    def lowest_common_taxonomy(taxopaths: list[str]) -> str:
        """Find lowest common taxonomy among set of taxopaths

        Args:
            taxopaths (list[str]): list of taxopath strings

        Returns:
            str: lowest common taxopaht
        """
        data = pd.DataFrame([t.split(";") for t in taxopaths])
        taxlevels = ["d__", "p__", "c__", "o__", "f__", "g__", "s__"]
        lowest_tax = []
        for n, tax_level in enumerate(taxlevels):
            taxa = data.iloc[:, n].drop_duplicates().values
            if len(taxa) == 1:
                lowest_tax.append(taxa[0])
            else:
                break
        return ";".join(lowest_tax)

    def _extract_genome_id_from_label(self, label: str) -> str:
        genome_id = LabelParser.extract_genome_id(label)
        return genome_id

    def assign_taxonomy_to_label(self, label: str) -> str:
        """Assign GTDB taxonomy to label based on genome ID

        Args:
            label (str): reference label containing genome ID

        Returns:
            str: GTDB taxonomy as a string taxopath
        """
        genome_id = self._extract_genome_id_from_label(label)
        if genome_id in self._taxodata.index:
            return self._taxodata.loc[genome_id].item()
        else:
            return "No_taxonomy_found"

    def assign_lowest_common_taxonomy_to_labels(self, labels: list[str]) -> str:
        """Assing taxonomy to set of labels and find lowest common taxonomy
        among them

        Args:
            labels (list[str]): list of reference labels containing genome IDs

        Returns:
            str: lowest common GTDB taxonomy
        """
        taxopaths = [
            taxopath
            for taxopath in map(self.assign_taxonomy_to_label, labels)
            if taxopath != "No_taxonomy_found"
        ]
        try:
            if taxopaths:
                return self.lowest_common_taxonomy(taxopaths)
            else:
                return "Unspecified"
        except Exception:
            return "Unspecified"

    def assign_lowest_common_taxonomy_to_clusters(
        self, clusters: dict, label_dict: dict = None
    ) -> dict:
        """Find lowest possible common taxonomy to reference labels in clusters
        If reference labels do not contain genome IDs, a dictionary, label_dict,
        of reference labels and genome ids (or labels with genome ids) must be passed

        Args:
            clusters (dict): dictionary with keys as cluster IDs and values as lists
                of reference labels in each cluster
            label_dict (dict, optional): dictionary with keys as short IDs and values
                as reference (full) labels. Defaults to None.

        Returns:
            dict: dictionary with keys as cluster IDs and values as lowest common
                taxopath for each cluster
        """
        clusters_taxopath = {}
        for cluster_id, cluster in clusters.items():
            if label_dict is not None:
                cluster_labels = [label_dict[ref_id] for ref_id in cluster]
            else:
                cluster_labels = cluster
            taxopath = self.assign_lowest_common_taxonomy_to_labels(cluster_labels)
            clusters_taxopath[cluster_id] = taxopath
        return clusters_taxopath

    def build_gappa_taxonomy_table(
        self, ref_id_dict: dict, output_file: Path = None
    ) -> None:
        """Build gappa-compatible taxonomy file as specified here:
        https://github.com/lczech/gappa/wiki/Subcommand:-assign
        Removes references without assigned taxonomy

        Args:
            ref_id_dict (dict): dictionary with keys as reference IDs and values as
                reference labels
            output_file (Path, optional): path to output file. Defaults to None.
        """
        if output_file is None:
            output_file = Path().resolve() / "gappa_taxonomy.tsv"
        else:
            output_file = Path(output_file).resolve()

        with open(output_file, "w") as outfile:
            lines = []
            for ref_id, label in ref_id_dict.items():
                taxon_str = self.assign_taxonomy_to_label(label)
                taxon_str = (
                    "Unspecified" if ("No_taxonomy_found" in taxon_str) else taxon_str
                )
                # if taxon_str != 'Unspecified':
                lines.append(f"{ref_id}\t{taxon_str}\n")
            outfile.writelines(lines)


class TaxonomyCounter:
    def __init__(self, taxopath_list: list[str]):
        """
        Tools to summarize taxonomical diversity in a list of taxopaths
        """
        taxodicts = [Taxopath(taxopath_str).taxodict for taxopath_str in taxopath_list]
        self._taxohits = pd.DataFrame(taxodicts).applymap(
            lambda x: "Unspecified" if x is None else x
        )

    def get_counts(
        self,
        taxlevel: str = "family",
        output_tsv: Path = None,
        plot_type: str = "bar",
        output_pdf: Path = None,
    ) -> None:
        """Compute counts and fraction at specified taxonomy levels

        Args:
            taxlevel (str, optional): tanoxomy level to perform counts at.
                Defaults to "family".
            output_tsv (Path, optional): path to output file. Defaults to None.
            plot_type (str, optional): choose either "bar" or "pie". Defaults to "bar".
            output_pdf (Path, optional): path to output pdf with figures.
                Defaults to None.

        Returns:
            pd.DataFrame: dataframe with counts and fraction at specified taxlevel
        """
        counts = self._taxohits[taxlevel].value_counts(normalize=False)
        # Merge counts from "Unspecified" and empty tax level, e.g., "f__"
        matched_row = counts.index[counts.index.str.fullmatch("[a-zA-Z]__")]
        if not matched_row.empty:
            empty_tax_level = matched_row.item()
            counts[counts.index == "Unspecified"] = (
                counts[counts.index == "Unspecified"] + counts[empty_tax_level].item()
            )
            counts = counts.drop(labels=empty_tax_level)

        counts.index.name = taxlevel
        df = counts.to_frame(name="counts")
        df["fraction"] = df.counts.apply(lambda x: x / df.counts.sum())
        if output_tsv is not None:
            df.to_csv(output_tsv, sep="\t")
        if output_pdf is not None:
            self.plot_counts(
                df,
                plot_type=plot_type,
                output_pdf=output_pdf,
                title=f"Taxonomy at level: {taxlevel}",
            )

    def plot_counts(
        self,
        count_data: pd.DataFrame,
        output_pdf: Path,
        plot_type: str = "bar",
        title: str = None,
    ) -> None:
        """Make (and optionally export) barplot ('bar') or pieplot ('pie')
        figure depicting counting results at specified taxonomic level

        Args:
            count_data (pd.DataFrame): dataframe with counts and fraction
                at specified taxonomy level as returned by get_counts()
            plot_type (str, optional): choose between "bar" and "pie".
                Defaults to "bar".
            output_pdf (Path, optional): path to output pdf containing figure.
                Defaults to None.
            title (str, optional): figure title. Defaults to None.
        """
        if title is None:
            title = ""
        if plot_type == "pie":
            count_data.counts.plot.pie(
                figsize=(15, 15), title=title, rotatelabels=True
            ).get_figure().savefig(output_pdf, format="pdf")
        else:
            count_data.counts.plot.bar(
                figsize=(15, 15), title=title, rot=90
            ).get_figure().savefig(output_pdf, format="pdf")
