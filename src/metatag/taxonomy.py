#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to assign taxonomy to reference and query (placed) sequences
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from metatag.database.labelparsers import LabelParser
from metatag.utils import read_from_pickle_file, set_default_output_path


class Taxopath:
    """
    Object to contain taxopath
    """

    def __init__(self, taxopath_str: str = None, delimiter: str = ";"):
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
        """
        Instantiate Taxopath object from dict
        """
        taxa = []
        for taxon in taxodict.values():
            if taxon is None:
                break
            else:
                taxa.append(taxon)
        taxostr = delimiter.join(taxa)
        return cls(taxopath_str=taxostr, delimiter=delimiter)

    # @classmethod
    # def get_lowest_common_taxopath(cls, taxopaths: list[str]) -> Taxopath:
    #     """
    #     compute lowest common taxopath (ancestor) of a list
    #     of taxopaths
    #     """
    #     taxopath_dicts = [cls(taxostr).taxodict for taxostr in taxopaths]
    #     common_taxodict = cls().taxodict
    #     for taxlevel in cls().taxlevels:
    #         taxa = set([taxdict[taxlevel] for taxdict in taxopath_dicts])
    #         if len(taxa) > 1:
    #             break
    #         else:
    #             common_taxodict[taxlevel] = list(taxa)[0]
    #     return cls.from_dict(common_taxodict)

    # @property
    # def taxostring(self):
    #     return self._taxopath

    # @property
    # def taxlevels(self):
    #     return self._tax_levels


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

    @staticmethod
    def lowest_common_taxonomy(taxopaths: list[str]) -> str:
        """
        Find lowest common taxonomy among set of taxopaths
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
        """
        Assign GTDB taxonomy to label based on genome ID
        """
        genome_id = self._extract_genome_id_from_label(label)
        if genome_id in self._taxodata.index:
            return self._taxodata.loc[genome_id].item()
        else:
            return "No_taxonomy_found"

    def assign_lowest_common_taxonomy_to_labels(self, labels: list[str]) -> str:
        """
        Assing taxonomy to set of labels and find lowest common taxonomy
        among them
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
        """
        Find lowest possible common taxonomy to reference labels in clusters
        If reference labels do not contain genome IDs, a dictionary, label_dict,
        of reference labels and genome ids (or labels with genome ids) must be passed
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
        """
        Build gappa-compatible taxonomy file as specified here:
        https://github.com/lczech/gappa/wiki/Subcommand:-assign
        Removes references without assigned taxonomy
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
    ) -> pd.DataFrame:
        """
        Compute counts and fraction at specified taxonomy levels
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
        # Add fraction
        df = counts.to_frame(name="counts")
        df["fraction"] = df.counts.apply(lambda x: x / df.counts.sum())
        if output_tsv is not None:
            df.to_csv(output_tsv, sep="\t")
        # Make figure
        fig = self.plot_counts(
            df,
            plot_type=plot_type,
            output_pdf=output_pdf,
            title=f"Taxonomy at level: {taxlevel}",
        )
        return df, fig

    def plot_counts(
        self,
        count_data: pd.DataFrame,
        plot_type: str = "bar",
        output_pdf: Path = None,
        title: str = None,
    ):
        """
        Make (and optionally export) barplot ('bar') or pieplot ('pie')
        figure depicting counting results at specified taxonomic level
        """
        if title is None:
            title = ""
        if plot_type == "pie":
            fig = count_data.counts.plot.pie(
                figsize=(15, 15), title=title, rotatelabels=True
            ).get_figure()
        else:
            fig = count_data.counts.plot.bar(
                figsize=(15, 15), title=title, rot=90
            ).get_figure()
        if output_pdf is not None:
            fig.savefig(output_pdf, format="pdf")
        return fig


# def evaluate_taxonomy_of_reference_database(
#     label_dict_pickle: Path = None,
#     taxlevels: list[str] = None,
#     output_dir: Path = None,
#     plot_results: bool = False,
# ) -> None:
#     """
#     Assign taxonomy to sequences and evaluate taxonomical representation
#     """
#     if taxlevels is None:
#         taxlevels = ["class", "order", "family", "genus"]
#     if output_dir is None:
#         output_dir = set_default_output_path(label_dict_pickle, only_dir_name=True)
#     else:
#         output_dir = Path(output_dir).resolve()

#     taxonomy = TaxonomyAssigner(taxo_file=Path("./data/merged_taxonomy.tsv"))
#     label_dict = read_from_pickle_file(label_dict_pickle)
#     taxopaths = [
#         taxonomy.assign_taxonomy_to_label(label) for label in label_dict.values()
#     ]
#     taxcounter = TaxonomyCounter(taxopaths)

#     for taxlevel in taxlevels:
#         if plot_results:
#             outpdf = output_dir / f"ref_taxonomy_counts_{taxlevel}.pdf"
#         else:
#             outpdf = None
#         taxcounter.get_counts(
#             taxlevel,
#             output_tsv=output_dir / f"ref_taxonomy_counts_{taxlevel}.tsv",
#             plot_type="bar",
#             output_pdf=outpdf,
#         )
