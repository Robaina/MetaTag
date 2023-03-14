#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import sys
from importlib import metadata

from metatag.scripts import (
    buildtree,
    countplacements,
    labelplacements,
    makedatabase,
    placesequences,
    plottree,
    preprocess,
    relabeltree,
)

meta = metadata.metadata("metatag")
__version__ = meta["Version"]
__author__ = meta["Author"]


class MetaTag:
    """Main command class"""

    def __init__(
        self, subcommand: str, subcommand_args: list[str], silent: bool = False
    ):
        """Initialize main command

        Args:
            subcommand (str): subcommand name.
            subcommand_args (list[str]): list of subcommand arguments and values.
            silent (bool): do not print Pynteny header to command line.
        """
        self._subcommand = subcommand
        if "--hmmsearch_args" in subcommand_args:
            hmm_arg_idx = subcommand_args.index("--hmmsearch_args") + 1
            subcommand_args[hmm_arg_idx] = " " + subcommand_args[hmm_arg_idx]
        self._subcommand_args = subcommand_args
        parser = argparse.ArgumentParser(
            description=(self._printLogo(silent)),
            usage=("metatag <subcommand> [-h] [args] \n"),
            epilog=(""),
            formatter_class=argparse.RawTextHelpFormatter,
        )
        parser._positionals.title = "subcommands"
        parser.add_argument(
            help=(
                "preprocess \n"
                "database \n"
                "tree \n"
                "place \n"
                "assign \n"
                "count \n"
                "plot \n"
                "relabel \n"
                "cite \n"
            ),
            dest="subcommand",
            metavar="",
        )
        parser.add_argument(
            "-v",
            "--version",
            help="show version and exit",
            action="version",
            version=__version__,
        )
        if len(sys.argv) < 2:
            parser.print_help()
            sys.exit(1)
        args = parser.parse_args(self._subcommand)
        input_subcommand = getattr(args, "subcommand")
        self.call_subcommand(subcommand_name=input_subcommand)

    def _printLogo(self, silent: str = False) -> None:
        if not silent:
            print(
                (
                    """
                 _     _______
                | |   |__   __|
  _ __ ___   ___| |_ __ _| | __ _  __ _
 | '_ ` _ \ / _ \ __/ _` | |/ _` |/ _` |
 | | | | | |  __/ || (_| | | (_| | (_| |
 |_| |_| |_|\___|\__\__,_|_|\__,_|\__, |
                                   __/ |
                                  |___/
"""
                    f"Metagenome function and taxonomy assignment through evolutionary placement, v{__version__}\n"
                    "Semidán Robaina Estévez (srobaina@ull.edu.es), 2022\n"
                    " \n"
                )
            )

    def call_subcommand(self, subcommand_name: str) -> None:
        subcommand = getattr(self, subcommand_name)
        subcommand()

    def preprocess(self):
        """Call preprocess subcommand."""
        args = preprocess.initialize_parser().parse_args(self._subcommand_args)
        preprocess.run(args)

    def database(self):
        """Call database subcommand."""
        args = makedatabase.initialize_parser().parse_args(self._subcommand_args)
        makedatabase.run(args)

    def tree(self):
        """Call tree subcommand."""
        args = buildtree.initialize_parser().parse_args(self._subcommand_args)
        buildtree.run(args)

    def place(self):
        """Call place subcommand."""
        args = placesequences.initialize_parser().parse_args(self._subcommand_args)
        placesequences.run(args)

    def assign(self):
        """Call place subcommand"""
        args = labelplacements.initialize_parser().parse_args(self._subcommand_args)
        labelplacements.run(args)

    def count(self):
        """Call count subcommand"""
        args = countplacements.initialize_parser().parse_args(self._subcommand_args)
        countplacements.run(args)

    def plot(self):
        """Call plot subcommand"""
        args = plottree.initialize_parser().parse_args(self._subcommand_args)
        plottree.run(args)

    def relabel(self):
        """Call relabel subcommand"""
        args = relabeltree.initialize_parser().parse_args(self._subcommand_args)
        relabeltree.run(args)

    def cite(self):
        """Print pynteny's citation string"""
        citation = (
            "Semidán Robaina Estévez (2022). MetaTag: Metagenome functional and taxonomical annotation through phylogenetic tree placement."
            f"(Version {__version__}). Zenodo."
        )
        print("If you use this software, please cite it as below: ")
        print(citation)
        return citation


def main():
    subcommand, subcommand_args = sys.argv[1:2], sys.argv[2:]
    MetaTag(subcommand, subcommand_args)


if __name__ == "__main__":
    main()
