#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import sys
import tempfile
from importlib import metadata
from pathlib import Path

from metatag.scripts import (
    preprocess,
    makedatabase,
    buildtree,
    placesequences,
    labelplacements,
    countplacements,
    relabeltree,
)

meta = metadata.metadata("metatag")
__version__ = meta["Version"]
__author__ = meta["Author"]


class MetaTag:
    """Main command class"""

    def __init__(self, subcommand: str, subcommand_args: list[str]):
        """Initialize main command

        Args:
            subcommand (str): subcommand name
            subcommand_args (list[str]): list of subcommand arguments and values.
        """
        self._subcommand = subcommand
        if "--hmmsearch_args" in subcommand_args:
            hmm_arg_idx = subcommand_args.index("--hmmsearch_args") + 1
            subcommand_args[hmm_arg_idx] = " " + subcommand_args[hmm_arg_idx]
        self._subcommand_args = subcommand_args
        parser = argparse.ArgumentParser(
            description=(self._printLogo()),
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
        self._call_subcommand(subcommand_name=input_subcommand)

    def _printLogo(self):
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

    def _call_subcommand(self, subcommand_name: str) -> None:
        subcommand = getattr(self, subcommand_name)
        subcommand()

    def preprocess(self):
        """Call preprocess subcommand."""
        preprocess.main()

    def database(self):
        """Call database subcommand."""
        makedatabase.main()

    def tree(self):
        """Call tree subcommand."""
        buildtree.main()

    def place(self):
        """Call place subcommand."""
        placesequences.main()

    def assign(self):
        """Call place subcommand"""
        labelplacements.main()

    def count(self):
        """Call count subcommand"""
        countplacements.main()

    def relabel(self):
        """Call relabel subcommand"""
        relabeltree.main()

    def cite(self):
        """Print pynteny's citation string"""
        citation = (
            "Semidán Robaina Estévez (2022). MetaTag: Metagenome functional and taxonomical annotation through phylogenetic tree placement."
            f"(Version {__version__}). Zenodo. https://doi.org/10.5281/zenodo.7048685"
        )
        print("If you use this software, please cite it as below: ")
        print(citation)
        return citation


class SubcommandParser:
    """Argparse parsers for Pynteny's subcommands"""

    @staticmethod
    def get_help_str(subcommand: str) -> str:
        """Get help string for subcommand.

        Args:
            subcommand (str): subcommand name.

        Returns:
            str: help string.
        """
        parser = getattr(SubcommandParser, subcommand)()
        with tempfile.NamedTemporaryFile(mode="w+") as file:
            parser.print_help(file)
            file.flush()
            with open(file.name, encoding="UTF-8") as help_file:
                help_str = help_file.read()
        return help_str


def main():
    subcommand, subcommand_args = sys.argv[1:2], sys.argv[2:]
    print(subcommand)
    MetaTag(subcommand, subcommand_args)


if __name__ == "__main__":
    main()
