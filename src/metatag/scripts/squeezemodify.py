#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Replace SqueezeMeta's taxonomy by placement-derived taxonomy
"""
import argparse
from pathlib import Path

from metatag.squeezemeta import (
    SqueezeMetaOutputParser,
    SqueezeMetaTaxonomyParser,
)

parser = argparse.ArgumentParser(
    description=("Replace SqueezeMeta's taxonomy by placement-derived taxonomy"),
    epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2022",
)

optional = parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
parser._action_groups.append(optional)

required.add_argument(
    "--squeeze_table",
    "-s",
    dest="squeeze_table",
    type=Path,
    required=True,
    help="path to SqueezeMeta's output file",
)
required.add_argument(
    "--taxplacement_table",
    "-p",
    dest="taxplacement_table",
    type=Path,
    required=True,
    help="path to placement assignments output file",
)
optional.add_argument(
    "--outfile",
    "-o",
    dest="outfile",
    default=None,
    type=Path,
    help="path to output (modified squeezemeta) file",
)
args = parser.parse_args()


def main():
    with open(args.squeeze_table) as f:
        next(f)
        dataline = next(f)
    if len(dataline.split("\t")) > 2:
        metaparser = SqueezeMetaOutputParser(args.squeeze_table)
    else:
        metaparser = SqueezeMetaTaxonomyParser(args.squeeze_table)
    metaparser.replace_by_placement_taxopath(
        placement_assignments=args.taxplacement_table, output_file=args.outfile
    )


if __name__ == "__main__":
    main()
    print("Done!")
