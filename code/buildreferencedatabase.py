#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Build common (base) database by merging:

1. MAR complete and partial genomes
2. Genomes in Oceans Microbiomics database (except MAR)
3. Genomes in OceanDNA

In all cases, keep genomes with available GTDB taxonomy and
remove duplicates
"""

import os
import shutil
import pyfastx
import pandas as pd

from phyloplacement.utils import terminalExecute, extractTarFile, fullPathListDir
from phyloplacement.database.parsers.mardb import MARdbLabelParser

data_dir = "/home/robaina/Documents/TRAITS/data/"
work_dir = "/home/robaina/Documents/TRAITS/code/"

# Raw databases
mar_partial = os.path.join(data_dir, "databases/aMARPartialHighQuality.tar.gz")
mar_complete = os.path.join(data_dir, "databases/aMARgenomesCompleteWithMMP.tar.gz")

oceandna_rep = os.path.join(data_dir, "databases/afasta_species-representatives.tar.gz")
oceandna_nonrep = os.path.join(data_dir, "databases/afasta_non-representatives.tar.gz")

paoli = '/home/robaina/Documents/TRAITS/data/databases/aPaoli.tar.gz'


# ******************************************************************
#                             Merge MAR
# ******************************************************************

# Extract files
extractTarFile(mar_partial, dest_dir=os.path.join(data_dir, "MAR"))
extractTarFile(mar_complete, dest_dir=os.path.join(data_dir, "MAR"))

all_mar_dir = os.path.join(data_dir, "MAR/allMAR/")
shutil.copytree(
    os.path.join(data_dir, "MAR/aMARgenomesComplete"),
    all_mar_dir,
    dirs_exist_ok=True
    )
shutil.copytree(
    os.path.join(data_dir, "MAR/aMARPartialHighQuality"),
    all_mar_dir,
    dirs_exist_ok=True
    )

all_mar_file = os.path.join(data_dir, "MAR/allMAR.faa")
cmd_str = (
    f"python preprocess.py --in {all_mar_dir} --outfile {all_mar_file}"
)
terminalExecute(cmd_str, work_dir=work_dir)

# Remove records without GTDB taxonomy
marparser = MARdbLabelParser()
mar_taxonomy = os.path.join(data_dir, "taxonomy/MAR_gtdb_taxonomy.tsv")
filtered_mar_file = os.path.join(data_dir, "preprocessedDatabases/mar.faa")
genomes_with_taxonomy = set(pd.read_csv(mar_taxonomy, sep="\t").genome.values)

with open(filtered_mar_file, "w") as outfile:
    for name, seq in pyfastx.Fasta(all_mar_file, build_index=False, full_name=True):
        genome_id = marparser.extractMMPid(name)
        if genome_id in genomes_with_taxonomy:
            outfile.write(f">{name}\n")
            outfile.write(f"{seq}\n")

shutil.rmtree(all_mar_dir)
shutil.remove(all_mar_file)


# ********************************************************************
#               Merge Paoli (Oceans Microbiomics Data)
# ********************************************************************

# Extract files
extractTarFile(paoli, dest_dir=os.path.join(data_dir, "Paoli/"))

# Remove files with MAR genomes
for file in fullPathListDir(os.path.join(data_dir, "Paoli/aPaoli/")):
    if 'MARD_' in file:
        os.remove(file)

all_paoli_dir = os.path.join(data_dir, "Paoli/aPaoli/")
all_paoli_file = os.path.join(data_dir, "Paoli/allPaoli.faa")
cmd_str = (
    f"python preprocess.py --in {all_paoli_dir} --outfile {all_paoli_file}"
)
terminalExecute(cmd_str, work_dir=work_dir)

# Remove records without GTDB taxonomy
paoli_taxonomy = os.path.join(data_dir, "taxonomy/Paoli_gtdb_taxonomy.tsv")
filtered_paoli_file = os.path.join(data_dir, "preprocessedDatabases/paoli.faa")
genomes_with_taxonomy = set(pd.read_csv(paoli_taxonomy, sep="\t").genome.values)

with open(filtered_paoli_file, "w") as outfile:
    for name, seq in pyfastx.Fasta(all_paoli_file, build_index=False, full_name=True):
        genome_id = name.split("__")[0]
        if genome_id in genomes_with_taxonomy:
            outfile.write(f">{name}\n")
            outfile.write(f"{seq}\n")

shutil.remove(all_paoli_file)


# ********************************************************************
#                           Merge OceanDNA
# ********************************************************************

# Extract files
extractTarFile(oceandna_rep, dest_dir=os.path.join(data_dir, "OceanDNA"))
extractTarFile(oceandna_nonrep, dest_dir=os.path.join(data_dir, "OceanDNA"))

all_oceandna_dir = os.path.join(data_dir, "OceanDNA/allOceanDNA/")
shutil.copytree(
    os.path.join(data_dir, "OceanDNA/afasta_non-representatives"),
    all_oceandna_dir,
    dirs_exist_ok=True
    )
shutil.copytree(
    os.path.join(data_dir, "OceanDNA/afasta_species-representatives"),
    all_oceandna_dir,
    dirs_exist_ok=True
    )

all_oceandna_file = os.path.join(data_dir, "OceanDNA/allOceanDNA.faa")
cmd_str = (
    f"python preprocess.py --in {all_oceandna_dir} --outfile {all_oceandna_file}"
)
terminalExecute(cmd_str, work_dir=work_dir)

# Remove records without GTDB taxonomy
oceandna_taxonomy = os.path.join(data_dir, "taxonomy/OceanDNA_gtdb_taxonomy.tsv")
filtered_oceandna_file = os.path.join(data_dir, "preprocessedDatabases/oceandna.faa")
genomes_with_taxonomy = set(pd.read_csv(oceandna_taxonomy, sep="\t").genome.values)

with open(filtered_oceandna_file, "w") as outfile:
    for name, seq in pyfastx.Fasta(all_oceandna_file, build_index=False, full_name=True):
        genome_id = name.split("__")[0]
        if genome_id in genomes_with_taxonomy:
            outfile.write(f">{name}\n")
            outfile.write(f"{seq}\n")

shutil.rmtree(all_oceandna_dir)
shutil.remove(all_oceandna_file)


# ********************************************************************
#                          Merge all filtered databases
# ********************************************************************

preprocessed_dir = os.path.join(data_dir, "preprocessedDatabases")
final_database_file = os.path.join(data_dir, "final_ref_database.faa")
cmd_str = (
    f"python preprocess.py --in {preprocessed_dir} --outfile {final_database_file}"
)
terminalExecute(cmd_str, work_dir=work_dir)