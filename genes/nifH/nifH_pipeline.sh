#!/bin/bash

# ***************************************************************** #
#                           nifH pipeline                           #
# ***************************************************************** #

# Notes: TIGR01287 only returns 221 hits with full database

rm -r /home/robaina/Documents/TRAITS/genes/nifH/results/; mkdir /home/robaina/Documents/TRAITS/genes/nifH/results/

# # Make database
python3 ./code/makedatabase.py \
 --in data/final_ref_database.fasta \
 --outdir genes/nifH/results/ \
 --hmm genes/nifH/hmms/TIGR01287.1.HMM \
 --max_size 800 \
 --relabel

# Alignment and tree
python3 ./code/buildtree.py \
 --in genes/nifH/results/ref_database.faa \
 --outdir genes/nifH/results/ \
 --msa_method "muscle" \
 --tree_model "LG+I+G4" \
 --tree_method "iqtree"

# Remove tree branch outliers
python3 ./code/removetreeoutliers.py \
 --tree genes/nifH/results/ref_database.newick \
 --outdir genes/nifH/results/ \
 --aln genes/nifH/results/ref_database.faln

# Classify nifH sequences according to CART model
python3 ./code/classifyNifHsequences.py \
 --seqs genes/nifH/results/ref_database.faa \
 --aln genes/nifH/results/ref_database.faln \
 --indict genes/nifH/results/ref_database_id_dict.pickle \
 --outdict genes/nifH/results/ref_database_id_dict_clustered.pickle \
 --clusters_file genes/nifH/data/clusters.tsv

# Relabel reference tree and msa
python3 ./code/relabeltree.py \
 --tree genes/nifH/results/ref_database_shrink.newick \
 --aln genes/nifH/results/ref_database.faln \
 --labels genes/nifH/results/ref_database_id_dict_clustered.pickle \

# # Commit to GitHub
# git add .
# git commit -m "Add nifH results, clustered"
# git push origin main

# Send notification
python3 ./code/notify.py --link https://github.com/Robaina/TRAITS/
