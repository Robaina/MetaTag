#!/bin/bash

# ***************************************************************** #
#                           nifH pipeline                           #
# ***************************************************************** #

# # Preprocess
# python3 ./code/preprocess.py \
#  --in data/MAR_HQ \
#  --outfile data/marhq_cleaned.faa

# # Make database
# python3 ./code/makedatabase.py \
#  --in data/marhq_cleaned.faa \
#  --outdir genes/nifH/results/ \
#  --hmm genes/nifH/hmms/TIGR01287.1.HMM \
#  --max_size 800 \
#  --relabel

# # Alignment and tree
# python3 ./code/buildtree.py \
#  --in genes/nifH/results/ref_database.faa \
#  --outdir genes/nifH/results/ \
#  --msa_method "muscle" \
#  --tree_model "LG+I+G4" \
#  --tree_method "iqtree"

# # Remove tree branch outliers
# python3 ./code/removetreeoutliers.py \
#  --tree genes/nifH/results/ref_database.contree \
#  --outdir genes/nifH/results/ \
#  --aln genes/nifH/results/ref_database.faln

# Classify nifH sequences according to CART model
python3 ./code/classifyNifHsequences.py \
 --seqs genes/nifH/results/ref_database.faa \
 --aln genes/nifH/results/ref_database.faln \
 --indict genes/nifH/results/ref_database_id_dict.pickle \
 --outdict genes/nifH/results/ref_database_id_dict_clustered.pickle \

# Relabel reference tree and msa
python3 ./code/relabeltree.py \
 --tree genes/nifH/results/ref_database_shrink.contree \
 --aln genes/nifH/results/ref_database.faln \
 --labels genes/nifH/results/ref_database_id_dict_clustered.pickle \

 # Commit to GitHub
git add .
git commit -m "Add nifH results, clustered"
git push origin main

# Send notification
# python3 ./code/notify.py --klink https://github.com/Robaina/TRAITS/
