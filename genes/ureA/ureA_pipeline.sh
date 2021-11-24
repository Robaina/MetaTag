#!/bin/bash

# ***************************************************************** #
#                           ureA pipeline                           #
# ***************************************************************** #

# Preprocess
python3 ./code/preprocess.py \
 --in data/MAR_HQ \
 --outfile data/marhq_cleaned.faa

# Make database
python3 ./code/makedatabase.py \
 --in data/marhq_cleaned.faa \
 --outdir genes/ureA/results/ \
 --hmm genes/ureA/hmms/alpha_TIGR01792.1.HMM \
 --max_size 800 \
 --relabel

# Alignment and tree
python3 ./code/buildtree.py \
 --in genes/ureA/results/ref_database.faa \
 --outdir genes/ureA/results/ \
 --msa_method "muscle" \
 --tree_model "TEST" \
 --tree_method "iqtree"

# Remove tree branch outliers
python3 ./code/removetreeoutliers.py \
 --tree genes/ureA/results/ref_database.contree \
 --outdir genes/ureA/results/ \
 --aln genes/ureA/results/ref_database.faln

# Relabel reference tree and msa
python3 ./code/relabeltree.py \
 --tree genes/ureA/results/ref_database_shrink.contree \
 --aln genes/ureA/results/ref_database.faln \
 --labels genes/ureA/results/ref_database_id_dict_clustered.pickle \

 # Commit to GitHub
# git add .
# git commit -m "Add ureA results, clustered"
# git push origin main

# Send notification
# python3 ./code/notify.py --klink https://github.com/Robaina/TRAITS/
