#!/bin/bash

# ***************************************************************** #
#                           nifH pipeline                            #
# ***************************************************************** #

# mkdir -p mardb_mc16

# Preprocess
python3 ./code/preprocess.py \
 --in data/MAR_HQ \
 --outfile data/marhq_cleaned.faa

# Make database
python3 ./code/makedatabase.py \
 --in data/marhq_cleaned.faa \
 --outdir genes/nifH/results/ \
 --hmm genes/nifH/hmms/TIGR01287.1.HMM \
 --max_size 800

# Alignment and tree
python3 ./code/buildtree.py \
 --in genes/nifH/results/ref_database.faa \
 --outdir genes/nifH/results/ \
 --msa_method "muscle" \
 --tree_model "TEST" \
 --tree_method "iqtree"

# Remove tree branch outliers
python3 ./code/removetreeoutliers.py \
 --tree genes/nifH/results/ref_database.contree \
 --outdir genes/nifH/results/ \
 --aln genes/nifH/results/ref_database.faln

# Relabel reference tree and msa
python3 ./code/relabeltree.py \
 --tree genes/nifH/results/ref_database_shrink.contree \
 --aln genes/nifH/results/ref_database.faln \
 --labels genes/nifH/results/ref_database_id_dict.pickle \

# # Classify nifH sequences according to CART model
# python3 ./code/classifyNifHsequences.py \
#  --seqs genes/nifH/results/ref_database.faa \
#  --aln genes/nifH/results/ref_database.faln \
#  --outfile genes/nifH/results/ref_database_classified.faa \

 # Commit to GitHub
# git add -f /home/robaina/Documents/TRAITS/genes/nifH/nxr_results/narG_molyb_nitz18_nitz21/ref_database_shrink_relabel.contree
# git commit -m "Update data"
# git push origin robaina_nxr

# Send notification
python3 ./code/notify.py --link https://github.com/Robaina/TRAITS/
