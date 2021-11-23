#!/bin/bash

# ***************************************************************** #
#                           Amt pipeline                            #
# ***************************************************************** #

mkdir -p mardb_mc16

# Preprocess
python3 ./code/preprocess.py \
 --in data/MAR_HQ \
 --outfile data/marhq_cleaned.faa

# Make database
python3 ./code/makedatabase.py \
 --in data/marhq_cleaned.faa \
 --outdir genes/amt/results/ \
 --hmm genes/amt/hmms/TIGR00836.1.HMM \
 --prefix "amt_" --relabel \
 --max_size 600

# Add amt-classified bacterial and archaeal sequences from mcdonald2016 
python3 ./code/preprocess.py \
 --in genes/amt/data/mcdonald2016_prokaryotes_parsed.fasta \
 --outfile genes/amt/results/mcdonald2016_prokaryotes_short_ids.faa \
 --idprefix "ref_mc16_" --relabel

# Add RH-classified sequences from mcdonald2016 
python3 ./code/preprocess.py \
 --in genes/amt/data/mcdonald2016_rh_parsed.fasta \
 --outfile genes/amt/results/mcdonald2016_RH_short_ids.faa \
 --idprefix "ref_mc16rh_" --relabel

# Move databases to directory to merge
mv genes/amt/results/amt_ref_database.faa genes/amt/results/mardb_mc16/
mv genes/amt/results/mcdonald2016_prokaryotes_short_ids.faa genes/amt/results/mardb_mc16/
mv genes/amt/results/mcdonald2016_RH_short_ids.faa genes/amt/results/mardb_mc16/

# Merge all four databases into final reference database
python3 ./code/preprocess.py \
 --in genes/amt/results/mardb_mc16/ \
 --outfile genes/amt/results/ref_database.faa

# Alignment and tree
python3 ./code/buildtree.py \
 --in genes/amt/results/ref_database.faa \
 --outdir genes/amt/results/ \
 --msa_method "muscle" \
 --tree_model "LG+F+I+G4" \
 --tree_method "iqtree"

# Remove tree branch outliers
python3 ./code/removetreeoutliers.py \
 --tree genes/amt/results/ref_database.contree \
 --outdir genes/amt/results/ \
 --aln genes/amt/results/ref_database.faln

# Relabel reference tree and msa
python3 ./code/relabeltree.py \
 --tree genes/amt/results/ref_database_shrink.contree \
 --aln genes/amt/results/ref_database.faln \
 --labels genes/amt/results/amt_ref_database_id_dict.pickle \
          genes/amt/results/mcdonald2016_prokaryotes_short_ids_id_dict.pickle \
          genes/amt/results/mcdonald2016_RH_short_ids_id_dict.pickle \
 --label_prefixes "amt_" "mc16_" "mc16rh_"

 # Commit to GitHub
# git add -f /home/robaina/Documents/TRAITS/genes/amt/nxr_results/narG_molyb_nitz18_nitz21/ref_database_shrink_relabel.contree
# git commit -m "Update data"
# git push origin robaina_nxr

# Send notification
# python3 ./code/notify.py --link https://github.com/Robaina/TRAITS/
