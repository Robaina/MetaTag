#!/bin/bash

# ***************************************************************** #
#                           ureC pipeline                           #
# ***************************************************************** #

Preprocess
python3 ./code/preprocess.py \
 --in data/MAR_HQ \
 --outfile data/marhq_cleaned.faa

# Make database
python3 ./code/makedatabase.py \
 --in data/marhq_cleaned.faa \
 --outdir genes/ureC/results/ \
 --hmm genes/ureC/hmms/alpha_TIGR01792.1.HMM \
 --prefix "ureC_" --relabel \
 --max_size 800

# Add Koper 2004 ureC sequences
python3 ./code/preprocess.py \
 --in genes/ureC/data/koper_2004_seqs.fasta \
 --outfile genes/ureC/results/koper2004_seqs_short_ids.faa \
 --idprefix "ref_ko04_" --relabel

# Move databases to directory to merge
mv genes/ureC/results/ureC_ref_database.faa genes/ureC/results/mardb_ko04/
mv genes/ureC/results/koper2004_seqs_short_ids.faa genes/ureC/results/mardb_ko04/

# Merge all four databases into final reference database
python3 ./code/preprocess.py \
 --in genes/ureC/results/mardb_ko04/ \
 --outfile genes/ureC/results/ref_database.faa

# Alignment and tree
python3 ./code/buildtree.py \
 --in genes/ureC/results/ref_database.faa \
 --outdir genes/ureC/results/ \
 --msa_method "muscle" \
 --tree_model "TEST" \
 --tree_method "iqtree"

# Remove tree branch outliers
python3 ./code/removetreeoutliers.py \
 --tree genes/ureC/results/ref_database.contree \
 --outdir genes/ureC/results/ \
 --aln genes/ureC/results/ref_database.faln

# Relabel reference tree and msa
python3 ./code/relabeltree.py \
 --tree genes/ureC/results/ref_database_shrink.contree \
 --aln genes/ureC/results/ref_database.faln \
 --labels genes/ureC/results/ureC_ref_database_id_dict.pickle \
          genes/ureC/results/koper2004_seqs_short_ids_id_dict.pickle \
 --label_prefixes "ureC_" "ko04_"

# Commit to GitHub
# git add .
# git commit -m "Add ureC results"
# git push origin main

# Send notification
# python3 ./code/notify.py --link https://github.com/Robaina/TRAITS/
