#!/bin/bash

# ***************************************************************** #
#                           ureC pipeline                           #
# ***************************************************************** #
# rm -r /home/robaina/Documents/TRAITS/genes/ureC/results
# mkdir /home/robaina/Documents/TRAITS/genes/ureC/results

# # Make database
# python3 ./code/makedatabase.py \
#  --in data/final_ref_database.fasta \
#  --outdir genes/ureC/results/ \
#  --hmms genes/ureC/hmms/alpha_TIGR01792.1.HMM \
#  --prefix "ureC_" --relabel \
#  --relabel_prefixes "ureC_" \
#  --max_sizes 800

# # Add Koper 2004 ureC sequences
# python3 ./code/preprocess.py \
#  --in genes/ureC/data/koper_2004_seqs.fasta \
#  --outfile genes/ureC/results/koper2004_seqs_short_ids.faa \
#  --idprefix "ko04_" --relabel

# # Add Holn 1997 amidohydrolases sequences (E.coli and M.jannaschii)
# python3 ./code/preprocess.py \
#  --in genes/ureC/data/holn_1997_amidohydrolases.fasta \
#  --outfile genes/ureC/results/holn1997_seqs_short_ids.faa \
#  --idprefix "ho97_" --relabel

# # Move databases to directory to merge
# mkdir -p genes/ureC/results/mardb_ko04_ho97/
# mv genes/ureC/results/ureC_ref_database.faa genes/ureC/results/mardb_ko04_ho97/
# mv genes/ureC/results/koper2004_seqs_short_ids.faa genes/ureC/results/mardb_ko04_ho97/
# mv genes/ureC/results/holn1997_seqs_short_ids.faa genes/ureC/results/mardb_ko04_ho97/

# # Merge all four databases into final reference database
# python3 ./code/preprocess.py \
#  --in genes/ureC/results/mardb_ko04_ho97/ \
#  --outfile genes/ureC/results/ref_database.faa

# Alignment and tree
# python3 ./code/buildtree.py \
#  --in genes/ureC/results/ref_database_no_outliers.faa \
#  --outdir genes/ureC/results/ \
#  --msa_method "muscle" \
#  --tree_model "LG+F+I+G4" \
#  --tree_method "iqtree"

# # Remove tree branch outliers
# python3 ./code/removetreeoutliers.py \
#  --tree genes/ureC/results/ref_database.newick \
#  --outdir genes/ureC/results/ \
#  --aln genes/ureC/results/ref_database.faln \
#  --additional_args "-q 0.15"

# Relabel reference tree and msa
python3 ./code/relabeltree.py \
 --tree genes/ureC/results/ref_database.newick \
 --aln genes/ureC/results/ref_database.faln \
 --labels genes/ureC/results/ureC_ref_database_id_dict.pickle \
          genes/ureC/results/koper2004_seqs_short_ids_id_dict.pickle \
          genes/ureC/results/holn1997_seqs_short_ids_id_dict.pickle \
 --label_prefixes "ureC_" "ko04_" "ho97_"

# # Commit to GitHub
# git add .
# git commit -m "Add ureC results"
# git push origin main

# Send notification
# python3 ./code/notify.py --link https://github.com/Robaina/TRAITS/
