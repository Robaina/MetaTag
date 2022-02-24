#!/bin/bash

# ***************************************************************** #
#                           nifH pipeline                           #
# ***************************************************************** #

# Notes: TIGR01287 only returns 221 hits with full database
# Adding azotobacter nifH sequence to database (ref_azo) and pickle

# rm -r /home/robaina/Documents/TRAITS/genes/nifH/results/; mkdir /home/robaina/Documents/TRAITS/genes/nifH/results/

# # Make database
# python3 ./code/makedatabase.py \
#  --in data/final_ref_database.fasta \
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

# # Multi to single line fasta aln
# # awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' genes/nifH/results/ref_database.faln > singleLine.faln
# # mv singleLine.faln genes/nifH/results/ref_database.faln

# # Remove tree branch outliers
# python3 ./code/removetreeoutliers.py \
#  --tree genes/nifH/results/ref_database.newick \
#  --outdir genes/nifH/results/ \
#  --aln genes/nifH/results/ref_database.faln

# # Classify nifH sequences according to CART model
# python3 ./code/classifyNifHsequences.py \
#  --seqs genes/nifH/results/ref_database.faa \
#  --aln genes/nifH/results/ref_database.faln \
#  --indict genes/nifH/results/ref_database_id_dict.pickle \
#  --outdict genes/nifH/results/ref_database_id_dict_clustered.pickle \
#  --clusters_file genes/nifH/data/clusters.tsv

# # Relabel reference tree and msa
# python3 ./code/relabeltree.py \
#  --tree genes/nifH/results/ref_database_shrink.newick \
#  --aln genes/nifH/results/ref_database.faln \
#  --labels genes/nifH/results/ref_database_id_dict_clustered.pickle \

# # Preprocess query sequences
# python3 code/preprocess.py \
#  --in data/environmental/remedios/P14401_101_S5_L002_P_001_sliced_1000.fasta \
#  --outfile genes/nifH/data/query_cleaned.faa \
#  --idprefix "query_" --relabel --translate

# # Place query sequences
# python3 code/placesequences.py \
#  --aln genes/nifH/results/ref_database_shrink.faln \
#  --tree genes/nifH/results/ref_database_shrink.newick \
#  --query genes/nifH/data/query_cleaned.faa \
#  --outdir genes/nifH/results/remedios/ \
#  --aln_method "papara" \
# --tree_model genes/nifH/results/ref_database.log

# # Assign taxonomy to placed sequences  (--outgroup tests/test_results/data/outliers_short_ids.faa \)
# python3 code/labelplacements.py \
#  --jplace genes/nifH/results/remedios/epa_result.jplace \
#  --labels genes/nifH/results/ref_database_id_dict_clustered.pickle \
#  --ref_clusters genes/nifH/data/clusters.tsv \
#  --ref_cluster_scores genes/nifH/data/cluster_scores.tsv \
#  --prefix "placed_tax_" \
#  --outdir genes/nifH/results/remedios/

# Count placements (filter by taxon, cluster id and quality score)
python3 code/countplacements.py \
 --taxtable genes/nifH/results/remedios/placed_tax_assignments.tsv \
 --taxlevel "family" \
 --cluster_ids "cluster_I" "cluster_II" "cluster_III" "cluster_IV" \
 --score_threshold 0.5 \
 --outfile genes/nifH/results/remedios/placed_family_tax_counts.tsv

# # Commit to GitHub
# git add .
# git commit -m "Add nifH results, clustered"
# git push origin main

# Send notification
# python3 ./code/notify.py --link https://github.com/Robaina/TRAITS/
