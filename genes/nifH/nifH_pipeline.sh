#!/bin/bash

# ***************************************************************** #
#                           nifH pipeline                           #
# ***************************************************************** #

# Notes: TIGR01287 only returns 221 hits with full database
# Adding azotobacter nifH sequence to database (ref_azo) and pickle
# Manually finding key amino acids positions in alignment
# Adding outgroup bchX (TIGR02016)

# rm -r /home/robaina/Documents/TRAITS/genes/nifH/results/; mkdir /home/robaina/Documents/TRAITS/genes/nifH/results/

# # Make database
# python3 ./code/makedatabase.py \
#  --in data/final_ref_database.fasta \
#  --outdir genes/nifH/results/ \
#  --hmms genes/nifH/hmms/TIGR01287.1.HMM \
#         genes/nifH/hmms/TIGR02016.1.HMM \
#  --max_sizes 1000 10 \
#  --relabel_prefixes "ref_" "out_" --relabel \
#  --remove_duplicates

# # Alignment and tree
# python3 ./code/buildtree.py \
#  --in genes/nifH/results/ref_database.faa \
#  --outdir genes/nifH/results/ \
#  --msa_method "muscle" \
#  --tree_model "iqtest" \
#  --tree_method "iqtree"

# # Classify nifH sequences according to CART model
# python3 ./code/classifyNifHsequences.py \
#  --seqs genes/nifH/results/ref_database.faa \
#  --aln genes/nifH/results/ref_database.faln \
#  --indict genes/nifH/results/ref_database_id_dict.pickle \
#  --outdict genes/nifH/results/ref_database_id_dict_clustered.pickle \
#  --out_clusters_file genes/nifH/data/clusters.tsv

# # Relabel reference tree and msa
# python3 ./code/relabeltree.py \
#  --tree genes/nifH/results/ref_database.newick \
#  --aln genes/nifH/results/ref_database.faln \
#  --labels genes/nifH/results/ref_database_id_dict_clustered.pickle \

# # Preprocess query sequences
# python3 code/preprocess.py \
#  --in genes/nifH/daniel/DanielNifH.fasta \
#  --outfile genes/nifH/daniel/query_cleaned.faa \
#  --idprefix "query_" --relabel

# # Place query sequences
# python3 code/placesequences.py \
#  --aln genes/nifH/results/ref_database.faln \
#  --tree genes/nifH/results/ref_database.newick \
#  --query genes/nifH/daniel/query_cleaned.faa \
#  --outdir genes/nifH/daniel/ \
#  --aln_method "papara" \
# --tree_model genes/nifH/results/ref_database.log

# Assign taxonomy to placed sequences
python3 code/labelplacements.py \
 --jplace genes/nifH/daniel/epa_result.jplace \
 --labels genes/nifH/results/ref_database_id_dict_clustered.pickle \
          genes/nifH/daniel/query_cleaned_id_dict.pickle \
 --ref_clusters genes/nifH/data/clusters.tsv \
 --ref_cluster_scores genes/nifH/data/cluster_scores.tsv \
 --outgroup "out_" \
 --prefix "placed_tax_" \
 --outdir genes/nifH/daniel/

# # Count placements (filter by taxon, cluster id and quality score)
# python3 code/countplacements.py \
#  --taxtable genes/nifH/daniel/placed_tax_assignments.tsv \
#  --taxlevels "order" "class" "family" "genus" \
#  --cluster_ids "cluster_I" "cluster_II" "cluster_III" \
#  --score_threshold 0.5 \
#  --outdir genes/nifH/daniel/

# # Relabel tree and alignment
# python3 code/relabeltree.py \
#  --tree genes/nifH/daniel/epa_result.newick \
#  --labels genes/nifH/results/ref_database_id_dict_clustered.pickle \
#           genes/nifH/daniel/query_cleaned_id_dict.pickle \
#  --taxonomy \
#  --label_prefixes "ref_" "query_"

# # Commit to GitHub
# git add .
# git commit -m "Add nifH results, clustered"
# git push origin main

# Send notification
# python3 ./code/notify.py --link https://github.com/Robaina/TRAITS/
