#! /bin/bash

# ************************************************************************************** #
#                                  Run minimal example                                   #
#                                                                                        #
# Run this bash script to test whether the installation has been successfull.            #
#                                                                                        #
# The script runs the entire pipeline, from peptide sequences, gene-specific database    #
# (narG in this example), reference tree reconstruction to placement of query sequences. #
#                                                                                        #
# ************************************************************************************** #

mkdir -p tests/test_results; mkdir -p tests/test_results/gappa/

# Preprocess
metatag preprocess \
 --in tests/test_data/database/ \
 --outfile tests/test_results/test_data_cleaned.faa \
 --export-duplicates

# Build reference database
metatag database \
 --in tests/test_results/test_data_cleaned.faa \
 --outdir tests/test_results \
 --hmms tests/test_data/TIGR01287.1.HMM \
        tests/test_data/TIGR02016.1.HMM \
 --max_sizes 20 5 \
 --min_seq_length 10 --max_seq_length 2000 \
 --relabel_prefixes "ref_" "out_" \
 --relabel --remove_duplicates \
 --hmmsearch_args " None, --cut_ga"

# Alignment and tree
metatag tree \
 --in tests/test_results/ref_database.faa \
 --outdir tests/test_results/ \
 --msa_method "muscle" \
 --tree_model "iqtest" \
 --tree_method "fasttree"

# Relabel reference tree and assign taxonomy
metatag relabel \
 --tree tests/test_results/ref_database.newick \
 --labels tests/test_results/ref_database_id_dict.pickle \
 --taxonomy

# Preprocess query sequences
metatag preprocess \
 --in tests/test_data/query.faa \
 --outfile tests/test_results/query_cleaned.faa \
 --idprefix "query_" --relabel --export-duplicates

# Place query sequences
metatag place \
 --aln tests/test_results/ref_database.faln \
 --tree tests/test_results/ref_database.newick \
 --query tests/test_results/query_cleaned.faa \
 --outdir tests/test_results/ \
 --aln_method "papara" \
 --tree_model "JTT"
 
# Assign taxonomy to placed sequences
metatag assign \
 --jplace tests/test_results/epa_result.jplace \
 --labels tests/test_results/ref_database_id_dict.pickle \
 --query_labels tests/test_results/query_cleaned_id_dict.pickle \
 --duplicated_query_ids tests/test_results/query_cleaned_duplicates.txt \
 --ref_clusters tests/test_data/clusters.tsv \
 --ref_cluster_scores tests/test_data/cluster_scores.tsv \
 --prefix "placed_tax_" \
 --outdir tests/test_results/gappa/ \
 --max_placement_distance 1.0 --distance_measure "pendant_diameter_ratio" \
 --min_placement_lwr 0.8

# Count placements (filter by taxon, cluster id and quality score)
metatag count \
 --taxtable tests/test_results/gappa/placed_tax_assignments.tsv \
 --taxlevels "genus" "family" "class" "order" \
 --cluster_ids "G1" "G2" \
 --score_threshold 0.6 \
 --outdir tests/test_results/gappa/ \
 --export_right_queries

# Relabel tree with placements 
# (note: tree retrieved with gappa examine graft, which selects placements with highest LWR)
metatag relabel \
 --tree tests/test_results/epa_result.newick \
 --labels tests/test_results/ref_database_id_dict.pickle \
          tests/test_results/query_cleaned_id_dict.pickle \
 --label_prefixes "ref_" "query_" --taxonomy