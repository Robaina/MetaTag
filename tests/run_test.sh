# ************************************************************************************** #
#                                  Run minimal example                                   #
#                                                                                        #
# Run this bash script to test whether the installation has been successfull.            #
#                                                                                        #
# The script runs the entire pipeline, from peptide sequences, gene-specific database    #
# (narG in this example), reference tree reconstruction to placement of query sequences. #
#                                                                                        #
# If everything works fine, you should obtain a final figure:                            #
# test/test_results/epa_result_tree.svg and                                              #
# a final tree: test/test_results/epa_result_tree_relabel.newick.                        #
#                                                                                        #
#                                                                                        #
# Remember that the script must be run within the "traits" conda environment, which      #
# can be installed from "environment.yml". See https://github.com/Robaina/TRAITS.        #
#                                                                                        #
# ************************************************************************************** #

rm -r tests/test_results; mkdir tests/test_results; mkdir tests/test_results/gappa/

# Preprocess
python3 code/preprocess.py \
 --in tests/test_data/database/ \
 --outfile tests/test_results/test_data_cleaned.faa \
 --export-duplicates

# Make database
python3 code/makedatabase.py \
 --in tests/test_results/test_data_cleaned.faa \
 --outdir tests/test_results \
 --hmms tests/test_data/TIGR01287.1.HMM \
        tests/test_data/TIGR02016.1.HMM \
 --max_sizes 20 None \
 --min_seq_length 10 --max_seq_length 2000 \
 --relabel_prefixes "ref_" "out_" \
 --relabel

# # Alignment and tree
# python3 code/buildtree.py \
#  --in tests/test_results/ref_database.faa \
#  --outdir tests/test_results/ \
#  --msa_method "muscle" \
#  --tree_model "iqtest" \
#  --tree_method "iqtree"

# # Remove tree branch outliers
# python3 code/removetreeoutliers.py \
#  --tree tests/test_results/ref_database.newick \
#  --outdir tests/test_results/ \
#  --aln tests/test_results/ref_database.faln

# # Preprocess query sequences
# python3 code/preprocess.py \
#  --in tests/test_data/query.faa \
#  --outfile tests/test_results/query_cleaned.faa \
#  --idprefix "query_" --relabel

# # Place query sequences
# python3 code/placesequences.py \
#  --aln tests/test_results/ref_database_shrink.faln \
#  --tree tests/test_results/ref_database_shrink.newick \
#  --query tests/test_results/query_cleaned.faa \
#  --outdir tests/test_results/ \
#  --aln_method "papara" \
# --tree_model tests/test_results/ref_database.log

# # Extract outgroup IDs for gappa examine assign (temporary fix)
# grep "out_" tests/test_results/ref_database.faa | cut -c 2- > tests/test_results/outgroup_ids.txt

# # Assign taxonomy to placed sequences
# python3 code/labelplacements.py \
#  --jplace tests/test_results/epa_result.jplace \
#  --labels tests/test_results/ref_database_id_dict.pickle \
#  --ref_clusters tests/test_data/clusters.tsv \
#  --ref_cluster_scores tests/test_data/cluster_scores.tsv \
#  --outgroup tests/test_results/outgroup_ids.txt \
#  --prefix "placed_tax_" \
#  --outdir tests/test_results/gappa/

# # Count placements (filter by taxon, cluster id and quality score)
# python3 code/countplacements.py \
#  --taxtable tests/test_results/gappa/placed_tax_assignments.tsv \
#  --taxlevel "family" \
#  --cluster_ids "G1" "G2" \
#  --score_threshold 0.6 \
#  --outfile tests/test_results/gappa/placed_family_tax_counts.tsv \
#  --outpdf tests/test_results/gappa/placed_family_tax_counts.pdf

# # Relabel tree and alignment
# python3 code/relabeltree.py \
#  --tree tests/test_results/epa_result.newick \
#  --labels tests/test_results/ref_database_id_dict.pickle \
#           tests/test_results/query_cleaned_id_dict.pickle \
#  --label_prefixes "ref_" "query_" \
#  --taxonomy