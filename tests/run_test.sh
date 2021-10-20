# *************************************************************************************
#                                    Run minimal example
#
# Run this bash script to test whether the installation has been successfull.
# 
# The script runs the entire pipeline, from peptide sequences, gene-specific database 
# (narG in this example), reference tree reconstruction to placement of query sequences.
#
# If everything works fine, you should obtain a final figure:
# test/test_results/epa_result_tree.svg and
# a final tree: test/test_results/epa_result_tree_relabel.newick.
# 
#
# Remember that the script must be run within the "traits" conda environment, which 
# can be install from "environment.yml". See https://github.com/Robaina/TRAITS.
#
# *************************************************************************************

rm -r tests/test_results; mkdir tests/test_results

# Preprocess
python3 ./code/preprocess.py --in ./tests/test_data/database/ --outfile ./tests/test_results/test_data_cleaned.faa


# Make database
python3 ./code/makedatabase.py --in ./tests/test_results/test_data_cleaned.faa --outdir ./tests/test_results/ --hmm ./tests/test_data/TIGR01580.1.HMM --reduce


# Alignment and tree
python3 ./code/buildtree.py --in ./tests/test_results/ref_database.faa --outdir ./tests/test_results/ --msa_method "muscle" --tree_model "TEST" --tree_method "iqtree" # "fasttree"


# Remove tree branch outliers
python3 ./code/removetreeoutliers.py --tree ./tests/test_results/ref_database.contree --outdir ./tests/test_results/ --aln ./tests/test_results/ref_database.faln


# Preprocess query sequences
python3 ./code/preprocess.py --in ./tests/test_data/query.faa --outfile ./tests/test_results/query_cleaned.faa --is_query


# Place query sequences
python3 code/placesequences.py --aln tests/test_results/ref_database.faln --tree tests/test_results/ref_database.contree --query tests/test_results/query_cleaned.faa --outdir tests/test_results/ --aln_method "papara" --tree_model ./tests/test_results/ref_database.log # JTT
