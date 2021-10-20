#                                   Run minimal example
mkdir -p tests/test_results
# rm -r tests/test_results; mkdir tests/test_results

# Preprocess
python3 ./code/preprocess.py --in ./tests/test_data/database/ --outfile ./tests/test_results/test_data_cleaned.faa


# Make database
python3 ./code/makedatabase.py --in ./tests/test_results/test_data_cleaned.faa --outdir ./tests/test_results/ --hmm ./tests/test_data/TIGR01580.1.HMM --reduce


# Alignment and tree
python3 ./code/buildtree.py --in ./tests/test_results/ref_database.faa --outdir ./tests/test_results/ --msa_method "muscle" --tree_method "fasttree"


# Remove tree branch outliers
python3 ./code/removetreeoutliers.py --tree ./tests/test_results/ref_database.fasttree --outdir ./tests/test_results/ --aln ./tests/test_results/ref_database.faln


# Preprocess query sequences
python3 ./code/preprocess.py --in ./tests/test_data/query.faa --outfile ./tests/test_results/query_cleaned.faa --is_query


# Place query sequences
python3 code/placesequences.py --aln tests/test_results/ref_database.faln --tree tests/test_results/ref_database.fasttree --query tests/test_results/query_cleaned.faa --outdir tests/test_results/ --aln_method "papara" --tree_model "JTT" # ./tests/test_results/ref_database.model.gz
