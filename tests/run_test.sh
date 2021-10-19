#                                   Run minimal example

# Preprocess
python3 ./code/preprocess.py --in ./tests/test_data/ --outfile ./tests/test_results/test_data_cleaned.faa


# Make database
python3 ./code/makedatabase.py --in ./tests/test_results/merged_data.faa --outdir ./tests/test_results/ --hmm ./tests/TIGR01580.1.HMM --reduce


# Alignment and tree
python3 ./code/buildtree.py --in ./tests/test_results/ref_database.faa --outdir ./tests/test_results/ --msa_method "muscle" --tree_method "fasttree"

# Remove tree branch outliers
python3 ./code/removetreeoutliers.py --tree ./tests/test_results/ref_database.fasttree --outdir ./tests/test_results/ --aln ./tests/test_results/ref_database.faa.aln

# Preprocess query sequences
#python3

# Place query sequences
#python3 

