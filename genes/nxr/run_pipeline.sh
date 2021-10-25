# ***************************************************************** #
#                           Nxr pipeline                            #
# ***************************************************************** #

# Preprocess
#python3 ./code/preprocess.py \
#  --in genes/nxr/data/MAR_HQ \
#  --outfile genes/nxr/nxr_results/data_cleaned.faa

# Make database
python3 ./code/makedatabase.py \
 --in genes/nxr/nxr_results/data_cleaned.faa \
 --outdir genes/nxr/nxr_results/ \
 --hmm genes/nxr/hmms/Molybdopterin.hmm \
 --prefix "molyb_" \
 --max_size 300

# python3 ./code/makedatabase.py \
#  --in genes/nxr/nxr_results/data_cleaned.faa \
#  --outdir genes/nxr/nxr_results/ \
#  --hmm genes/nxr/hmms/narGTIGR01580.1.hmm \
#  --prefix "narG_" \
#  --max_size 300

# Merge narG, molyb and remove possible duplicated seqs
# mkdir nxr_results/molyb_narG/
# mv nxr_results/narG_ref_database.faa nxr_results/molyb_narG/narG_ref_database.faa
# mv nxr_results/molyb_ref_database.faa nxr_results/molyb_narG/molyb_ref_database.faa

# python3 ./code/preprocess.py \
#  --in genes/nxr/nxr_results/molyb_narG/ \
#  --outfile genes/nxr/nxr_results/molyb_narG_ref_database.faa

# Add nxr-classified sequences from Nitzinger, 2021




# # Alignment and tree
# python3 ./code/buildtree.py \
#  --in ./tests/test_results/ref_database.faa \
#  --outdir ./tests/test_results/ \
#  --msa_method "muscle" \
#  --tree_model "TEST" \
#  --tree_method "iqtree"

# # Remove tree branch outliers
# python3 ./code/removetreeoutliers.py \
#  --tree ./tests/test_results/ref_database.contree \
#  --outdir ./tests/test_results/ \
#  --aln ./tests/test_results/ref_database.faln

# # Preprocess query sequences
# python3 ./code/preprocess.py \
#  --in ./tests/test_data/query.faa \
#  --outfile ./tests/test_results/query_cleaned.faa \
#  --is_query

# # Place query sequences
# python3 code/placesequences.py \
#  --aln tests/test_results/ref_database.faln \
#  --tree tests/test_results/ref_database.contree \
#  --query tests/test_results/query_cleaned.faa \
#  --outdir tests/test_results/ \
#  --aln_method "hmmalign" \
#  --tree_model ./tests/test_results/ref_database.log \
#  --open_in_browser
 
 # ./tests/test_results/ref_database.log # JTT