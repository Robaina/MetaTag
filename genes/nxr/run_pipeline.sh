# ***************************************************************** #
#                           Nxr pipeline                            #
# ***************************************************************** #

# Preprocess
#python3 ./code/preprocess.py \
#  --in genes/nxr/data/MAR_HQ \
#  --outfile genes/nxr/nxr_results/data_cleaned.faa

# Make database
# python3 ./code/makedatabase.py \
#  --in genes/nxr/nxr_results/data_cleaned.faa \
#  --outdir genes/nxr/nxr_results/ \
#  --hmm genes/nxr/hmms/narGTIGR01580.1.hmm \
#  --prefix "narG_" --relabel \
#  --max_size 300

# python3 ./code/makedatabase.py \
#  --in genes/nxr/nxr_results/data_cleaned.faa \
#  --outdir genes/nxr/nxr_results/ \
#  --hmm genes/nxr/hmms/Molybdopterin.hmm \
#  --prefix "molyb_" --relabel \
#  --max_size 300

# Merge narG, molyb and remove possible duplicated seqs
# mkdir -p genes/nxr/nxr_results/molyb_narG/
# mv genes/nxr/nxr_results/narG_ref_database.faa genes/nxr/nxr_results/molyb_narG/
# mv genes/nxr/nxr_results/molyb_ref_database.faa genes/nxr/nxr_results/molyb_narG/

# Add nxr-classified sequences from Nitzinger, 2021
# python3 ./code/preprocess.py \
#  --in genes/nxr/data/Nxr_kitzinger_2021.fasta \
#  --outfile genes/nxr/nxr_results/Nitzinger_short_ids.faa \
#  --idprefix "ref_Nitz_" --relabel

# mv genes/nxr/nxr_results/Nitzinger_short_ids.faa genes/nxr/nxr_results/molyb_narG/

# # Merge all three databases into final reference database
# python3 ./code/preprocess.py \
#  --in genes/nxr/nxr_results/molyb_narG/ \
#  --outfile genes/nxr/nxr_results/molyb_narG_Nitz_ref_database.faa

# Alignment and tree
python3 ./code/buildtree.py \
 --in genes/nxr/nxr_results/molyb_narG_Nitz_ref_database.faa \
 --outdir genes/nxr/nxr_results/ \
 --msa_method "muscle" \
 --tree_model "TEST" \
 --tree_method "iqtree"

# Remove tree branch outliers
python3 ./code/removetreeoutliers.py \
 --tree genes/nxr/nxr_results/molyb_narG_Nitz_ref_database.contree \
 --outdir genes/nxr/nxr_results/ \
 --aln genes/nxr/nxr_results/molyb_narG_Nitz_ref_database.faln

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