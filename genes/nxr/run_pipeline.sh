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
#  --hmm genes/nxr/hmms/Molybdopterin.hmm \
#  --prefix "molyb_" \
#  --reduce

# python3 ./code/makedatabase.py \
#  --in genes/nxr/nxr_results/data_cleaned.faa \
#  --outdir genes/nxr/nxr_results/ \
#  --hmm genes/nxr/hmms/narGTIGR01580.1.hmm \
#  --prefix "narG_" \
#  --reduce

# Merge narG, molyb and remove possible duplicated seqs
mkdir genes/nxr/nxr_results/molyb_narG/
mv genes/nxr/nxr_results/molyb_ref_database.faa genes/nxr/nxr_results/molyb_narG/
mv genes/nxr/nxr_results/narG_ref_database.faa genes/nxr/nxr_results/molyb_narG/

python3 ./code/preprocess.py \
 --in genes/nxr/nxr_results/molyb_narG \
 --outfile genes/nxr/nxr_results/molyb_narG_ref_database.faa

# Run repset on narG and molyb reference databases (these databases have been prevously cd-hited, and hmmsearch with --cut_nc)
source /home/ubuntu/anaconda3/etc/profile.d/conda.sh
conda deactivate; conda activate repset



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