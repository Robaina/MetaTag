# !/bin/bash

# ***************************************************************** #
                        #   Nxr pipeline                            #
# ***************************************************************** #
rm -r /home/robaina/Documents/TRAITS/genes/nxr/results
mkdir /home/robaina/Documents/TRAITS/genes/nxr/results

# Make database
python3 ./code/makedatabase.py \
 --in data/final_ref_database.fasta \
 --outdir genes/nxr/results/ \
 --hmms genes/nxr/hmms/narGTIGR01580.1.hmm \
 --relabel_prefixes "narG_" --relabel \
 --prefix "narG_" \
 --max_sizes 400

 python3 ./code/makedatabase.py \
 --in data/final_ref_database.fasta \
 --outdir genes/nxr/results/ \
 --hmms genes/nxr/hmms/Molybdopterin.hmm \
 --relabel_prefixes "molyb_" --relabel \
 --prefix "molyb_" \
 --max_sizes 400

# Merge narG, molyb and remove possible duplicated seqs
mkdir -p genes/nxr/results/molyb_narG/
mv genes/nxr/results/narG_ref_database.faa genes/nxr/results/molyb_narG/
mv genes/nxr/results/molyb_ref_database.faa genes/nxr/results/molyb_narG/

# Add nxr-classified sequences from Nitzinger, 2021
python3 ./code/preprocess.py \
 --in genes/nxr/data/Nxr_kitzinger_2021.fasta \
 --outfile genes/nxr/results/Nitzinger21_short_ids.faa \
 --idprefix "Nitz21_" --relabel

# Add nxr-classified sequences from Nitzinger, 2018
python3 ./code/preprocess.py \
 --in genes/nxr/data/Nxr_kitzinger_2018.fasta \
 --outfile genes/nxr/results/Nitzinger18_short_ids.faa \
 --idprefix "Nitz18_" --relabel

mv genes/nxr/results/Nitzinger21_short_ids.faa genes/nxr/results/molyb_narG/
mv genes/nxr/results/Nitzinger18_short_ids.faa genes/nxr/results/molyb_narG/

# Merge all four databases into final reference database
python3 ./code/preprocess.py \
 --in genes/nxr/results/molyb_narG/ \
 --outfile genes/nxr/results/molyb_narG_Nitz18_Nitz21_ref_database.faa

# Alignment and tree
python3 ./code/buildtree.py \
 --in genes/nxr/results/molyb_narG_Nitz18_Nitz21_ref_database.faa \
 --outdir genes/nxr/results/ \
 --msa_method "muscle" \
 --tree_model "LG+F+G4" \
 --tree_method "iqtree"

# Remove tree branch outliers
python3 ./code/removetreeoutliers.py \
 --tree genes/nxr/results/ref_database.newick \
 --outdir genes/nxr/results/ \
 --aln genes/nxr/results/ref_database.faln

# Relabel reference tree and msa
python3 ./code/relabeltree.py \
 --tree genes/nxr/results/ref_database_shrink.newick \
 --aln genes/nxr/results/ref_database_shrink.faln \
 --labels /home/robaina/Documents/TRAITS/genes/nxr/results/molyb_ref_database_id_dict.pickle \
          /home/robaina/Documents/TRAITS/genes/nxr/results/narG_ref_database_id_dict.pickle \
          /home/robaina/Documents/TRAITS/genes/nxr/results/Nitzinger18_short_ids_id_dict.pickle \
          /home/robaina/Documents/TRAITS/genes/nxr/results/Nitzinger21_short_ids_id_dict.pickle \
 --label_prefixes "molyb_" "narG_" "Nitz18_" "Nitz21_"

#  Commit to GitHub
git add -f /home/robaina/Documents/TRAITS/genes/nxr/results/narG_molyb_nitz18_nitz21/ref_database_shrink_relabel.contree
git commit -m "Update data"
git push origin robaina_nxr

# Send notification
python3 ./code/notify.py --link https://github.com/Robaina/TRAITS/tree/robaina_nxr/genes/nxr/results/narG_molyb_nitz18_nitz21
