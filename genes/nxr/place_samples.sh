# Place environmental samples onto NXR tree

# Preprocess query sequences
# python3 code/preprocess.py \
#  --in ../../remedios_merged.fasta \
#  --outfile query_cleaned.faa \
#  --idprefix "query_" --dna --translate --relabel

# Place query sequences
# python3 code/placesequences.py \
#  --aln genes/nxr/results/ref_database_shrink.faln \
#  --tree genes/nxr/results/ref_database_shrink.contree \
#  --query genes/nxr/results/remedios_placement_results/query_cleaned.faa \
#  --outdir genes/nxr/results/remedios_placement_results/ \
#  --aln_method "papara" \
#  --tree_model genes/nxr/results/ref_database.log

# Assign taxonomy to placed sequences
# python3 code/labelplacements.py \
#  --jplace genes/nxr/results/remedios_placement_results/epa_result.jplace \
#  --labels genes/nxr/results/remedios_placement_results/molyb_ref_database_id_dict.pickle \
#           genes/nxr/results/remedios_placement_results/narG_ref_database_id_dict.pickle \
#           genes/nxr/results/remedios_placement_results/query_cleaned_id_dict.pickle \
#  --outgroup genes/nxr/results/remedios_placement_results/data/outliers_short_ids.faa \
#  --prefix "test_placed_tax_" \
#  --outdir genes/nxr/results/remedios_placement_results/gappa/

# # Relabel tree and alignment
# python3 code/relabeltree.py \
#  --tree genes/nxr/results/remedios_placement_results/epa_result.newick \
#  --labels genes/nxr/results/remedios_placement_results/test_ref_database_id_dict.pickle \
#           genes/nxr/results/remedios_placement_results/outliers_short_ids_id_dict.pickle \
#           genes/nxr/results/remedios_placement_results/query_cleaned_id_dict.pickle \
#  --label_prefixes "ref_" "out_" "query_" \
#  --taxonomy