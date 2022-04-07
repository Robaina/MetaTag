# Preprocess query sequences
python3 code/preprocess.py \
 --in genes/nifH/envdata/TIGR01287_arctic.faa \
 --outfile genes/nifH/envdata/query_cleaned.faa \
 --idprefix "query_" --relabel

# Place query sequences
python3 code/placesequences.py \
 --aln genes/nifH/results/ref_database.faln \
 --tree genes/nifH/results/ref_database.newick \
 --query genes/nifH/envdata/query_cleaned.faa \
 --outdir genes/nifH/envdata/ \
 --aln_method "papara" \
--tree_model genes/nifH/results/ref_database.log

# Assign taxonomy to placed sequences
python3 code/labelplacements.py \
 --jplace genes/nifH/envdata/epa_result.jplace \
 --labels genes/nifH/results/ref_database_id_dict_clustered.pickle \
 --ref_clusters genes/nifH/data/clusters.tsv \
 --ref_cluster_scores genes/nifH/data/cluster_scores.tsv \
 --outgroup "out_" \
 --prefix "placed_tax_" \
 --outdir genes/nifH/envdata/

# Count placements (filter by taxon, cluster id and quality score)
python3 code/countplacements.py \
 --taxtable genes/nifH/envdata/placed_tax_assignments.tsv \
 --taxlevels "order" "class" "family" "genus" \
 --cluster_ids "cluster_I" "cluster_II" "cluster_III" \
 --score_threshold 0.5 \
 --outdir genes/nifH/envdata/

# Relabel tree and alignment
python3 code/relabeltree.py \
 --tree genes/nifH/envdata/epa_result.newick \
 --labels genes/nifH/results/ref_database_id_dict_clustered.pickle \
          genes/nifH/envdata/query_cleaned_id_dict.pickle \
 --taxonomy \
 --label_prefixes "ref_" "query_"