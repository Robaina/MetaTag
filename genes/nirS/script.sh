python3 ./code/makedatabase.py \
 --in data/final_ref_database.fasta \
 --outdir genes/nirS/results \
 --hmms genes/nirS/hmms/K15864.hmm \
        genes/nirS/hmms/PF02239.hmm \
 --hmmsearch_args " -T 413.67, None" \
 --prefix "nirS_DB_" \
 --relabel \
 --relabel_prefixes "refKO_" "refPF_" \
 --remove_duplicates \
 --max_sizes 800