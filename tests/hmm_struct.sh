# rm -r tests/test_results; mkdir tests/test_results

# Make database
python3 src/makedatabase_struct.py \
 --in /home/robaina/Documents/TRAITS/data/marhq_cleaned.faa \
 --outdir /home/robaina/Documents/TRAITS/tests/test_results/ \
 --hmms /home/robaina/Documents/tigrfams/hmm_PGAP/TIGR01282.1.HMM /home/robaina/Documents/tigrfams/hmm_PGAP/TIGR01286.1.HMM /home/robaina/Documents/tigrfams/hmm_PGAP/TIGR01287.1.HMM \
 --target_hmm "TIGR01287.1" \
 --hmm_struc "TIGR01287.1 10 TIGR01282.1 10 <TIGR01286.1" \
 --max_size 800 \
 --prefix "TIGR01287.1_" \
 --hmmsearch_args "--cut_nc, None, --cut_nc"  # or comma-separated arguments for each hmm. e.g.: "--cut_nc, --cut_nc, --cut_ga"