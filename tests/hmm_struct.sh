# rm -r tests/test_results; mkdir tests/test_results

# Make database
python3 code/makedatabase_struct.py \
 --in /home/robaina/Documents/TRAITS/data/marhq_cleaned.faa \
 --outdir /home/robaina/Documents/TRAITS/tests/test_results/ \
 --hmms /home/robaina/Documents/tigrfams/hmm_PGAP/TIGR01282.1.HMM /home/robaina/Documents/tigrfams/hmm_PGAP/TIGR01286.1.HMM /home/robaina/Documents/tigrfams/hmm_PGAP/TIGR01287.1.HMM \
 --hmm_struc "TIGR01287.1 10 TIGR01282.1 10 <TIGR01286.1" \
 --max_size 800