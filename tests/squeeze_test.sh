# Test squeezemodify.py

python3 code/squeezemodify.py \
 --squeeze_table tests/test_data/squeeze_out_example.txt \
 --taxplacement_table tests/test_data/placement_example.tsv \
 --outfile tests/test_data/squeeze_out_example_replaced.txt


python3 code/squeezemodify.py \
 --squeeze_table /home/robaina/Documents/TRAITS/08.artico022_DNA.fun3.tax.noidfilter.wranks \
 --taxplacement_table tests/test_data/placement_example.tsv \
#  --outfile tests/test_data/squeeze_out_example_replaced.txt