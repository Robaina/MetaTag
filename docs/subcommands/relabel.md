# Description


Relabel tree and msa based on input label dictionaries
# Epilog


Semidán Robaina Estévez (srobaina@ull.edu.es), 2021
# Usage:


```bash
usage: metatag relabel [-h] [args] 

```
# Arguments

|short|long|default|help|
| :--- | :--- | :--- | :--- |
|`-h`|`--help`||show this help message and exit|
||`--tree`|`None`|path to tree in newick format|
||`--labels`|`None`|path to label dict in pickle format. More than one space-separated path can be input|
||`--label_prefixes`|`None`|prefix(es) to be added to sequences in each label dict,input in same order as labels.More than one space-separated prefix can be input|
||`--taxonomy`||assign GTDB taxonomy to labels containing MMP identifiers|
||`--aln`|`None`|path to fasta alignment file to be relabelled|
||`--outdir`|`None`|path to output directory|
||`--taxonomy_file`|`None`|path to tsv containing taxonomy, formated like GTDB taxopaths, for each genome ID in reference database. Defaults to None, in which case a custom GTDB taxonomy database of marine prokaryotes is used.|
