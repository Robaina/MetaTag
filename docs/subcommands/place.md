# Description


Place query sequences onto reference tree
# Epilog


Semidán Robaina Estévez (srobaina@ull.edu.es), 2021
# Usage:


```bash
usage: metatag place [-h] [args] 

```
# Arguments

|short|long|default|help|
| :--- | :--- | :--- | :--- |
|`-h`|`--help`||show this help message and exit|
||`--aln`|`None`|path to reference fasta alignment|
||`--tree`|`None`|path to reference tree|
||`--query`|`None`|path to query peptide sequences.  Query sequences should be already preprocessed to handle illegal symbols|
||`--tree_model`|`None`|provide subsitution model employed to infer tree. Can be: 1) a valid model name or 2) a path to the log file returned by iqtree|
||`--outdir`|`None`|path to output directory|
||`--aln_method`|`papara`|choose method to align query sequences to reference alignment|