# Description


MSA on reference database and infer reference tree
# Epilog


Semidán Robaina Estévez (srobaina@ull.edu.es), 2021
# Usage:


```bash
usage: metatag tree [-h] [args] 

```
# Arguments

|short|long|default|help|
| :--- | :--- | :--- | :--- |
|`-h`|`--help`||show this help message and exit|
||`--in`|`None`|path to reference database|
||`--outdir`|`None`|path to output directory|
||`--msa_method`|`muscle`|choose method for msa|
||`--tree_method`|`iqtree`|choose method for tree inference|
||`--tree_model`|`modeltest`|choose substitution model for iqtree inference.  Choices=["iqtest", "modeltest", "a valid model name"].  Defaults to optimal per modeltest-ng.|