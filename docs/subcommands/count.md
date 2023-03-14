# Description


Count placed sequences based on taxon level, function and quality score
# Epilog


Semidán Robaina Estévez (srobaina@ull.edu.es), 2021
# Usage:


```bash
usage: metatag count [-h] [args] 

```
# Arguments

|short|long|default|help|
| :--- | :--- | :--- | :--- |
|`-h`|`--help`||show this help message and exit|
||`--taxtable`|`None`|path to placements taxonomy table file|
||`--taxlevels`|`None`|specify space-separated tax levels to count hits|
||`--cluster_ids`|`None`|list of space-separated target cluster ids of the reference tree corresponding to the selected function to filter counts. If not provided, then all clusters in the tree are considered for counting.|
||`--score_threshold`|`None`|cluster score threshold value to filter placement results|
||`--outdir`|`None`|path to output directory|
||`--prefix`|`None`|prefix to be added to output files|
||`--export_right_queries`||export a tsv file containing queries with cluster assginments that passed provided filters (--score_threshold and/or --cluster_ids)|
