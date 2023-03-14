# Description


Assgin taxonomy and function to placed query sequences
# Epilog


Semidán Robaina Estévez (srobaina@ull.edu.es), 2021
# Usage:


```bash
usage: metatag assign [-h] [args] 

```
# Arguments

|short|long|default|help|
| :--- | :--- | :--- | :--- |
|`-h`|`--help`||show this help message and exit|
||`--jplace`|`None`|path to placements jplace file|
||`--labels`|`None`|path to label dict in pickle format. More than one space-separated path can be input|
||`--query_labels`|`None`|path to query label dict in pickle format. More than one space-separated path can be input|
||`--ref_clusters`|`None`|path to tsv file containing cluster assignment to each reference sequence id. Must contain one column named "id" and another (tab-separated) column named "cluster"|
||`--ref_cluster_scores`|`None`|path to tsv file containing cluster quality scores assigned to each cluster in the reference tree. Contains one column named "cluster" and another (tab-separated) column named "score"|
||`--outgroup`|`None`|path to text file containing IDs of sequences to be considered as an outgroup to root the tree. It can also be a fasta file from which sequence names will be extracted. It can also be a string containing a tag to filter record labels by it. The outgroup will be used to recover missing taxonomic infomation by gappa examine assign. |
||`--prefix`|`placed_tax_`|prefix to be added to output files|
||`--outdir`|`None`|path to output directory|
||`--max_placement_distance`|`None`|Maximum allowed pendant distance to consider a placement as valid. Change distance measure with parameter: "distance_measure" (defaults to pendant length)|
||`--distance_measure`|`pendant`|Choose distance measure to remove placements with distance larger than "max_placement_distance". Choose among: 1. "pendant": corresponding to pendant length of placement 2. "pendant_distal_ratio": ratio between pendant and distal distances 3. "pendant_diameter_ratio": ratio between pendant and tree diameter (largest pairwise distance) ratio. See https://github.com/lczech/gappa/wiki for a description of distal and pendant lengths.|
||`--min_placement_lwr`|`None`|Minimum allowed placement LWR to consider a placement as valid. Values between 0 and 1.|
||`--duplicated_query_ids`|`None`|path to text file containing duplicated query ids as output by seqkit rmdup|
||`--taxonomy_file`|`None`|path to tsv containing taxonomy, formated like GTDB taxopaths, for each genome ID in reference database. Defaults to None, in which case a custom GTDB taxonomy database of marine prokaryotes is used.|
