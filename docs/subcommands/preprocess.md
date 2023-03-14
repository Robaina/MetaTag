# Description


Database preprocessing: removal of duplicated sequences and of sequences with illegal symbols. To preferentially keep one duplicate sequence over another, place preferred sequences first.
# Epilog


Semidán Robaina Estévez (srobaina@ull.edu.es), 2021
# Usage:


```bash
usage: metatag preprocess [-h] [args] 

```
# Arguments

|short|long|default|help|
| :--- | :--- | :--- | :--- |
|`-h`|`--help`||show this help message and exit|
||`--in`|`None`|path to fasta file or directory containing fasta files|
||`--dna`||declare if sequences are nucleotides. Defaults to peptide sequences.|
||`--translate`||choose whether nucleotide sequences are translated with prodigal|
||`--export-duplicates`||choose whether duplicated sequences are exported to file (same directory than outfile)|
||`--outfile`|`None`|path to output fasta file|
||`--relabel`||relabel record IDs with numeral ids|
||`--relabel_prefix`|`None`|prefix to be added to sequence IDs|