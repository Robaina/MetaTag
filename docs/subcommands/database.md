# Description


Build peptide reference database
# Epilog


Semidán Robaina Estévez (srobaina@ull.edu.es), 2021
# Usage:


```bash
usage: metatag database [-h] [args] 

```
# Arguments

|short|long|default|help|
| :--- | :--- | :--- | :--- |
|`-h`|`--help`||show this help message and exit|
||`--hmms`|`None`|a single or a list of space-separated paths to tigrfam or pfam models|
||`--in`|`None`|path to peptide database|
||`--outdir`|`None`|path to output directory|
||`--prefix`|``|prefix to be added to output files|
||`--relabel_prefixes`|`None`|List of space-separated prefixes to be added to each hmm-derived set of sequences after relabelling. Only used if "--relabel" is set. Label prefixes set to "None" are assigned "ref_" by default.|
||`--max_sizes`|`None`|maximum size of representative set of sequences for each hmm model. Each (space-separated) integer corresponds to a hmm model inputed in "--hmms", thus, sorted in the same order. A value of "None" may be given to a hmm model in the list, in which case the maximum number of sequences is unlimited for that hmm.Defaults to full set of sequences for all hmm modells inputed.|
||`--min_seq_length`|`None`|minimum sequence length in reference database. Defaults to zero|
||`--max_seq_length`|`None`|maximum sequence length in reference database. Defaults to inf|
||`--relabel`||relabel record IDs with numerical ids. Unrequired to build database, but highly recommended to avoid possible conflicts downstream the pipeline.|
||`--nocdhit`||do not run cd-hit on peptide database|
||`--remove_duplicates`||remove duplicated sequences from final (merged) database|
||`--hmmsearch_args`|`None`|a string of comma-separated additional arguments for each hmm passed to hmmsearch. e.g. inputing 3 hmms: " --cut_ga --cpu 4, --cut_nc, None". IMPORTANT: the string must be preceded by a white space. A single string may be provided, in which case the same additinal arguments will be passed for each hmm. Defaults to additional arguments string: "--cut_nc". If no additional arguments are needed, provide the value "None"|
