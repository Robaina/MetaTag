# TRAITS

Repository of project TRAITS.

---
## Installation
1. Fork git repo into local machine (click on fork) and clone, or simply clone main branch with
```
git clone https://github.com/Robaina/TRAITS.git
```
2. CD to project traits and set conda environment if not already set:
```
conda env create -f environment.yml
```
3. [Install papara](https://cme.h-its.org/exelixis/web/software/papara/index.html) executable and add path to executable in module phyloplacement.wrappers.py:
```python
"""
Simple CLI wrappers to several tools
"""

import os
import shutil
import tempfile
from phyloplacement.utils import terminalExecute, setDefaultOutputPath

papara_exec = 'path/to/papara.exe'
```

## Installation test
To check that everything is working properly, you can run a test what will perform the entire workflow on a minimal dataset. To run it, call bash script:
```
bash tests/run_test.sh
```
It should produce a final tree with query sequences placed on it, as well as a bunch of intermediary files without any errors.

## Obtaining reference tree

1. Preprocess database. Only run once to remove duplicates and sequences containing illegal symbols:
```
python3 code/preprocess.py --help
```

2. Make reference (protein-specific) peptide database. Run to filter origial database by given tigrfam or pfam. Also reduces redundancy of peptide database
```
python3 code/makedatabase.py --help
```

3. Perform multiple sequence alignment (MSA) on reference database and infer reference tree
```
python3 code/buildtree.py --help
```
## Curating reference tree

* Remove tree and reference msa outliers
```
python3 code/removetreeoutliers.py --help
```
* Manual curation and classification of clusters

## Placement of query short reads on reference tree

1. Preprocess query sequences. Remove illegal symbols from nucleotide or peptide sequence (if already translated). Can also translate nucleotide sequence with prodigal if not already translated.
```
python3 code/preprocess.py --help
```
2. Place sequences onto tree. Run either papara or hmmalign to align query sequences to reference MSA. Call epa-ng to place sequences onto tree. Call gappa to make new tree newick file with placements.
```
python3 code/placesequences.py --help
```




