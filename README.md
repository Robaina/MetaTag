# TRAITS

This repository contains:

1. Tools to annotate short reads taxonomy and function via pylogenetic tree evolutionary placement

---
## Installation
1. Fork git repo into local machine (click on fork) and clone, or simply clone main branch with
```
git clone https://github.com/Robaina/TRAITS.git
```
2. CD to project traits and set conda environment if not already set:
```
conda env create -n traits -f environment.yml
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

## Working with GitHub branches
The best way to safely interact with the codebase is through a personal git branch. It may happen that the main branch is updated (for instace with a new script added) and one needs to integrate (merge) these new changes in their own branch without affecting the main branch. The best way to do this is by making use of the following command:

```
git pull --rebase origin main
```

or if accesing git within Visual Studio Code, one can also run the command by clicking on the Source Control options (three dots) within the Git tab and cliking on Branch / Rebase branch. All of this making sure you have selected your branch first.

__IMPORTANT NOTES__: 

1. Before pulling new updated from main, commit all the you changes you made in yout local branch
2. The command above must be exectued within your git branch, so first select the branch, the run the command.


## Obtaining reference tree

1. Preprocess database. Only run once to remove duplicates and sequences containing illegal symbols:
```
python3 code/preprocess.py --help
```

2. Make reference (protein-specific) peptide database. Run to filter origial database by given tigrfam or pfam. Also reduces redundancy of peptide database and allows filtering database sequences by minimum and/or maximum sequence length
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




