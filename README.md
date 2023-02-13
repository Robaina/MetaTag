# TRAITS
[![tests](https://github.com/Robaina/TRAITS/actions/workflows/tests.yml/badge.svg)](https://github.com/Robaina/TRAITS/actions/workflows/tests.yml)
![python](https://img.shields.io/badge/Python-3.8%2B-blue)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![license](https://img.shields.io/github/license/Robaina/Pynteny)
![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4)

This repository contains tools to assign taxonomy and function annotations to short reads through pylogenetic tree evolutionary placement

---
## Setup
The easiest way to use TRAITS is through the provided [docker](https://www.docker.com/) container. To use it, pull the image:

```
docker pull ghcr.io/robaina/traits:latest
```

Then run the container interactively:

```
docker run -i ghcr.io/robaina/traits:latest
```

Otherwise, you can install TRAITS like this:

1. Fork git repo into local machine (click on fork) and clone, or simply clone main branch with
```
git clone https://github.com/Robaina/TRAITS.git
```
2. CD to project traits and set conda environment if not already set:
```
conda env create -n traits -f envs/metatag.yml
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
python3 src/preprocess.py --help
```

2. Make reference (protein-specific) peptide database. Run to filter origial database by given tigrfam or pfam. Also reduces redundancy of peptide database and allows filtering database sequences by minimum and/or maximum sequence length
```
python3 src/makedatabase.py --help
```

3. Perform multiple sequence alignment (MSA) on reference database and infer reference tree
```
python3 src/buildtree.py --help
```
## Curating reference tree

* Remove tree and reference msa outliers
```
python3 src/removetreeoutliers.py --help
```
* Manual curation and classification of clusters

## Placement of query short reads on reference tree

1. Preprocess query sequences. Remove illegal symbols from nucleotide or peptide sequence (if already translated). Can also translate nucleotide sequence with prodigal if not already translated.
```
python3 src/preprocess.py --help
```
2. Place sequences onto tree. Run either papara or hmmalign to align query sequences to reference MSA. Call epa-ng to place sequences onto tree. Call gappa to make new tree newick file with placements.
```
python3 src/placesequences.py --help
```
3. Assign taxonomy and function to placed sequences. Assign taxonomy and function to placed query sequences.
```
python3 src/labelplacements.py --help
```
4. Quantify placed queries. Filter placed sequences according to cluster/placement quality score and quantify remaining placed sequences according to function and taxonomy
```
python3 src/countplacements.py --help
```




