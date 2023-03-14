![logo](https://user-images.githubusercontent.com/21340147/223135028-afba2744-767f-4275-a60e-93d0cf957a70.png)

MetaTag is a Python-based pipeline to assign functional and taxonomical annotations to unannotated short metagenomic sequence data sets.

These are the available subcommands, run as ```metatag <subcommand> <options>```:

- [preprocess](subcommands/preprocess.md)
- [database](subcommands/database.md)
- [tree](subcommands/tree.md)
- [place](subcommands/place.md)
- [assign](subcommands/assign.md)
- [count](subcommands/count.md)
- [plot](subcommands/plot.md)
- [relabel](subcommands/relabel.md)

## Setup

The easiest way to use MetaTag is through the provided [docker](https://www.docker.com/) container. To use it, pull the image:

```bash
docker pull ghcr.io/robaina/metatag:latest
```

Then run the container interactively:

```bash
docker run -i ghcr.io/robaina/metatag:latest
```

Otherwise, you can install MetaTag like this:

1. Fork git repo into local machine (click on fork) and clone, or simply clone main branch with
```bash
git clone https://github.com/Robaina/MetaTag.git
```

2. CD to project MetaTag and set conda environment if not already set:
```bash
conda env create -n metatag -f envs/metatag-dev.yml
```

3. Build and install MetaTag:
```bash
conda activate metatag
(metatag) poetry build && pip install dist/metatag*.whl
```

### Installing MetaTag on Windows

MetaTag is designed to run on Linux machines. However, it can be installed within the [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install) via conda or docker.

## General usage


## Getting started and Examples

## Dependencies

MetaTag would not work without these awesome projects:

- [hmmer](https://github.com/EddyRivasLab/hmmer)
- [prodigal](https://github.com/hyattpd/Prodigal)
- [pyfastx](https://github.com/lmdu/pyfastx)
- [biopython](https://github.com/biopython/biopython)
- [papara](https://cme.h-its.org/exelixis/web/software/papara/index.html)
- [epa-ng](https://github.com/Pbdas/epa-ng)
- [gappa](https://github.com/lczech/gappa)
- [empress](https://github.com/biocore/empress)
- [muscle](https://github.com/biocore/empress)
- [mafft](https://github.com/biocore/empress)
- [numpy](https://github.com/numpy/numpy)
- [pandas](https://github.com/pandas-dev/pandas)
- [psutil](https://github.com/giampaolo/psutil)

Thanks!


## Contributing

Contributions are always welcome! If you don't know where to start, you may find an interesting [issue to work in here](https://github.com/Robaina/MetaTag/issues). Please, read our [contribution guidelines](https://github.com/Robaina/MetaTag/blob/main/CONTRIBUTING.md) first.

## Citation

If you use this software, please cite it as below:

Semidán Robaina Estévez. (2022). MetaTag: (Version 0.0.2). Metagenome functional and taxonomical annotation through phylogenetic tree placement. Zenodo. https://doi.org/10.5281/zenodo.7048685