![logo](docs/imgs/meta_trans.png)

[![tests](https://github.com/Robaina/MetaTag/actions/workflows/tests.yml/badge.svg)](https://github.com/Robaina/MetaTag/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/Robaina/MetaTag/branch/main/graph/badge.svg?token=PY48LGP84S)](https://codecov.io/gh/Robaina/MetaTag)

![python](https://img.shields.io/badge/Python-3.10-blue)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

![license](https://img.shields.io/github/license/Robaina/MetaTag)
![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4)

---
## :bulb: What is MetaTag?
This repository contains tools to assign taxonomy and function annotations to short reads through pylogenetic tree evolutionary placement

## :wrench: Setup
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
### Installation test
To check that everything is working properly, you can run a test what will perform the entire workflow on a minimal dataset. To run it, call bash script:
```bash
conda activate metatag
(metatag) bash tests/run_test.sh
```
It should produce a final tree with query sequences placed on it, as well as a bunch of intermediary files without any errors.


## :arrows_counterclockwise: Dependencies
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

## :octocat: Contributing

Contributions are always welcome! If you don't know where to start, you may find an interesting [issue to work in here](https://github.com/Robaina/MetaTag/issues). Please, read our [contribution guidelines](CONTRIBUTING.md) first.
