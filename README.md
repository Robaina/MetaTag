![logo](docs/imgs/meta_trans.png)

[![tests](https://github.com/Robaina/TRAITS/actions/workflows/tests.yml/badge.svg)](https://github.com/Robaina/TRAITS/actions/workflows/tests.yml)
![python](https://img.shields.io/badge/Python-3.8%2B-blue)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![license](https://img.shields.io/github/license/Robaina/Pynteny)
![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4)

---
## 1. :bulb: What is MetaTag?
This repository contains tools to assign taxonomy and function annotations to short reads through pylogenetic tree evolutionary placement

## :wrench: Setup
The easiest way to use TRAITS is through the provided [docker](https://www.docker.com/) container. To use it, pull the image:

```
docker pull ghcr.io/robaina/metatag:latest
```

Then run the container interactively:

```
docker run -i ghcr.io/robaina/metatag:latest
```

Otherwise, you can install MetaTag like this:

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


## :arrows_counterclockwise: Dependencies
Pynteny would not work without these awesome projects:

- [hmmer](https://github.com/EddyRivasLab/hmmer)
- [prodigal](https://github.com/hyattpd/Prodigal)
- [pyfastx](https://github.com/lmdu/pyfastx)
- [biopython](https://github.com/biopython/biopython)
- [numpy](https://github.com/numpy/numpy)
- [pandas](https://github.com/pandas-dev/pandas)
- [psutil](https://github.com/giampaolo/psutil)

Thanks!

## :octocat: Contributing

Contributions are always welcome! If you don't know where to start, you may find an interesting [issue to work in here](https://github.com/Robaina/TRAITS/issues). Please, read our [contribution guidelines](CONTRIBUTING.md) first.