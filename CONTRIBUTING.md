# Contributing to MetaTag

First of all, thanks for taking the time to contribute! :tada::+1:

Here you will find a set of guidelines for contributing to MetaTag. Feel free to propose changes to this document in a pull request.

## Code of conduct

This project and everyone participating in it is governed by the [Contributor Covenant, v2.0](CODE_OF_CONDUCT.md) code of conduct. By participating, you are expected to uphold this code.

## I have a question!

If you only have a question about all things related to MetaTag, the best course of actions for you is to open a new [discussion](https://github.com/Robaina/TRAITS/discussions).

## How can I contribute?

### 1. Reporting bugs

We all make mistakes, and the developers behind MetaTag are no exception... So, if you find a bug in the source code, please open an [issue](https://github.com/Robaina/TRAITS/issues) and report it. Please, first search for similar issues that are currrently open.

### 2. Suggesting enhancements

Are you missing some feature that would like MetaTag to have? No problem! You can contribute by suggesting an enhancement, just open a new issue and tag it with the [```enhancement```](https://github.com/Robaina/TRAITS/labels/enhancement) label. Please, first search for similar issues that are currrently open.

### 3. Improving the documentation

Help is always needed at improving the [documentation](https://robaina.github.io/TRAITS/). Either adding more detailed docstrings, usage explanations or new examples.

## First contribution

Unsure where to begin contributing to MetaTag? You can start by looking for issues with the label [```good first issue```](https://github.com/Robaina/TRAITS/labels/good%20first%20issue). If you are unsure about how to set a developer environment for MetaTag, do take a look at the section below. Thanks!

## Setting up a local developer environment

MetaTag depends on packages that are not available in pip, namely [HMMER](https://github.com/EddyRivasLab/hmmer) and [Prodigal](https://github.com/hyattpd/Prodigal). These can be installed from the bioconda channel. Hence, to setup up a developer environment for MetaTag:

1. Fork and download repo, cd to downloaded directory. You should create a new branch to work on your issue.

2. Create conda environment with required dependencies:

The file `envs/metatag-dev.yml` contains all dependencies required to use MetaTag. If conda is very slow solving the environment you can try installing [mamba])(https://github.com/mamba-org/mamba), just replace "conda" by "mamba" in the command below:

```bash
conda env create -n metatag-dev -f envs/metatag-dev.yml
conda activate metatag-dev
```

3. Build package

```bash
(metatag-dev) poetry build
```

4. Install MetaTag

```bash
(metatag-dev) pip install dist/metatag*.whl
```

5. Run tests

```bash
(metatag-dev) python -m unittest discover tests
```

## Working with GitHub branches
The best way to safely interact with the codebase is through a personal git branch. It may happen that the main branch is updated (for instance with a new script) and one needs to integrate (merge) these new changes in their own branch without affecting the main branch. The best way to do this is by making use of the following command:

```
git pull --rebase origin main
```

or if accesing git within Visual Studio Code, one can also run the command by clicking on the Source Control options (three dots) within the Git tab and cliking on Branch / Rebase branch. All of this making sure you have selected your branch first.

__IMPORTANT NOTES__: 

1. Before pulling new updated from main, commit all the you changes you made in yout local branch
2. The command above must be exectued within your git branch, so first select the branch, the run the command.

## Enabling GitHub's https authentication in Ubuntu/Debian
1) Create a Personal Access Token (PAT) in GitHub, copy the token (will be invisible thereafter)

2) Install GitHub CLI to securely store the token and avoid entering the token each time we need to access the CLI: 

a) https://docs.github.com/en/get-started/getting-started-with-git/caching-your-github-credentials-in-git
b) https://github.com/cli/cli/blob/trunk/docs/install_linux.md

## Using GitHub codespaces

Alternatively, you can directly work on a developer environment in the browser within GitHub's codespaces. To set this environment up:

1. Fork repository and create new branch for your issue.

2. Start a new Codespace for the new branch. GitHub will build an environment with all required dependencies as well as MetaTag the first time (it will take a couple of minutes). MetaTag will be installed in a conda environment named "metatag-dev".


## Tests on push and pull request to main

MetaTag's repo contains a [GitHub Action](https://github.com/features/actions) to perform build and integration tests which is triggered automatically on push and pull request events to the main brach. Currently the tests include building and installing MetaTag in Ubuntu and running the [test](tests) suit.