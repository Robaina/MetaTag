FROM mcr.microsoft.com/devcontainers/miniconda:0-3
USER root

WORKDIR /MetaTag
# Copy repo to docker container
COPY src src/
COPY tests tests/
COPY envs envs/
COPY README.md .
COPY pyproject.toml .
COPY LICENSE .

# Make conda environment and activate
RUN conda env create -n metatag -f envs/metatag.yml
SHELL ["conda", "run", "-n", "metatag", "/bin/bash", "-c"]
# Build and install
RUN poetry build && pip install dist/metatag*.whl && rm -r dist && metatag --help

# WHEN run locally as root, use this:
# Initialize conda for default user: root
RUN conda init
# Activate environment by default
RUN echo "conda activate metatag" >> ../root/.bashrc
RUN source ../root/.bashrc

# WHEN run in GitHub Action, use this:
# Initialize conda for default user: vscode
# RUN runuser -l vscode -c 'conda init'
# # Activate environment by default
# RUN echo "conda activate metatag" >> /home/vscode/.bashrc
# RUN source /home/vscode/.bashrc