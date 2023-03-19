FROM mcr.microsoft.com/devcontainers/miniconda:0-3

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
SHELL ["conda", "run", "-n", "metatag", "/bin/bash", "-c", "poetry build && pip install dist/metatag*.whl && rm -r dist && metatag --help && conda deactivate"]
# Build and install
# RUN poetry build && pip install dist/metatag*.whl && rm -r dist && metatag --help

# Initialize conda for default user
RUN conda init
# Activate environment by default
RUN echo "conda activate metatag" >> ../root/.bashrc
RUN source ../root/.bashrc