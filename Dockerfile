FROM mcr.microsoft.com/devcontainers/miniconda:0-3

WORKDIR /TRAITS
# Copy repo to docker container
COPY code code/
COPY tests tests/
COPY envs envs/
COPY README.md .
# COPY pyproject.toml .
COPY LICENSE .

# Make conda environment and activate
# RUN conda install mamba -n base -c conda-forge
RUN conda env create -f envs/traits-test.yml
SHELL ["conda", "run", "-n", "traits", "/bin/bash", "-c"]
# Build and install Pynteny
# RUN poetry build && pip install dist/pynteny*.whl && pynteny --help
# Give read/write permissions to install directory (needed to make config.json)
RUN python --version
RUN chmod ugo+rw /opt/conda/envs/traits/lib/python3.8/site-packages

# Initialize conda for default user
RUN conda init
# # Activate pynteny environment by default
RUN echo "conda activate traits" >> ../root/.bashrc
RUN source ../root/.bashrc