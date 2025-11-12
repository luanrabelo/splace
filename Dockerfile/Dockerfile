FROM ubuntu:latest

LABEL maintainer="Luan Rabelo and Renato Oliveira"

# Non-interactive apt and basic locale
ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PATH=/opt/conda/bin:$PATH

RUN set -eux; \
    apt-get update; \
    apt-get dist-upgrade -y; \
    apt-get install -y --no-install-recommends \
      wget bzip2 ca-certificates curl git; \
    rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN set -eux; \
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh; \
    bash /tmp/miniconda.sh -b -p /opt/conda; \
    rm -f /tmp/miniconda.sh; \
    conda config --system --set always_yes yes; \
    conda config --system --set changeps1 no; \
    conda update -n base -c defaults conda

# Create the conda environment from environment.yml
ARG ENV_NAME=splace_env
WORKDIR /tmp
COPY environment.yml /tmp/environment.yml
RUN set -eux; \
    conda env create -f /tmp/environment.yml; \
    conda clean -afy

# Put env on the PATH for subsequent layers
ENV CONDA_DEFAULT_ENV=${ENV_NAME}
ENV PATH=/opt/conda/envs/${ENV_NAME}/bin:/opt/conda/bin:$PATH

# Copy project
WORKDIR /app
COPY . /app

# Run SPLACE using the conda environment
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "splace_env", "python", "splace.py"]
CMD ["--help"]