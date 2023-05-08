FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="2e9fee8dc236c82b61942a3c9bec4108b093e9d2d907ef6a99fe8b732fd04ac5"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/R.yaml
#   prefix: /conda-envs/50eb9a3e12a1bbfafcbb643868be3ca4
#   name: R
#   channels:
#     - conda-forge
#   dependencies:
#     - r-rmarkdown ==2.14
#     - r-tibble ==3.2.1
#     - r-ggplot2 ==3.4.1
#     - r-heatmaply ==1.4.2
#     - r-ggdendro ==0.1.23
#     - r-reshape2 ==1.4.4
#     - r-RColorBrewer ==1.1_3
#     - r-cowplot ==1.1.1
RUN mkdir -p /conda-envs/50eb9a3e12a1bbfafcbb643868be3ca4
COPY workflow/envs/R.yaml /conda-envs/50eb9a3e12a1bbfafcbb643868be3ca4/environment.yaml

# Conda environment:
#   source: workflow/envs/bash.yaml
#   prefix: /conda-envs/92c542c01e3b7bfb410fd03a54e00875
#   name: bash
#   channels:
#     - conda-forge
#   dependencies:
#     - python ==3.11.0
RUN mkdir -p /conda-envs/92c542c01e3b7bfb410fd03a54e00875
COPY workflow/envs/bash.yaml /conda-envs/92c542c01e3b7bfb410fd03a54e00875/environment.yaml

# Conda environment:
#   source: workflow/envs/bwa-mem2.yaml
#   prefix: /conda-envs/6e8595c2869a59bcfb9d3251dcc7705f
#   name: bwa-mem2
#   channels:
#     - bioconda
#   dependencies:
#     - bwa-mem2 ==2.2.1
#     - samtools ==1.16.1
RUN mkdir -p /conda-envs/6e8595c2869a59bcfb9d3251dcc7705f
COPY workflow/envs/bwa-mem2.yaml /conda-envs/6e8595c2869a59bcfb9d3251dcc7705f/environment.yaml

# Conda environment:
#   source: workflow/envs/fastp.yaml
#   prefix: /conda-envs/f3b3ff3311107c3df70e1c86e328ec26
#   name: fastp
#   channels:
#     - bioconda
#   dependencies:
#     - fastp ==0.23.2
RUN mkdir -p /conda-envs/f3b3ff3311107c3df70e1c86e328ec26
COPY workflow/envs/fastp.yaml /conda-envs/f3b3ff3311107c3df70e1c86e328ec26/environment.yaml

# Conda environment:
#   source: workflow/envs/fastqc.yaml
#   prefix: /conda-envs/fa067ed38305bc1f79c0a5ec8cb9fa4e
#   name: fastqc
#   channels:
#     - bioconda
#   dependencies:
#     - fastqc ==0.11.9
RUN mkdir -p /conda-envs/fa067ed38305bc1f79c0a5ec8cb9fa4e
COPY workflow/envs/fastqc.yaml /conda-envs/fa067ed38305bc1f79c0a5ec8cb9fa4e/environment.yaml

# Conda environment:
#   source: workflow/envs/multiqc.yaml
#   prefix: /conda-envs/46a5165f902c6961e35090006ff1ab42
#   name: multiqc
#   channels:
#     - bioconda
#   dependencies:
#     - multiqc ==1.14
RUN mkdir -p /conda-envs/46a5165f902c6961e35090006ff1ab42
COPY workflow/envs/multiqc.yaml /conda-envs/46a5165f902c6961e35090006ff1ab42/environment.yaml

# Conda environment:
#   source: workflow/envs/qualimap.yaml
#   prefix: /conda-envs/9fba6491a3109390677fdbe34aa1e978
#   name: qualimap
#   channels:
#     - bioconda
#   dependencies:
#     - qualimap ==2.2.2a
RUN mkdir -p /conda-envs/9fba6491a3109390677fdbe34aa1e978
COPY workflow/envs/qualimap.yaml /conda-envs/9fba6491a3109390677fdbe34aa1e978/environment.yaml

# Conda environment:
#   source: workflow/envs/samtools.yaml
#   prefix: /conda-envs/c4c27e0d847cd651d475210d250344ad
#   name: samtools
#   channels:
#     - bioconda
#   dependencies:
#     - samtools ==1.16.1
RUN mkdir -p /conda-envs/c4c27e0d847cd651d475210d250344ad
COPY workflow/envs/samtools.yaml /conda-envs/c4c27e0d847cd651d475210d250344ad/environment.yaml

# Conda environment:
#   source: workflow/envs/star.yaml
#   prefix: /conda-envs/01443c9159c2e2c675f060065c7a970e
#   name: star
#   channels:
#     - bioconda
#   dependencies:
#     - star ==2.7.10b
#     - samtools ==1.16.1
RUN mkdir -p /conda-envs/01443c9159c2e2c675f060065c7a970e
COPY workflow/envs/star.yaml /conda-envs/01443c9159c2e2c675f060065c7a970e/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/50eb9a3e12a1bbfafcbb643868be3ca4 --file /conda-envs/50eb9a3e12a1bbfafcbb643868be3ca4/environment.yaml && \
    mamba env create --prefix /conda-envs/92c542c01e3b7bfb410fd03a54e00875 --file /conda-envs/92c542c01e3b7bfb410fd03a54e00875/environment.yaml && \
    mamba env create --prefix /conda-envs/6e8595c2869a59bcfb9d3251dcc7705f --file /conda-envs/6e8595c2869a59bcfb9d3251dcc7705f/environment.yaml && \
    mamba env create --prefix /conda-envs/f3b3ff3311107c3df70e1c86e328ec26 --file /conda-envs/f3b3ff3311107c3df70e1c86e328ec26/environment.yaml && \
    mamba env create --prefix /conda-envs/fa067ed38305bc1f79c0a5ec8cb9fa4e --file /conda-envs/fa067ed38305bc1f79c0a5ec8cb9fa4e/environment.yaml && \
    mamba env create --prefix /conda-envs/46a5165f902c6961e35090006ff1ab42 --file /conda-envs/46a5165f902c6961e35090006ff1ab42/environment.yaml && \
    mamba env create --prefix /conda-envs/9fba6491a3109390677fdbe34aa1e978 --file /conda-envs/9fba6491a3109390677fdbe34aa1e978/environment.yaml && \
    mamba env create --prefix /conda-envs/c4c27e0d847cd651d475210d250344ad --file /conda-envs/c4c27e0d847cd651d475210d250344ad/environment.yaml && \
    mamba env create --prefix /conda-envs/01443c9159c2e2c675f060065c7a970e --file /conda-envs/01443c9159c2e2c675f060065c7a970e/environment.yaml && \
    mamba clean --all -y
