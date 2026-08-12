# PARALUS - Dependencies Installation Guide

This document provides detailed instructions for the installation of dependencies on Linux local machines and Apple Silicon Macs. If you are working on an HPC cluster (e.g. Biowulf), the required tools are likely already installed centrally and only need to be loaded before running PARALUS - refer to your cluster's documentation for instructions (e.g. `module load STAR`). For a quick start, see the main [README.md](README.md).

---

## Installing STAR

STAR is the most environment-sensitive dependency. We recommend installing it from the official pre-compiled static binary rather than via conda, as this avoids version conflicts and has been tested to work reliably.

### Linux

Download and install the pre-compiled static binary from the STAR GitHub releases page:

```bash
wget https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR_2.7.11b.zip
unzip STAR_2.7.11b.zip
sudo cp STAR_2.7.11b/Linux_x86_64_static/STAR /usr/local/bin/
```

Verify:

```bash
STAR --version
```

### Apple Silicon Mac (M1/M2/M3/M4) - compile from source



1. Install Homebrew and GNU gcc

Install Homebrew from https://brew.sh if you don't have it, then:

```bash
brew install gcc
ls /opt/homebrew/bin/g++*   # note the version, e.g. g++-15
```

2. Clone and compile STAR

```bash
cd ~
git clone https://github.com/alexdobin/STAR.git
cd STAR
git checkout 2.7.11b
cd source
make STARforMacStatic CXX=/opt/homebrew/bin/g++-15 CXXFLAGS_SIMD=
```

Replace `g++-15` with the version you found above. `CXXFLAGS_SIMD=` suppresses the `-mavx2` flag, which is x86-only and does not exist on ARM64 - this does not affect performance. Compilation takes a few minutes; warnings are harmless.

3. Install and verify:

```bash
sudo cp ~/STAR/source/STAR /usr/local/bin/
STAR --version
```

---

## Installing all other dependencies via conda

### Installing Miniconda

If you don't have conda, install Miniconda first.

**Linux:**

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
shasum -a 256 Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

**Mac (Apple Silicon):**

```bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
shasum -a 256 Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-arm64.sh
source ~/.zshrc
```

Verify the installer hash against the one listed at https://docs.anaconda.com/miniconda/miniconda-hashes/ before running it.

### Create a conda environment and install dependencies

```bash
conda create -n rnaseq -c bioconda -c conda-forge \
    samtools bedtools cutadapt python matplotlib numpy
```

> **Note:** Specific version pins (e.g. `samtools=1.17`) may not be available for ARM64. Drop version pins and let conda pick the latest compatible versions if you get a `PackagesNotFoundError`. Including both `-c bioconda` and `-c conda-forge` is important - some packages may fail to resolve dependencies with bioconda alone on Apple Silicon.

### Activate the conda environment before running PARALUS

```bash
conda activate rnaseq
```

---

## HPC clusters (e.g. Biowulf/NIH) - module load

On clusters where the required tools are already centrally installed, load them before running PARALUS. The exact module names vary by cluster - the following is an example for Biowulf/NIH:

```bash
module load star/2.7.11b
module load samtools
module load bedtools2
module load cutadapt
module load python
```

Run these lines before running PARALUS or add them to the top of your submission script.

---

## Installing PARALUS

Clone the repository to any location on your machine - PARALUS is a bash script and does not need to be installed in the conda environment:

```bash
git clone https://github.com/KonstantinidouP/PARALUS.git
cd PARALUS
chmod +x paralus_1.7.1.sh
```

To view all available options:

```bash
./paralus_1.7.1.sh --help
```

Before running PARALUS, make sure all dependencies are available in your environment. If using conda:

```bash
conda activate rnaseq
```

On HPC clusters, load the required modules as described above instead.
