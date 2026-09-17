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

## Installing HISAT2 (alternative to STAR - for Macs with 16 GB RAM or less)

STAR needs its whole genome index in RAM, both to build it and to align against it - roughly 30 GB for human/mouse. That's out of reach on most Mac laptops with 16 GB of RAM or less. HISAT2 is a drop-in alternative aligner that uses an FM-index instead, needing only ~4-5 GB for the same genomes, both to build and to align against.

If your machine can't comfortably run STAR, install HISAT2 instead and pass PARALUS the `-H` flag - see `./paralus_1.8.5.sh --help` for what changes (`-x`/`-e` then take a HISAT2 index basename, built with `hisat2-build`, instead of a STAR genomeDir).

HISAT2 installs cleanly via conda on both Linux and Apple Silicon, so no separate compile-from-source path is needed here the way it is for STAR:

```bash
conda install -n rnaseq -c bioconda -c conda-forge hisat2
```

(or add `hisat2` to the package list in the `conda create` command below if you're setting the environment up for the first time.)

Verify:

```bash
hisat2 --version
```

Building a HISAT2 index is a one-time step per genome, similar to STAR's `--genomeGenerate`, but with far lower memory and disk requirements:

```bash
mkdir -p /path/to/genome_index
hisat2-build /path/to/hg38.fa /path/to/genome_index/genome
```

Pass the resulting basename (`/path/to/genome_index/genome`, not a directory) to PARALUS with `-x -H`. The same applies to the optional rRNA index used with `-e` - see the HISAT2 note in [README.md](README.md#step-3---download-the-reference-genome-and-build-a-star-index-run-once-per-genome) for the equivalent rRNA-index command.

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

> **Tip:** if you're using HISAT2 instead of STAR (see above, recommended for Macs with ≤16 GB RAM), add `hisat2` to this package list so it's installed in the same environment: `conda create -n rnaseq -c bioconda -c conda-forge samtools bedtools cutadapt python matplotlib numpy hisat2`.

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

Run these lines before running PARALUS or add them to the top of your submission script. If using HISAT2 instead of STAR (`-H`), load `module load hisat2` in place of the `star` module - exact module name varies by cluster.

---

## Installing PARALUS

Clone the repository to any location on your machine - PARALUS is a bash script and does not need to be installed in the conda environment:

```bash
git clone https://github.com/KonstantinidouP/PARALUS.git
cd PARALUS
chmod +x paralus_1.8.5.sh
```

To view all available options:

```bash
./paralus_1.8.5.sh --help
```

Before running PARALUS, make sure all dependencies are available in your environment. If using conda:

```bash
conda activate rnaseq
```

On HPC clusters, load the required modules as described above instead.
