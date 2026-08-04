# PARALUS — Detailed Installation Guide

This document provides detailed installation instructions for different environments.

---

PARALUS requires the following tools:

> **Optional — PDF plots:** PARALUS can generate k-mer enrichment and metagene plots as PDF files if `python3` with `matplotlib` and `numpy` are available in your environment. If they are not found, the pipeline skips the plots gracefully and notes it in the `_stats.txt` file — the core analysis is unaffected. To enable plots, install them in your conda environment:
> ```bash
> conda install -c conda-forge matplotlib numpy
> ```

- [STAR](https://github.com/alexdobin/STAR) (tested with v2.7.2b)
- [samtools](http://www.htslib.org/) (tested with v1.17)
- [bedtools](https://bedtools.readthedocs.io/) (tested with v2.29.2)
- [cutadapt](https://cutadapt.readthedocs.io/)
- python
- perl (standard installation, used for k-mer generation)

Follow the instructions below for your environment.

---

### A — Apple Silicon Mac (M1/M2/M3/M4)

If you don't have conda, install Miniconda. Download the **arm64** installer:

```bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
shasum -a 256 Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-arm64.sh
source ~/.zshrc
```

Create a conda environment with all dependencies **except STAR**:

```bash
conda create -n rnaseq -c bioconda -c conda-forge     samtools bedtools cutadapt python perl
conda activate rnaseq
```

> **Note:** Specific version pins (e.g. `samtools=1.17`) may not be available for ARM64. Drop version pins and let conda pick the latest compatible versions if you get a `PackagesNotFoundError`.

**STAR requires special handling on Apple Silicon — the conda build does not work.**

The bioconda ARM64 build of STAR has a confirmed bug where it silently reads 0 input reads. The root cause is that the binary uses `vfork()` to spawn child processes for reading input, and macOS restricts this in ad-hoc signed binaries. STAR must be compiled from source:

**1. Install Homebrew and GNU gcc**

Install Homebrew from https://brew.sh if you don't have it, then:

```bash
brew install gcc
ls /opt/homebrew/bin/g++*   # note the version, e.g. g++-15
```

**2. Clone and compile STAR**

```bash
cd ~
git clone https://github.com/alexdobin/STAR.git
cd STAR/source
make STARforMacStatic CXX=/opt/homebrew/bin/g++-15 CXXFLAGS_SIMD=
```

Replace `g++-15` with the version you found above. `CXXFLAGS_SIMD=` suppresses the `-mavx2` flag, which is x86-only and does not exist on ARM64 — this does not affect performance.

**3. Install into conda environment**

```bash
cp ~/STAR/source/STAR ~/miniconda3/envs/rnaseq/bin/STAR
```

**4. Verify**

```bash
conda activate rnaseq
STAR --version
```

---

### B — Linux local machine (Intel/AMD x86)

If you don't have conda, first install Miniconda. Verify the integrity of the installer by comparing the SHA256 hash shown in your terminal with the one listed at https://docs.anaconda.com/miniconda/miniconda-hashes/:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
shasum -a 256 Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

Follow the prompts: press Space to scroll through the license, type `yes` to accept, press Enter to confirm the default install location, and type `yes` to initialize conda.

Create the environment with all dependencies including STAR:

```bash
conda create -n rnaseq -c bioconda -c conda-forge     star=2.7.2b     samtools=1.17     bedtools=2.29.2     cutadapt     python     perl
conda activate rnaseq
```

Activate before each run:

```bash
conda activate rnaseq
```

---

### C — Shared server (e.g. pinbot), install with conda

If conda is not already available, install Miniconda in your home directory:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

Then create and activate the environment:

```bash
conda create -n rnaseq -c bioconda -c conda-forge     star=2.7.2b     samtools=1.17     bedtools=2.29.2     cutadapt     python     perl
conda activate rnaseq
```

Activate before each run:

```bash
conda activate rnaseq
```

---

### D — HPC cluster (e.g. Biowulf/NIH), use module load

On clusters that use the `module` system, load the required tools before running PARALUS:

```bash
module load star/2.7.2b
module load samtools/1.17
module load bedtools2/2.29.2
module load cutadapt
module load python
```

Add these `module load` lines to the top of your submission script (see `RunParalus_for_hg38.sh` for an example). The exact module names may vary by cluster — check available versions with `module avail star`.