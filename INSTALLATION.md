# PARALUS - Installation Guide

This guide explains how to install PARALUS and its dependencies. Pick the section for your system and follow its steps in order. For a quick start, see the main [README.md](README.md).

## Contents

- [Linux](#linux)
- [Mac (Apple Silicon)](#mac-apple-silicon)
- [HPC clusters](#hpc-clusters)

PARALUS needs the following tools: STAR (or HISAT2 as the low-RAM option), samtools, bedtools and cutadapt. python3 with matplotlib and numpy is optional and only needed for the PDF plots. PARALUS checks for the required tools at startup and stops with a list of any that are missing.

---

## Linux

If the required tools are already installed on your machine or server (e.g. a shared lab server), skip to [step 4](#4-install-paralus). To check, run `which STAR samtools bedtools cutadapt` (use `hisat2` instead of `STAR` for the low-RAM option): each tool that prints a path is available. If some are missing, a conda environment (steps 1-2) is the easiest way to add them. It installs the tools in your own account, needs no admin rights and keeps them separate from other software on the system.

### 1. Install Miniconda

Skip this step if you already have conda.

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

### 2. Create a conda environment with the dependencies

```bash
conda create -n rnaseq -c bioconda -c conda-forge \
    samtools bedtools cutadapt python matplotlib numpy
```

### 3. Install STAR

We recommend the official pre-compiled static binary rather than conda, as it avoids version conflicts. The commands below install STAR 2.7.11b, the version PARALUS was tested with. To use a newer release, replace `2.7.11b` with its version number (see https://github.com/alexdobin/STAR/releases). If you change STAR versions, you may need to rebuild your genome index.

Download and unzip STAR:

```bash
wget https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR_2.7.11b.zip
unzip STAR_2.7.11b.zip
```

Then install it in one of two ways:

**With admin (`sudo`) rights**, install STAR system-wide:

```bash
sudo cp STAR_2.7.11b/Linux_x86_64_static/STAR /usr/local/bin/
STAR --version
```

**Without admin rights**, as is common on shared servers, install STAR in your home folder:

```bash
mkdir -p ~/bin
cp STAR_2.7.11b/Linux_x86_64_static/STAR ~/bin/
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
STAR --version
```

**Low-RAM option with HISAT2 (for local machines with less than 32 GB RAM).** Install HISAT2 into the same environment instead of STAR:

```bash
conda install -n rnaseq -c bioconda -c conda-forge hisat2
hisat2 --version
```

HISAT2 needs ~8 GB RAM to build the hg38 index and ~5 GB to align. Run PARALUS with the `-H` flag, and pass `-x`/`-e` the HISAT2 index basename (e.g. `/path/to/genome_index/genome`, not the folder). See [README.md](README.md) for building the HISAT2 genome and rRNA indexes.

### 4. Install PARALUS

Clone the repository to any location on your machine. PARALUS is a bash script and does not need to be installed in the conda environment. The version number may be newer than the one shown here (`paralus_1.8.7.sh`). Use the file name you find in the repository:

```bash
git clone https://github.com/KonstantinidouP/PARALUS.git
cd PARALUS
chmod +x paralus_1.8.7.sh
./paralus_1.8.7.sh --help
```

### 5. Before each run

If you installed the dependencies in a conda environment (steps 1-2), activate it before each run:

```bash
conda activate rnaseq
```

---

## Mac (Apple Silicon)

These steps are for Macs with M1, M2, M3 or M4 chips. Most Macs don't come with these tools, so follow all steps unless you've installed them before.

### 1. Install Miniconda

Skip this step if you already have conda.

```bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-arm64.sh
source ~/.zshrc
```

### 2. Create a conda environment with the dependencies

```bash
conda create -n rnaseq -c bioconda -c conda-forge \
    samtools bedtools cutadapt python matplotlib numpy
```

> **Note:** The command above doesn't pin versions, so conda picks the latest compatible ones. You can pin a specific version if you need one (e.g. `samtools=1.17`), but if that version isn't available for ARM64, conda will stop with a `PackagesNotFoundError`.

### 3. Install STAR (compile from source)

STAR has no pre-compiled binary for Apple Silicon, so it needs to be compiled. The commands below compile STAR 2.7.11b, the version PARALUS was tested with. To use a newer release, replace `2.7.11b` with its version number (see https://github.com/alexdobin/STAR/releases). If you change STAR versions, you may need to rebuild your genome index.

**a. Install Homebrew and GNU gcc.**

Install Homebrew from https://brew.sh (required to install gcc). Skip this if `brew --version` already works.

Then install GNU gcc:

```bash
brew install gcc
ls /opt/homebrew/bin/g++*   # note the version in the output, e.g. g++-15
```

**b. Clone and compile STAR:**

```bash
cd ~
git clone https://github.com/alexdobin/STAR.git
cd STAR
git checkout 2.7.11b
cd source
make STARforMacStatic CXX=/opt/homebrew/bin/g++-15 CXXFLAGS_SIMD=
```

Replace `g++-15` with the version you found above. `CXXFLAGS_SIMD=` removes the `-mavx2` flag, which is x86-only and does not exist on ARM64. Compilation takes a few minutes, and the warnings it prints are harmless.

**c. Install and verify:**

```bash
sudo cp ~/STAR/source/STAR /usr/local/bin/
STAR --version
```

**Low-RAM option with HISAT2 (for local machines with less than 32 GB RAM).** Install HISAT2 into the same environment instead of STAR. It installs through conda, so no compiling is needed:

```bash
conda install -n rnaseq -c bioconda -c conda-forge hisat2
hisat2 --version
```

HISAT2 needs ~8 GB RAM to build the hg38 index and ~5 GB to align. Run PARALUS with the `-H` flag, and pass `-x`/`-e` the HISAT2 index basename (e.g. `/path/to/genome_index/genome`, not the folder). See [README.md](README.md) for building the HISAT2 genome and rRNA indexes.

### 4. Install PARALUS

Clone the repository to any location on your Mac. PARALUS is a bash script and does not need to be installed in the conda environment. The version number may be newer than the one shown here (`paralus_1.8.7.sh`). Use the file name you find in the repository:

```bash
git clone https://github.com/KonstantinidouP/PARALUS.git
cd PARALUS
chmod +x paralus_1.8.7.sh
./paralus_1.8.7.sh --help
```

### 5. Before each run

If you installed the dependencies in a conda environment (steps 1-2), activate it before each run:

```bash
conda activate rnaseq
```

---

## HPC clusters

On most clusters (e.g. Biowulf/NIH), the required tools are already installed centrally and only need to be loaded.

### 1. Load the dependencies

Module names vary by cluster. Example for Biowulf/NIH:

```bash
module load star/2.7.11b
module load samtools
module load bedtools2
module load cutadapt
module load python
```

To use HISAT2 instead of STAR (`-H`), load `hisat2` in place of the `star` module. If a tool isn't available as a module on your cluster, install it with conda as described in the [Linux](#linux) section.

### 2. Install PARALUS

Clone the repository to your home or project directory. The version number may be newer than the one shown here (`paralus_1.8.7.sh`). Use the file name you find in the repository:

```bash
git clone https://github.com/KonstantinidouP/PARALUS.git
cd PARALUS
chmod +x paralus_1.8.7.sh
./paralus_1.8.7.sh --help
```

### 3. Before each run

Load the modules from step 1 in your session, or add the `module load` lines to the top of your job submission script.
