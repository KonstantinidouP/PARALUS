# PARALUS - Installation Guide

This guide explains how to install PARALUS and its dependencies. The quickest way is the installer (`install_dependencies.sh`), which checks your system and installs everything that's needed. The manual steps for each system follow below. For a quick start, see the main [README.md](README.md).

## Contents

- [Quick install (recommended)](#quick-install-recommended)
- [Manual installation](#manual-installation): [Linux](#linux) | [Mac (Apple Silicon)](#mac-apple-silicon) | [HPC clusters](#hpc-clusters)

PARALUS needs the following tools: STAR (or HISAT2 as the low-RAM option), samtools, bedtools and cutadapt. python3 with matplotlib and numpy is optional and only needed for the PDF plots. PARALUS checks for the required tools at startup and stops with a list of any that are missing.

---

## Quick install (recommended)

### 1. Install PARALUS

Clone the repository to any location on your machine and add it to your `PATH`, so the `paralus` command works from any folder. The last three lines add two lines to `~/.bashrc`, a label and the `PATH` setting. Run them from inside the PARALUS folder, as shown below:

```bash
git clone https://github.com/KonstantinidouP/PARALUS.git
cd PARALUS
echo "# Added for PARALUS - delete this line and the next one to undo" >> ~/.bashrc
echo "export PATH=\"$PWD:\$PATH\"" >> ~/.bashrc
source ~/.bashrc
```

On a Mac, use `~/.zshrc` instead of `~/.bashrc`.

Check that it worked:

```bash
paralus --help
```

### 2. Install the dependencies

From inside the PARALUS folder, run:

```bash
bash install_dependencies.sh
```

The installer above works on Linux and Apple Silicon Macs. It checks which tools you already have and asks before installing anything. It then creates a conda environment called `paralus` with all PARALUS dependencies, installing Miniconda first if you don't have conda. On a Mac, it offers to reuse STAR if you've already compiled it, or compiles it for you. On computers with less than 32 GB RAM, it suggests HISAT2 instead of STAR.

Useful options:
- `--dry-run` shows what the installer would do, without installing anything
- `--hisat2` installs HISAT2 instead of STAR
- `--prefix /path/to/env` - by default, the environment is saved in your home folder. Use this option to save it elsewhere if your home folder has limited disk space (common on HPC clusters). Then activate it with `conda activate /path/to/env` instead of `conda activate paralus`

### 3. Before each run

```bash
conda activate paralus
```

---

## Manual installation

Follow the section for your system if you prefer to install the dependencies yourself instead of using the installer.

### Linux

If the required tools are already installed on your machine or server (e.g. a shared lab server), skip to [step 4](#4-install-paralus). To check, run `which STAR samtools bedtools cutadapt` (use `hisat2` instead of `STAR` for the low-RAM option): each tool that prints a path is available. If some are missing, a conda environment (steps 1-2) is the easiest way to add them. It installs the tools in your own account, needs no admin rights and keeps them separate from other software on the system.

#### 1. Install Miniconda

Skip this step if you already have conda.

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

#### 2. Create a conda environment with the dependencies

```bash
conda create -n paralus -c bioconda -c conda-forge \
    samtools bedtools cutadapt python matplotlib numpy
```

#### 3. Install STAR

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
conda install -n paralus -c bioconda -c conda-forge hisat2
hisat2 --version
```

HISAT2 needs ~8 GB RAM to build the hg38 index and ~5 GB to align. Run PARALUS with the `-H` flag, and pass `-x`/`-e` the HISAT2 index basename (e.g. `/path/to/genome_index/genome`, not the folder). See [REFERENCE_SETUP.md](REFERENCE_SETUP.md) for building the HISAT2 genome and rRNA indexes.

#### 4. Install PARALUS

Clone the repository to any location on your machine and add it to your `PATH`, so the `paralus` command works from any folder. The last three lines add two lines to `~/.bashrc`, a label and the `PATH` setting. Run them from inside the PARALUS folder, as shown below:

```bash
git clone https://github.com/KonstantinidouP/PARALUS.git
cd PARALUS
echo "# Added for PARALUS - delete this line and the next one to undo" >> ~/.bashrc
echo "export PATH=\"$PWD:\$PATH\"" >> ~/.bashrc
source ~/.bashrc
```

Check that it worked:

```bash
paralus --help
```

#### 5. Before each run

If you installed the dependencies in a conda environment (steps 1-2), activate it before each run:

```bash
conda activate paralus
```

---

### Mac (Apple Silicon)

These steps are for Macs with M1, M2, M3 or M4 chips. Most Macs don't come with these tools, so follow all steps unless you've installed them before.

#### 1. Install Miniconda

Skip this step if you already have conda.

```bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-arm64.sh
source ~/.zshrc
```

#### 2. Create a conda environment with the dependencies

```bash
conda create -n paralus -c bioconda -c conda-forge \
    samtools bedtools cutadapt python matplotlib numpy
```

> **Note:** The command above doesn't pin versions, so conda picks the latest compatible ones. You can pin a specific version if you need one (e.g. `samtools=1.17`), but if that version isn't available for ARM64, conda will stop with a `PackagesNotFoundError`.

#### 3. Install STAR (compile from source)

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
conda install -n paralus -c bioconda -c conda-forge hisat2
hisat2 --version
```

HISAT2 needs ~8 GB RAM to build the hg38 index and ~5 GB to align. Run PARALUS with the `-H` flag, and pass `-x`/`-e` the HISAT2 index basename (e.g. `/path/to/genome_index/genome`, not the folder). See [REFERENCE_SETUP.md](REFERENCE_SETUP.md) for building the HISAT2 genome and rRNA indexes.

#### 4. Install PARALUS

Clone the repository to any location on your Mac and add it to your `PATH`, so the `paralus` command works from any folder. The last three lines add two lines to `~/.zshrc`, a label and the `PATH` setting. Run them from inside the PARALUS folder, as shown below:

```bash
git clone https://github.com/KonstantinidouP/PARALUS.git
cd PARALUS
echo "# Added for PARALUS - delete this line and the next one to undo" >> ~/.zshrc
echo "export PATH=\"$PWD:\$PATH\"" >> ~/.zshrc
source ~/.zshrc
```

Check that it worked:

```bash
paralus --help
```

#### 5. Before each run

If you installed the dependencies in a conda environment (steps 1-2), activate it before each run:

```bash
conda activate paralus
```

---

### HPC clusters

On most clusters (e.g. Biowulf/NIH), the required tools are already installed centrally and only need to be loaded.

#### 1. Load the dependencies

Module names vary by cluster. Example for Biowulf/NIH:

```bash
module load star/2.7.11b
module load samtools
module load bedtools2
module load cutadapt
module load python
```

To use HISAT2 instead of STAR (`-H`), load `hisat2` in place of the `star` module. If some tools aren't available as modules on your cluster, use the [installer](#quick-install-recommended) instead. With `--prefix`, you can put the environment on project storage if your home folder has limited space.

#### 2. Install PARALUS

Clone the repository to your home or project directory and add it to your `PATH`, so the `paralus` command works from any folder. The last three lines add two lines to `~/.bashrc`, a label and the `PATH` setting. Run them from inside the PARALUS folder, as shown below:

```bash
git clone https://github.com/KonstantinidouP/PARALUS.git
cd PARALUS
echo "# Added for PARALUS - delete this line and the next one to undo" >> ~/.bashrc
echo "export PATH=\"$PWD:\$PATH\"" >> ~/.bashrc
source ~/.bashrc
```

Check that it worked:

```bash
paralus --help
```

#### 3. Before each run

Load the modules from step 1 in your session, or add the `module load` lines to the top of your job submission script.
