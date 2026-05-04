# PARALUS

> **[PARALUS stands for: placeholder — fill in when confirmed]**

PARALUS is a bash-based pipeline for the analysis of PAR-CLIP and CLIP sequencing data. It performs end-to-end processing from raw reads to annotated crosslink clusters, and includes k-mer enrichment analysis to identify sequence motifs at RNA-protein binding sites.

---

## Overview

PARALUS takes CLIP or PAR-CLIP sequencing data as input and outputs high-confidence RNA-protein binding clusters with full annotation, crosslink site information, and binding statistics. The pipeline handles adapter trimming, read collapsing, genome alignment, T-to-C crosslink site extraction, peak calling with read distribution modeling, and genomic annotation — all in a single run.

A companion script, **masterchef**, builds the custom annotation BED file required by PARALUS from any Gencode/Ensembl GTF file.

---

## Features

- Handles multiple input formats: FASTQ, FASTQ.gz, and FASTA
- Optional UMI trimming and adapter removal via cutadapt
- Optional rRNA depletion by pre-mapping to an rRNA genome index
- T-to-C crosslink site extraction from PAR-CLIP data
- Two-stage peak calling based on read distribution modeling (normal and non-normal clusters)
- Genomic annotation with priority-ranked feature assignment (3ss > 5ss > start codon > stop codon > 3'UTR > 5'UTR > CDS > exon > intron)
- 5-mer enrichment analysis at binding sites vs. background
- Per-sample statistics report

---

## Performance & Accessibility

PARALUS is designed to be fast, self-contained, and accessible — running on both Linux servers and local Mac machines with a single command and no workflow manager, containers, or cluster required.

On an Apple MacBook Pro (M4 Pro, 48 GB RAM) using 4 threads, a full PAR-CLIP run on a 512 MB FASTQ file (~8.8 million reads) completes in approximately **15 minutes**, performing the following steps end to end:

- Read collapsing, UMI trimming, and adapter removal
- Optional rRNA depletion by pre-mapping to an rRNA genome
- Genome alignment against hg38
- T-to-C crosslink site extraction at single-nucleotide resolution
- Two-stage peak calling with read distribution modeling
- Full genomic annotation with priority-ranked feature assignment
- Sequence extraction at binding sites
- 5-mer enrichment analysis at binding sites vs. background

Many comparable CLIP analysis tools require multiple separate steps, different software environments, and are often limited to Linux clusters. PARALUS performs all of these steps in a single command, making it practical for everyday use on a local machine.

Runtime is consistent and reproducible across runs. Alignment and sorting steps benefit from additional threads (set with `-p`), but some steps such as peak calling and 5-mer counting are inherently serial. On local machines, disk I/O speed and available RAM can be as important as CPU core count — a fast SSD and sufficient RAM to cache large files like the genome FASTA in memory will improve performance significantly. For large datasets or batch processing of multiple samples, running on a Linux server is recommended.

## Dependencies

PARALUS requires the following tools:

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
## Installation

No installation is required. Clone the repository and run the scripts directly:

```bash
git clone https://github.com/KonstantinidouP/PARALUS.git
cd PARALUS
```

---

## Quick Start

### Step 1 — Build the annotation file (run once per genome)

masterchef requires two input files: a GTF annotation file and a tRNA BED file. Both are freely available from public databases.

> **Quick option:** If you want to skip the manual downloads, a helper script is provided that downloads all required reference files automatically:
>
> ```bash
> bash download_references.sh
> ```
>
> If any link has become outdated, the script will tell you where to go to download the file manually. You can also skip the script and follow the steps below.

#### 1a — Download the GTF file

Download the comprehensive gene annotation GTF from **GENCODE** (https://www.gencodegenes.org). GENCODE is strongly recommended over other sources (e.g. Ensembl, RefSeq) because masterchef relies on the `gene_type` biotype field, which is formatted consistently in GENCODE GTFs.

Always use the latest release for your assembly. For human (hg38), as of 2025 this is release 49:

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
gunzip gencode.v49.annotation.gtf.gz
```

For mouse (mm10), check the current release at https://www.gencodegenes.org/mouse/ and download accordingly.

> **Note:** Check the GENCODE release history page for the latest version before downloading, as new releases are issued regularly.

#### 1b — Download the tRNA BED file

Download the tRNA annotation BED file from the **Genomic tRNA Database (GtRNAdb)** (https://gtrnadb.ucsc.edu), which provides tRNA gene predictions generated by tRNAscan-SE for a wide range of genomes.

For human (hg38):

1. Navigate to https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/
2. Click **Download tRNAscan-SE Results** — this downloads a `.tar.gz` archive containing the BED file along with other result files
3. Extract it:

```bash
tar -xzf hg38-tRNAs.tar.gz
```

The BED file (`hg38-tRNAs.bed`) will be inside the extracted folder.

#### 1c — Run masterchef

```bash
bash masterchef_v2.sh -g gencode.v49.annotation.gtf -t hg38-tRNAs.bed
```

This produces a `*_Paralus.bed` annotation file to use with PARALUS.

### Step 2 — Download the reference genome and build a STAR index (run once per genome)

You need the reference genome of your organism of interest in FASTA format. If you don't have it yet, download it from UCSC. Note that these are large files (~1 GB compressed) — make sure you have enough disk space before downloading.

For human (hg38):

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
```

For mouse (mm10):

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/latest/mm10.fa.gz
```

Once you have the FASTA file, decompress it first — STAR requires an uncompressed FASTA for genome index generation:

```bash
gunzip -k hg38.fa.gz
```

The `-k` flag keeps the original compressed file. The decompressed file will be ~3.5 GB for hg38.

Then build the STAR genome index. This is a one-time step per genome and requires ~30 GB of disk space and significant RAM (≥32 GB recommended). On a Mac with 48 GB RAM (Apple M4 Pro) and 6 cores, indexing hg38 took ~44 minutes:

```bash
mkdir /path/to/genome_index
STAR --runMode genomeGenerate \
     --genomeDir /path/to/genome_index \
     --genomeFastaFiles /path/to/hg38.fa \
     --runThreadN 6 \
     --limitGenomeGenerateRAM 38000000000
```

> **Note:** The `--limitGenomeGenerateRAM` flag caps memory usage in bytes. 38 GB works well on a 48 GB machine, leaving enough headroom for the OS. Adjust this value based on your available RAM. On a server with ≥64 GB RAM you can omit this flag entirely.

Pass the genome index folder to PARALUS with the `-x` flag, and the uncompressed FASTA file with the `-f` flag.

### Step 2b — Download rRNA sequences and build an rRNA index (optional)

This step is only needed if you want to use the rRNA depletion option in PARALUS (the `-e` flag). It removes rRNA-mapping reads before genome alignment, which can significantly improve the signal-to-noise ratio in your data. For this you will need to build and index an rRNA genome using STAR.

Below are instructions for where to find the human rRNA sequences. Both files are available from the **NCBI Nucleotide database** (https://www.ncbi.nlm.nih.gov/nucleotide):

**1. The 45S rDNA repeat unit** — contains the 18S, 5.8S, and 28S rRNA sequences:
- Go to https://www.ncbi.nlm.nih.gov/nuccore/U13369.1
- Click **Send to** (top right) → **File** → Format: **FASTA** → **Create File**

**2. The 5S rRNA** — encoded separately in the human genome and not part of the 45S repeat:
- Go to https://www.ncbi.nlm.nih.gov/nuccore/X12811.1
- Click **Send to** (top right) → **File** → Format: **FASTA** → **Create File**

> **Note:** The 5S rRNA is located on a separate chromosomal locus from the 45S rDNA cluster, which is why they are downloaded separately.

Once both files are downloaded, combine them and build a STAR index:

```bash
# Combine into a single rRNA FASTA
cat U13369.1.fasta X12811.1.fasta > rRNA.fa

# Build STAR index
mkdir /path/to/rRNA_index
STAR --runMode genomeGenerate      --genomeDir /path/to/rRNA_index      --genomeFastaFiles rRNA.fa      --runThreadN 6      --genomeSAindexNbases 6
```

> **Note:** The `--genomeSAindexNbases` flag is required for small genomes. According to the STAR manual, the correct value should be scaled down using `min(14, log2(GenomeLength)/2 - 1)`. For the combined human rRNA sequences (~43 kb total), this gives a value of 6. If you use different rRNA sequences, adjust accordingly. Because the rRNA genome is so small, indexing completes almost instantly regardless of the number of threads used.

Pass the rRNA index folder to PARALUS with the `-e` flag.

### Step 3 — Run PARALUS

```bash
bash paralus1.5.9.sh \
  -i sample.fastq.gz \
  -x /path/to/genome_index \
  -A /path/to/genes_Paralus.bed \
  -f /path/to/genome.fa \
  -a NNTGACTGTGGAATTCTCGGGTGCCAAGG \
  -u 4 \
  -p 12
```

An example submission script for hg38 is provided in `RunParalus_for_hg38.sh`.

---

## Parameters

### Required

| Flag | Description |
|------|-------------|
| `-i` | Input file (FASTQ.gz / FASTQ / FASTA) |
| `-x` | Path to STAR genome index |
| `-A` | Annotation BED file (generated by masterchef) |
| `-f` | Genome FASTA file |

### Optional — Preprocessing

| Flag | Description | Default |
|------|-------------|---------|
| `-a` | 3' adapter sequence | none |
| `-d` | 5' adapter sequence | `GTTCAGAGTTCTACAGTCCGACGATC` |
| `-u` | Bases to trim from 5' end (UMIs) | `0` |
| `-n` | Minimum read length after trimming | `19` |
| `-p` | Threads | `8` |
| `-e` | rRNA genome index for pre-depletion | none |

### Optional — Cluster calling

| Flag | Description | Default |
|------|-------------|---------|
| `-m` | Minimum reads to define a group | `3` |
| `-c` | Minimum T-to-C conversion locations per cluster | `1` |
| `-r` | Minimum fraction of T-to-C reads per cluster | `0.3` |

---

## Output files

| File | Description |
|------|-------------|
| `*.clusters.csv` | Main output — annotated binding clusters with crosslink statistics, gene name, annotation, biotype, and binding sequence |
| `*.clusters.bed` | BED format version of the clusters file |
| `*.aligned_TtoC.bam` | BAM file containing only T-to-C crosslink-supporting reads |
| `*.aligned_TtoC.bam.bai` | Index for the TtoC BAM file |
| `*.aligned_TtoC_loc.bed` | BED file of individual T-to-C crosslink site positions |
| `*.Aligned.sortedByCoord.out.bam` | Full genome-aligned BAM file |
| `*_stats.txt` | Per-sample run statistics (read counts, cluster counts by feature type) |
| `*_5mers.txt` | 5-mer counts in bound vs. background sequences |
| `*_tri_3mers.txt` | Tri-3-mer counts in bound vs. background sequences |

---

## Repository structure

```
PARALUS/
├── paralus1.5.9.sh          # Main pipeline script
├── masterchef_v2.sh         # Annotation BED file builder
├── RunParalus_for_hg38.sh   # Example run script for hg38
└── README.md
```

---

## Citation

> **[Placeholder — manuscript in preparation]**
> Anastasakis D, Konstantinidou P, et al. PARALUS: a pipeline for PAR-CLIP and CLIP data analysis. *Journal TBD*, year.

If you use PARALUS in your work, please cite the above once published.

---

## Authors

- **Dimitrios Anastasakis** — lead developer
- **Parthena Konstantinidou** — co-developer

---

## License

> **[Placeholder — license to be decided]**
> 
> For now, this repository is private. License will be added prior to public release.