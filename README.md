# PARALUS

> **A Lightweight Unix Shell Pipeline for PAR-CLIP and CLIP Sequencing Data Analysis**

PARALUS is a bash-based pipeline for the analysis of PAR-CLIP and CLIP sequencing data. It performs end-to-end processing from raw reads to annotated crosslink clusters, and includes k-mer enrichment analysis to identify sequence preferences at RNA-protein binding sites.

---

## Overview

PARALUS takes CLIP or PAR-CLIP sequencing data as input and outputs high-confidence RNA-protein binding clusters with full annotation, crosslink site information, and binding statistics. The pipeline handles adapter trimming, read collapsing, genome alignment, T-to-C crosslink site extraction, peak calling with read distribution modeling, and genomic annotation — all in a single run.

PARALUS includes **Trierarch**, a built-in annotation builder that generates the custom annotation BED file directly from a Gencode/Ensembl GTF file. Trierarch is embedded in the same script and can be run standalone using the `--trierarch` mode, or invoked automatically during a PARALUS run by providing a GTF file with the `-G` flag instead of a pre-built annotation BED file.

---

## Features

- Handles multiple input formats: FASTQ, FASTQ.gz, and FASTA
- Optional UMI trimming and adapter removal via cutadapt
- Optional rRNA depletion by pre-mapping to an rRNA genome index
- T-to-C crosslink site extraction from PAR-CLIP data at single-nucleotide resolution
- Two-stage peak calling based on read distribution modeling (normal and non-normal clusters)
- Gene name, feature, and biotype information per cluster
- Resolution of overlapping annotations by priority-ranked feature assignment
- Built-in annotation builder (Trierarch) — provide a GTF file and the annotation BED is built automatically; supports canonical transcripts only (`-C`), optional tRNA (`-T`) and miRNA (`-M`) databases
- 5-mer enrichment analysis at binding sites vs. background, with optional PDF scatter plot output
- Metagene analysis plots for UTR/CDS and intron binding distribution (optional PDF output)
- Full run logging to a timestamped log file
- Per-sample statistics report including binding breakdown by gene feature type

---

## Performance & Accessibility

PARALUS is designed to be fast, self-contained, and accessible — running on both Linux servers and local Mac machines with a single command.

On an Apple MacBook Pro (M4 Pro, 48 GB RAM) using 4 threads, PARALUS was benchmarked on two datasets of different sizes, both completing on a local machine:

| Dataset size | Approx. reads | Runtime | CPU utilization |
|---|---|---|---|
| 512 MB FASTQ.gz | ~8.8 million | ~7 minutes | 171% |
| 4.7 GB FASTQ.gz | ~10x larger | ~100 minutes | 111% |

Runtime is consistent and reproducible across runs. Alignment and sorting steps benefit from additional threads (set with `-p`), but some steps such as peak calling and 5-mer counting are inherently serial. On local machines, disk I/O speed and available RAM can be as important as CPU core count — a fast SSD and sufficient RAM to cache large files like the genome FASTA in memory will improve performance significantly. For large datasets or batch processing of multiple samples, running on a Linux server is recommended.

---

## Dependencies

PARALUS requires the following tools installed in your environment:

- [STAR](https://github.com/alexdobin/STAR) (tested with v2.7.2b)
- [samtools](http://www.htslib.org/) (tested with v1.17)
- [bedtools](https://bedtools.readthedocs.io/) (tested with v2.29.2)
- [cutadapt](https://cutadapt.readthedocs.io/)
- python3
- matplotlib (OPTIONAL — required for PDF plot output)
- numpy (OPTIONAL — required for PDF plot output)

For Mac, Linux local machines and shared Linux servers, all dependencies can be installed via conda:

```bash
conda create -n rnaseq -c bioconda -c conda-forge \
    star samtools bedtools cutadapt python matplotlib numpy
conda activate rnaseq
```

> **Apple Silicon Mac (M1/M2/M3/M4):** The conda build of STAR does not work on ARM64 due to a process spawning bug. STAR must be compiled from source. See [INSTALLATION.md](INSTALLATION.md) for detailed instructions.

For detailed installation instructions across different environments (Apple Silicon Mac, Linux, shared servers, HPC clusters), see [INSTALLATION.md](INSTALLATION.md).

---

## Quick Start

### Step 1 — Clone the repository and activate the environment

Clone the repository to any location on your machine — PARALUS is a bash script and does not need to be installed in the conda environment:

```bash
git clone https://github.com/KonstantinidouP/PARALUS.git
cd PARALUS
```

Before running PARALUS, make sure all dependencies are available in your environment. If using conda:

```bash
conda activate rnaseq
```

On HPC clusters (e.g. Biowulf), load each dependency using your cluster's module system instead (see [INSTALLATION.md](INSTALLATION.md)).

### Step 2 — Build the annotation file (run once per genome)

Trierarch requires a gene annotation GTF file and optionally a tRNA BED file and/or a miRNA GFF3 file. The GTF can be from GENCODE or Ensembl. All reference files are freely available from public databases.

We recommend downloading the latest comprehensive gene annotation GTF for your organism from **GENCODE** (https://www.gencodegenes.org) or **Ensembl** (https://www.ensembl.org), and optionally the tRNA annotation BED file from the **Genomic tRNA Database (GtRNAdb)** (https://gtrnadb.ucsc.edu) and/or the miRNA annotation GFF3 from **miRBase** (https://www.mirbase.org/download/). Download links and instructions for various organisms can be found in the PARALUS help page. Run the following command to view:

```bash
bash paralus_1.7.1.sh --trierarch --help
```

#### Build the annotation BED file with Trierarch

Trierarch is built into PARALUS and produces the annotation BED file from the GTF. Using canonical transcripts only (the `-C` option) is recommended for cleaner annotation:

```bash
bash paralus_1.7.1.sh --trierarch \
  -G gencode.v49.annotation.gtf \
  -T hg38-tRNAs.bed \
  -C \
  -O gencode.v49.annotation_Trierarch_canonical
```

This produces `gencode.v49.annotation_Trierarch_canonical_Paralus.bed` to use with PARALUS via the `-A` flag. Alternatively, use `-G` directly in your PARALUS run and Trierarch will build the annotation automatically.

### Step 3 — Download the reference genome and build a STAR index (run once per genome)

You need the reference genome of your organism of interest in FASTA format. If you don't have it yet, download it from UCSC. These are large files (~1 GB compressed) — make sure you have enough disk space.

For human (hg38):

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
```

Decompress the FASTA — STAR requires an uncompressed file for genome index generation:

```bash
gunzip -k hg38.fa.gz
```

The decompressed file will be ~3.5 GB for hg38.

Then build the STAR genome index. This is a one-time step per genome and requires ~30 GB of disk space and significant RAM (≥32 GB recommended):

```bash
mkdir /path/to/genome_index
STAR --runMode genomeGenerate \
     --genomeDir /path/to/genome_index \
     --genomeFastaFiles /path/to/hg38.fa \
     --runThreadN 6 \
     --limitGenomeGenerateRAM 38000000000
```

> **Note:** The `--limitGenomeGenerateRAM` flag caps memory usage in bytes. Adjust based on your available RAM — on a server with ≥64 GB RAM you can omit this flag entirely.

Pass the genome index folder to PARALUS with the `-x` flag, and the uncompressed FASTA file with the `-f` flag.

### Step 3b — Download rRNA sequences and build an rRNA index (OPTIONAL — run once per genome)

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
STAR --runMode genomeGenerate \
     --genomeDir /path/to/rRNA_index \
     --genomeFastaFiles rRNA.fa \
     --runThreadN 6 \
     --genomeSAindexNbases 6
```

> **Note:** The `--genomeSAindexNbases` flag is required for small genomes. According to the STAR manual, the correct value should be scaled down using `min(14, log2(GenomeLength)/2 - 1)`. For the combined human rRNA sequences (~43 kb total), this gives a value of 6. If you use different rRNA sequences, adjust accordingly. Because the rRNA genome is so small, indexing completes almost instantly regardless of the number of threads used.

Pass the rRNA index folder to PARALUS with the `-e` flag.

### Step 4 — Run PARALUS

**Option A — using a pre-built annotation BED file:**

```bash
bash paralus_1.7.1.sh \
  -i sample.fastq.gz \
  -x /path/to/genome_index \
  -A /path/to/annotation_Paralus.bed \
  -f /path/to/genome.fa \
  -e /path/to/rRNA_index \
  -a NNTGACTGTGGAATTCTCGGGTGCCAAGG \
  -u 4 \
  -p 12
```

**Option B — using a GTF directly (Trierarch runs automatically):**

```bash
bash paralus_1.7.1.sh \
  -i sample.fastq.gz \
  -x /path/to/genome_index \
  -G /path/to/gencode.gtf \
  -T /path/to/hg38-tRNAs.bed \
  -C \
  -f /path/to/genome.fa \
  -e /path/to/rRNA_index \
  -a NNTGACTGTGGAATTCTCGGGTGCCAAGG \
  -u 4 \
  -p 12
```


---

## Parameters

### Required

| Flag | Description |
|------|-------------|
| `-i` | Input file (FASTQ.gz / FASTQ / FASTA) |
| `-x` | Path to STAR genome index |
| `-A` | Annotation BED file generated by Trierarch — required unless `-G` is used |
| `-f` | Genome FASTA file |

### Optional — Preprocessing

| Flag | Description | Default |
|------|-------------|---------|
| `-a` | 3' adapter sequence | none |
| `-d` | 5' adapter sequence | `GTTCAGAGTTCTACAGTCCGACGATC` |
| `-u` | Bases to trim from 5' end (UMIs) | `0` |
| `-n` | Minimum read length after trimming | `19` |
| `-p` | Threads | `8` |
| `-e` | Path to rRNA genome index for pre-depletion | none |

### Optional — Annotation

| Flag | Description |
|------|-------------|
| `-G` | Path to GTF file — Trierarch will build the annotation BED before alignment (use instead of `-A`) |
| `-T` | tRNA BED file from GtRNAdb (passed to Trierarch) |
| `-C` | Use canonical transcripts only (passed to Trierarch, recommended) |

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
| `*.clusters.csv` | Main output — annotated binding clusters with crosslink statistics, gene name, feature, biotype, and binding sequence |
| `*.clusters.bed` | BED format version of the clusters file |
| `*.aligned_TtoC.bam` | BAM file containing only T-to-C crosslink-supporting reads |
| `*.aligned_TtoC.bam.bai` | Index for the TtoC BAM file |
| `*.aligned_TtoC.bed` | BED file of all T-to-C crosslink-supporting reads |
| `*.Aligned.sortedByCoord.out.bam` | Full genome-aligned BAM file |
| `*_stats.txt` | Per-sample run statistics (read counts, cluster counts by feature type, k-mer enrichment summary, timing) |
| `*_5mers.txt` | 5-mer counts in bound vs. background sequences |
| `*_5mer_enrichment.pdf` | Scatter plot of 5-mer enrichment at binding sites vs. background (requires matplotlib) |
| `*_ttoc_utr_cds.pdf` | Metagene plot of T-to-C crosslink site distribution across 5'UTR, CDS, and 3'UTR (requires matplotlib) |
| `*_ttoc_intron.pdf` | Metagene plot of T-to-C crosslink site distribution across introns (requires matplotlib) |

---

## Repository structure

```
PARALUS/
├── paralus_1.7.1.sh         # Main pipeline script (includes Trierarch)
├── INSTALLATION.md          # Detailed dependencies installation instructions for all environments
└── README.md
```

---

## Citation

> **[Manuscript in preparation]**
> Konstantinidou P, Benhalevy D, Sarshad AA, Anastasakis DG. PARALUS: A Lightweight Unix Shell Pipeline for PAR-CLIP and CLIP Sequencing Data Analysis. *Journal TBD*, year.

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
