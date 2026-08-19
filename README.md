# PARALUS

> **A Lightweight Unix Shell Pipeline for PAR-CLIP and CLIP Sequencing Data Analysis**

PARALUS is a bash-based pipeline for the analysis of PAR-CLIP and CLIP sequencing data. It performs end-to-end processing from raw reads to annotated crosslink clusters, and includes k-mer enrichment analysis to identify sequence preferences at RNA-protein binding sites.

---

## Overview

PARALUS takes CLIP or PAR-CLIP sequencing data as input and outputs high-confidence RNA-protein binding clusters with full annotation, crosslink site information, and binding statistics. The pipeline handles adapter trimming, read collapsing, genome alignment, T-to-C crosslink site extraction (relevant for PAR-CLIP data), peak calling with read distribution modeling, and genomic annotation - all in a single run.

PARALUS includes **Trierarch**, a built-in annotation builder that generates the custom annotation BED file directly from a Gencode/Ensembl GTF file. Trierarch is embedded in the same script and can be run standalone using the `--trierarch` mode, or triggered within a PARALUS run by providing a GTF file with the `-G` flag.

---

## Features

- Handles multiple input formats: FASTQ, FASTQ.gz, and FASTA
- Optional UMI trimming and adapter removal via cutadapt
- Optional rRNA depletion by pre-mapping to an rRNA genome index
- T-to-C crosslink site extraction from PAR-CLIP data at single-nucleotide resolution
- Two-stage peak calling based on read distribution modeling (normal and non-normal clusters)
- Gene name, feature, biotype, read counts, and T-to-C conversion metrics per cluster
- Resolution of overlapping annotations by priority-ranked feature assignment
- Built-in annotation builder (Trierarch); supports canonical transcripts only (`-C`), optional tRNA (`-T`) and miRNA (`-M`) databases
- 5-mer enrichment analysis at binding sites vs. background, with optional PDF scatter plot output
- Metagene analysis plots for UTR/CDS and intron binding distribution (optional PDF output)
- Full run logging to a timestamped log file
- Per-sample statistics report including binding breakdown by gene feature type

---

## Performance & Accessibility

PARALUS is designed to be fast, self-contained, and accessible - running on both Linux servers and local Mac machines with a single command.

On an Apple MacBook Pro (M4 Pro, 48 GB RAM) using 4 threads, PARALUS was benchmarked on two datasets of different sizes:

| Dataset size | Approx. reads | Runtime | CPU utilization |
|---|---|---|---|
| 512 MB FASTQ.gz | ~8.8 million | ~7 minutes | 171% |
| 4.7 GB FASTQ.gz | ~10x larger | ~100 minutes | 111% |

Runtime is consistent and reproducible across runs. Alignment and sorting steps benefit from additional threads (set with `-p`), but some steps such as peak calling and 5-mer counting are inherently serial. On local machines, disk I/O speed and available RAM can be as important as CPU core count - a fast SSD and sufficient RAM to cache large files like the genome FASTA in memory will improve performance significantly. For large datasets or batch processing of multiple samples, running on a Linux server is recommended.

---

## Dependencies

PARALUS requires the following tools installed in your environment:

- [STAR](https://github.com/alexdobin/STAR) (tested with v2.7.2b)
- [samtools](http://www.htslib.org/) (tested with v1.17)
- [bedtools](https://bedtools.readthedocs.io/) (tested with v2.29.2)
- [cutadapt](https://cutadapt.readthedocs.io/)
- python3 (OPTIONAL - required for PDF plot output)
- matplotlib (OPTIONAL - required for PDF plot output)
- numpy (OPTIONAL - required for PDF plot output)

For step-by-step installation instructions across different environments (Apple Silicon Mac, Linux machines/servers, HPC clusters), see [INSTALLATION.md](INSTALLATION.md).

---

## Quick Start

### Step 1 - Clone the repository and activate the environment

Clone the repository to any location on your machine - PARALUS is a bash script and does not need to be installed in the conda environment:

```bash
git clone https://github.com/KonstantinidouP/PARALUS.git
cd PARALUS
chmod +x paralus_1.7.3.sh
```

To view all available options:

```bash
./paralus_1.7.3.sh --help
```

Before running PARALUS, make sure all dependencies are available in your environment. See [INSTALLATION.md](INSTALLATION.md) for details on how to install and activate the dependencies for different systems.

### Step 2 - Build the annotation file (run once per genome)

Trierarch requires a gene annotation GTF file and optionally a tRNA BED file and/or a miRNA GFF3 file. The GTF can be from GENCODE or Ensembl. All reference files are freely available from public databases.

We recommend downloading the latest comprehensive gene annotation GTF for your organism from **GENCODE** (https://www.gencodegenes.org) or **Ensembl** (https://www.ensembl.org), and optionally the tRNA annotation BED file from the **Genomic tRNA Database (GtRNAdb)** (https://gtrnadb.ucsc.edu) and/or the miRNA annotation GFF3 from **miRBase** (https://www.mirbase.org/download/). Download links and instructions for various organisms can be found in the PARALUS help page. Run the following command to view:

```bash
bash paralus_1.7.3.sh --trierarch --help
```

#### Build the annotation BED file with Trierarch

Trierarch is built into PARALUS and produces the annotation BED file from the GTF. Using canonical transcripts only (the `-C` option) is recommended for cleaner annotation:

```bash
bash paralus_1.7.3.sh --trierarch \
  -G gencode.v49.annotation.gtf \
  -T hg38-tRNAs.bed \
  -C \
  -O gencode.v49.annotation_Trierarch_canonical
```

This produces `gencode.v49.annotation_Trierarch_canonical_Paralus.bed` to use with PARALUS via the `-A` flag. Alternatively, use `-G` directly in your PARALUS run and Trierarch will build the annotation automatically.

### Step 3 - Download the reference genome and build a STAR index (run once per genome)

You need the reference genome of your organism of interest in FASTA format. If you don't have it yet, download it from UCSC. These are large files (~1 GB compressed) - make sure you have enough disk space.

For human (hg38) — Linux:

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
```

Mac (curl is pre-installed):

```bash
curl -o hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
```

Decompress the FASTA - STAR requires an uncompressed file for genome index generation:

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

> **Note:** The `--limitGenomeGenerateRAM` flag caps memory usage in bytes. Adjust based on your available RAM - on a server with ≥64 GB RAM you can omit this flag entirely.

Pass the genome index folder to PARALUS with the `-x` flag, and the uncompressed FASTA file with the `-f` flag.

### Step 3b - Download rRNA sequences and build an rRNA index (OPTIONAL - run once per genome)

This step is only needed if you want to use the rRNA depletion option in PARALUS (the `-e` flag). It removes rRNA-mapping reads before genome alignment, which can significantly improve the signal-to-noise ratio in your data.

Download the human rRNA sequences and build a STAR index. On Linux, `wget` is available by default. On Mac, use `curl` (pre-installed) or install `wget` via `brew install wget`:

Linux:

```bash
wget -O U13369.1.fasta "https://www.ebi.ac.uk/ena/browser/api/fasta/U13369.1?download=true"
wget -O X12811.1.fasta "https://www.ebi.ac.uk/ena/browser/api/fasta/X12811.1?download=true"
```

Mac:

```bash
curl -o U13369.1.fasta "https://www.ebi.ac.uk/ena/browser/api/fasta/U13369.1?download=true"
curl -o X12811.1.fasta "https://www.ebi.ac.uk/ena/browser/api/fasta/X12811.1?download=true"
```

Alternatively, both sequences can be downloaded manually from the **NCBI Nucleotide database** (https://www.ncbi.nlm.nih.gov/nucleotide) - search for accessions `U13369.1` (45S rDNA repeat unit) and `X12811.1` (5S rRNA), download each as FASTA, then run:

```bash
cat U13369.1.fasta X12811.1.fasta > rRNA.fa
mkdir /path/to/rRNA_index
STAR --runMode genomeGenerate \
     --genomeDir /path/to/rRNA_index \
     --genomeFastaFiles rRNA.fa \
     --runThreadN 6 \
     --genomeSAindexNbases 6
```

> **Note:** The `--genomeSAindexNbases 6` flag is required for small genomes - refer to the STAR manual for details. Pass the rRNA index folder to PARALUS with the `-e` flag.

### Step 4 - Run PARALUS

Make sure all dependencies are loaded or activated before running. See [INSTALLATION.md](INSTALLATION.md) for instructions tailored to your system.

**Option A - using a pre-built annotation BED file:**

```bash
bash paralus_1.7.3.sh \
  -i sample.fastq.gz \
  -x /path/to/genome_index \
  -A /path/to/annotation_Paralus.bed \
  -f /path/to/genome.fa \
  -e /path/to/rRNA_index \
  -a NNTGACTGTGGAATTCTCGGGTGCCAAGG \
  -u 4 \
  -p 12
```

**Option B - using a GTF directly (Trierarch runs automatically):**

```bash
bash paralus_1.7.3.sh \
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
| `-x` | Path to STAR genome index directory |
| `-A` | Annotation BED file generated by Trierarch - required unless `-G` is used |
| `-f` | Genome FASTA file |

### Preprocessing

| Flag | Description | Default |
|------|-------------|---------|
| `-a` | 3' adapter sequence | none |
| `-d` | 5' adapter sequence | `GTTCAGAGTTCTACAGTCCGACGATC` |
| `-u` | Bases to trim from 5' end (UMIs) | `0` |
| `-n` | Minimum read length after trimming | `19` |
| `-p` | Threads | `8` |
| `-e` | Path to rRNA genome index for pre-depletion | none |

### Annotation

The following flags are only relevant when `-G` is set:

| Flag | Description |
|------|-------------|
| `-G` | Path to GTF file - Trierarch will build the annotation BED before alignment (replaces `-A`) |
| `-T` | tRNA BED file from GtRNAdb - used with `-G` to include tRNA annotations (optional) |
| `-M` | miRNA GFF3 file from miRBase, or a BED6 file in the same format as `-T` - used with `-G` to include miRNA annotations (optional) |
| `-C` | Use canonical/MANE transcripts only - used with `-G`, recommended for GENCODE/Ensembl |

### Cluster calling

| Flag | Description | Default |
|------|-------------|---------|
| `-m` | Minimum reads to define a group | `3` |
| `-c` | Minimum T-to-C conversion locations per cluster | `1` |
| `-r` | Minimum fraction of T-to-C reads per cluster | `0.3` |

> **CLIP mode:** Set both `-c 0` and `-r 0` to run PARALUS in CLIP mode (no T-to-C signal expected). In this mode, T-to-C extraction is skipped and all groups are promoted to clusters. Both flags must be set to `0` together — setting only one to `0` will produce an error.

---

## Output files

| File | Description |
|------|-------------|
| `*.Aligned.sortedByCoord.out.bam` (+`.bai`) | Full genome alignment. |
| `*.aligned_TtoC.bam` (+`.bai`) | Subset of aligned reads carrying a T-to-C conversion. Empty in CLIP mode. |
| `*.aligned_TtoC.bed` | T-to-C crosslink-supporting reads in BED format. Empty in CLIP mode. |
| `*.groups.csv` | Every candidate group, before cluster-threshold filtering (`-c`/`-r`). |
| `*.clusters.csv` | Final annotated PAR-CLIP / CLIP binding sites: genomic coordinates, annotated gene name, feature type (e.g. 3'UTR, CDS, intron), biotype (e.g. protein_coding, lncRNA), binding sequence, read counts, and T-to-C conversion metrics, per row. |
| `*.clusters.bed` | BED version of the clusters file. |
| `*_run_<timestamp>.log` | Full stdout/stderr of the run, plus the exact invoked command. |
| `*_stats.txt` | Run parameters, total runtime, and summary counts (reads per biotype/feature, cluster counts, etc). |
| `*_5mers.txt` | Raw bound and background 5-mer counts. |
| `*_5mer_ratios.tsv` | Bound vs. background enrichment ratio per 5-mer, ranked from most to least enriched. |
| `*_5mer_enrichment.pdf` | Scatter plot illustrating 5-mer enrichment. Only generated if python3 with numpy and matplotlib are installed. |
| `*_ttoc_utr_cds.pdf` | Metagene plot of crosslinked read distribution across 5'UTR → CDS → 3'UTR. In CLIP mode, this plot shows the distribution of all aligned reads instead. Only generated if python3 with numpy and matplotlib are installed. |
| `*_ttoc_intron.pdf` | Metagene plot of crosslinked read distribution across introns (5' splice site → 3' splice site). In CLIP mode, this plot shows the distribution of all aligned reads instead. Only generated if python3 with numpy and matplotlib are installed. |
| `rRNA_aligned_*.bam` | Only if `-e` was used: reads that mapped to the rRNA reference. Retained for cases where rRNA binding is of biological interest. |

---

## Example Output

The following example shows PARALUS output from a PAR-CLIP experiment with CPSF1 (Cleavage and Polyadenylation Specificity Factor 1), using publicly available data (SRR13313335; Anastasakis et al., *Nucleic Acids Research* 2021, PMID: 33503264). CPSF1 is a component of the cleavage and polyadenylation complex known to recognize the canonical AAUAAA polyadenylation signal sequence on precursor messenger RNA (pre-mRNA). The predominant binding at 3'UTRs is consistent with its known function, and the enrichment of AATAAA (the DNA equivalent of AAUAAA) is clearly visible in the 5-mer scatter plot.

**5-mer enrichment scatter plot** - each of the 1024 pentamers is plotted by its count in bound sequences (y-axis) vs. background sequences (x-axis). Pentamers enriched at binding sites fall above the diagonal; depleted ones fall below. The AATAAA polyadenylation signal is represented among the top enriched 5-mers.

![5-mer enrichment](images/FH-CPSF1_fPARCLIP_0.5Uul_RNaseI_5mer_enrichment.png)

**Metagene plot** - distribution of T-to-C crosslink sites across the 5'UTR, CDS, and 3'UTR of protein-coding genes, normalized by region length. A distinct spike at the 3' end of the 3'UTR reflects CPSF1's binding at the cleavage and polyadenylation site.

![Metagene UTR/CDS](images/FH-CPSF1_fPARCLIP_0.5Uul_RNaseI_ttoc_utr_cds.png)

---

## Repository structure

```
PARALUS/
├── paralus_1.7.3.sh         # Main pipeline script (includes Trierarch)
├── INSTALLATION.md          # Detailed dependencies installation instructions for all environments
├── PARALUS_manual.md        # Extended manual with full parameter reference and troubleshooting
└── README.md
```

---

## Citation

> **[Manuscript in preparation]**
> Konstantinidou P, Benhalevy D, Sarshad AA, Anastasakis DG. PARALUS: A Lightweight Unix Shell Pipeline for PAR-CLIP and CLIP Sequencing Data Analysis. *Journal TBD*, year.

If you use PARALUS in your work, please cite the above once published.

---

## Authors

- **Dimitrios Anastasakis** - lead developer
- **Parthena Konstantinidou** - co-developer

---

## License

> **[Placeholder - license to be decided]**
>
> For now, this repository is private. License will be added prior to public release.
