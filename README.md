<h1 align="center">PARALUS</h1>

<p align="center"><b>A Lightweight Unix Shell Pipeline for PAR-CLIP and CLIP Sequencing Data Analysis</b></p>

<p align="center">
  <a href="https://github.com/KonstantinidouP/PARALUS/releases/latest"><img src="https://img.shields.io/github/v/release/KonstantinidouP/PARALUS?color=2a9d8f" alt="Latest release"></a>
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS-blue" alt="Platform">
  <img src="https://img.shields.io/badge/shell-bash-green" alt="Shell">
  <img src="https://img.shields.io/badge/data-PAR--CLIP%20%7C%20CLIP-e76f51" alt="Data">
</p>

PARALUS is a bash-based pipeline for the analysis of PAR-CLIP and CLIP sequencing data. It performs end-to-end processing from raw reads to annotated RNA-binding protein (RBP) binding sites, and includes k-mer enrichment analysis to identify sequence preferences.

---

## Overview

PARALUS takes CLIP or PAR-CLIP sequencing data as input and outputs RBP binding sites (clusters) with full annotation, crosslink site information, and binding statistics. The pipeline handles adapter trimming, duplicate read collapsing, genome alignment, T-to-C crosslink site extraction (relevant for PAR-CLIP data), peak calling with read distribution modeling, and genomic annotation - all in a single run.

PARALUS includes **Trierarch**, a companion annotation module that generates a custom annotation BED file directly from a Gencode/Ensembl GTF file. Trierarch is embedded in the same script and can be run standalone using the `--trierarch` mode, or triggered within a PARALUS run by providing a GTF file with the `-G` flag.

<p align="center">
  <img src="images/paralus_workflow.svg" width="800" alt="PARALUS workflow: input, preprocessing, alignment, crosslink sites, cluster calling and outputs, with Trierarch annotation">
</p>

---

## Features

- Handles multiple input formats: FASTQ, FASTQ.gz, and FASTA
- Optional UMI trimming and adapter removal via cutadapt
- Optional rRNA depletion by pre-mapping to an rRNA genome index
- T-to-C crosslink site extraction from PAR-CLIP data at single-nucleotide resolution
- Two-stage peak calling based on read distribution modeling that filters out artifacts and background
- Gene name, feature, biotype, read counts, and T-to-C conversion metrics per cluster
- Resolution of overlapping annotations by priority-ranked feature assignment
- Built-in annotation module (Trierarch); supports canonical transcripts only (`-C`), optional tRNA (`-T`) and miRNA (`-M`) databases
- 5-mer enrichment analysis at binding sites vs. background, with optional PDF scatter plot output
- Metagene analysis plots for UTR/CDS and intron binding distribution (optional PDF output)
- Full run logging to a timestamped log file
- Per-sample statistics report including binding breakdown by gene feature type

---

## Performance & Accessibility

PARALUS is designed to be fast, self-contained, and accessible - running on any Linux or Mac machine with a single command. Genome alignment with STAR requires sufficient RAM to load the genome index (~30 GB for human hg38). A minimum of 32 GB RAM is therefore recommended for human or mouse genome analysis on a local machine with STAR.

**Low-RAM option.** On machines with limited RAM (less than 32 GB), use the `-H` flag to align with HISAT2 instead of STAR. See [INSTALLATION.md](INSTALLATION.md) for installing HISAT2 and building a HISAT2 index, and the `-H` entry in `paralus --help` for details. HISAT2 needs far less memory and recovers the same biological signal: on CPSF1, both aligners showed the 3'UTR binding spike and enrichment of the AATAAA polyadenylation signal. HISAT2 reports fewer clusters (845 vs. 1,241 with STAR), so STAR is recommended when memory allows and maximum sensitivity matters.

PARALUS was benchmarked on three publicly available fPAR-CLIP datasets on an Apple MacBook Pro (M4 Pro, 48 GB RAM) using 4 threads (`-p 4`):

| Dataset | Input FASTQ.gz size | Raw reads | Unique reads (after collapsing) | Runtime (minutes) | CPU |
|---|---|---|---|---|---|
| CPSF1 (SRR13313335) | 238 MB | 9,135,633 | 2,142,567 | 2:29.85 (~2.5 min) | 113% |
| KHSRP (SRR13313340) | 559 MB | 19,992,492 | 8,878,651 | 7:06.72 (~7.1 min) | 173% |
| U2AF1 (SRR8606876) | 1.4 GB | 54,284,359 | 10,680,775 | 14:21.47 (~14.4 min) | 147% |

Runtime scales with input size and library complexity (unique reads remaining after collapsing duplicates), ranging from ~2.5 minutes for a 238 MB dataset to ~14.4 minutes for a 1.4 GB dataset. CPU utilization of 113-173% shows that part of the run uses multiple cores. All analyses completed on a local Mac laptop without requiring access to an HPC cluster.

> **Note:** Alignment and sorting steps benefit from additional threads (set with `-p`), but some steps such as peak calling and 5-mer counting are inherently serial. On local machines, disk speed and available RAM can be as important as CPU core count. Running PARALUS on a Linux server can be beneficial when analyzing large datasets or batch processing multiple samples.

---

## Dependencies

PARALUS requires the following tools installed in your environment:

- [STAR](https://github.com/alexdobin/STAR) (tested with v2.7.11b) **or** [HISAT2](http://daehwankimlab.github.io/hisat2/) (alternative aligner, selected with `-H` - recommended if you have less than 32 GB RAM)
- [samtools](http://www.htslib.org/) (tested with v1.22.1)
- [bedtools](https://bedtools.readthedocs.io/) (tested with v2.31.1)
- [cutadapt](https://cutadapt.readthedocs.io/)
- python3 (OPTIONAL - required for PDF plot output)
- matplotlib (OPTIONAL - required for PDF plot output)
- numpy (OPTIONAL - required for PDF plot output)

For step-by-step installation instructions across different environments (Apple Silicon Mac, Linux machines/servers, HPC clusters), see [INSTALLATION.md](INSTALLATION.md).

---

## Quick Start

### Step 1 - Install PARALUS and its dependencies

Clone the repository to any location on your machine, then add the PARALUS folder to your `PATH` (the list of folders your terminal searches for programs), so you can run PARALUS as `paralus` from any folder. The last three lines add two lines to your shell settings file, a label and the `PATH` setting. Run them from inside the PARALUS folder:

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

Then install the dependencies into a conda environment called `paralus`. The installer checks which tools you already have and asks before installing anything:

```bash
bash install_dependencies.sh
```

Activate the environment before each run:

```bash
conda activate paralus
```

For manual installation, HPC clusters and other details, see [INSTALLATION.md](INSTALLATION.md).

### Step 2 - Build the annotation file (run once per genome)

Trierarch requires a gene annotation GTF file and optionally a tRNA BED file and/or a miRNA GFF3 file. The GTF can be from GENCODE or Ensembl. All reference files are freely available from public databases.

We recommend downloading the latest comprehensive gene annotation GTF for your organism from **GENCODE** (https://www.gencodegenes.org) or **Ensembl** (https://www.ensembl.org), and optionally the tRNA annotation BED file from the **Genomic tRNA Database (GtRNAdb)** (https://gtrnadb.ucsc.edu) and/or the miRNA annotation GFF3 from **miRBase** (https://www.mirbase.org/download/). Download links and instructions for various organisms can be found in the PARALUS help page. Run the following command to view:

```bash
paralus --trierarch --help
```

#### Build the annotation BED file with Trierarch

Trierarch is built into PARALUS and produces the annotation BED file from the GTF. Using canonical transcripts only (the `-C` option) is recommended for cleaner annotation:

```bash
paralus --trierarch \
  -G gencode.v49.annotation.gtf \
  -T hg38-tRNAs.bed \
  -C \
  -O gencode.v49.annotation_Trierarch_canonical
```

This produces `gencode.v49.annotation_Trierarch_canonical_Paralus.bed` to use with PARALUS via the `-A` flag. Alternatively, use `-G` directly in your PARALUS run and Trierarch builds the annotation before alignment.

### Step 3 - Download the reference genome and build a STAR index (run once per genome)

You need the reference genome of your organism of interest in FASTA format. If you don't have it yet, download it from UCSC. These are large files (~1 GB compressed) - make sure you have enough disk space. Below is an example for the human genome (hg38).

Linux:

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
```

Mac (curl is pre-installed):

```bash
curl -o hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
```

Decompress the FASTA (decompressed file will be ~3.5 GB for hg38) and build the STAR genome index. This is a one-time step per genome and requires ~30 GB of disk space and significant RAM (≥32 GB recommended):

```bash
gunzip -k hg38.fa.gz
mkdir /path/to/genome_index
STAR --runMode genomeGenerate \
     --genomeDir /path/to/genome_index \
     --genomeFastaFiles /path/to/hg38.fa \
     --runThreadN 6 \
     --limitGenomeGenerateRAM 38000000000
```

> **Note:** The `--limitGenomeGenerateRAM` flag caps memory usage in bytes. Adjust based on your available RAM - on a server with ≥64 GB RAM you can omit this flag entirely.

Pass the genome index folder to PARALUS with the `-x` flag, and the uncompressed FASTA file with the `-f` flag.

<details>
<summary>💡 <b>Low-RAM option with HISAT2 (for local machines with less than 32 GB RAM)</b></summary>

Instead of the STAR index above, you can build a HISAT2 index. It needs ~8 GB RAM to build, ~5 GB RAM to align and ~4.5 GB of disk for hg38:

```bash
mkdir -p /path/to/genome_index
hisat2-build /path/to/hg38.fa /path/to/genome_index/genome
```

Pass the resulting index basename (`/path/to/genome_index/genome`, not a directory) to PARALUS with `-x`, together with the `-H` flag to switch PARALUS to HISAT2 mode. See [INSTALLATION.md](INSTALLATION.md) for installing HISAT2.

</details>

### Step 3b - Download rRNA sequences and build an rRNA index (OPTIONAL - run once per genome)

This step is only needed if you want to use the rRNA depletion option in PARALUS (the `-e` flag). It removes rRNA-mapping reads before genome alignment, which can significantly improve the signal-to-noise ratio in your data.

Download the human rRNA sequences:

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

Alternatively, both sequences can be downloaded manually from the **NCBI Nucleotide database** (https://www.ncbi.nlm.nih.gov/nucleotide) - search for accessions `U13369.1` (45S rDNA repeat unit) and `X12811.1` (5S rRNA), download each as FASTA.

Then concatenate the two files and build a STAR index:

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

<details>
<summary>💡 <b>Low-RAM option with HISAT2 (for local machines with less than 32 GB RAM)</b></summary>

Instead of the STAR rRNA index above, you can build a HISAT2 one:

```bash
mkdir -p /path/to/rRNA_index
hisat2-build rRNA.fa /path/to/rRNA_index/rRNA
```

Pass the resulting basename (`/path/to/rRNA_index/rRNA`) to PARALUS with `-e`.

</details>

### Step 4 - Run PARALUS

Make sure all dependencies are loaded or activated before running. See [INSTALLATION.md](INSTALLATION.md) for instructions tailored to your system.

**Option A - using a pre-built annotation BED file:**

```bash
paralus \
  -i /path/to/sample.fastq.gz \
  -x /path/to/genome_index \
  -A /path/to/annotation_Paralus.bed \
  -f /path/to/genome.fa \
  -e /path/to/rRNA_index \
  -a NNTGACTGTGGAATTCTCGGGTGCCAAGG \
  -u 4 \
  -p 12
```

**Option B - using a GTF directly (Trierarch builds the annotation first):**

```bash
paralus \
  -i /path/to/sample.fastq.gz \
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
| `-i` | Path to the input file (FASTQ.gz / FASTQ / FASTA). Output files are written to the same folder. |
| `-x` | Path to genome index (STAR genomeDir by default; a HISAT2 index basename if `-H` is set) |
| `-A` | Path to the annotation BED file generated by Trierarch - required unless `-G` is used |
| `-f` | Path to the genome FASTA file (uncompressed) |

### Preprocessing

| Flag | Description | Default |
|------|-------------|---------|
| `-a` | 3' adapter sequence | none |
| `-d` | 5' adapter sequence | `GTTCAGAGTTCTACAGTCCGACGATC` |
| `-u` | Bases to trim from 5' end (UMIs) | `0` |
| `-n` | Minimum read length after trimming | `19` |
| `-p` | Threads | `8` |
| `-e` | Path to rRNA genome index for pre-depletion (STAR genomeDir, or HISAT2 basename if `-H` is set) | none |
| `-H` | Use HISAT2 instead of STAR for all alignment steps. Option for low-RAM machines. HISAT2 needs ~8 GB RAM to build the index and ~5 GB to align. `-x`/`-e` must then be HISAT2 index basenames (see [INSTALLATION.md](INSTALLATION.md)) | off |

### Annotation

The following flags are only relevant when `-G` is set:

| Flag | Description |
|------|-------------|
| `-G` | Path to GTF file - Trierarch will build the annotation BED before alignment (replaces `-A`) |
| `-T` | tRNA BED file from GtRNAdb - used with `-G` to include tRNA annotations (optional) |
| `-M` | miRNA GFF3 file from miRBase, or a BED6 file in the same format as `-T` - used with `-G` to include miRNA annotations (optional) |
| `-C` | Use canonical/MANE transcripts only - used with `-G`, recommended for GENCODE/Ensembl |
| `-U` | Keep genes with no `gene_name` attribute, using their Ensembl accession (e.g. `ENSG...`) as the display name instead - used with `-G`. Off by default; without it, unnamed genes are excluded from the annotation BED and any clusters overlapping them appear as unannotated. |

### Cluster calling

| Flag | Description | Default |
|------|-------------|---------|
| `-m` | Minimum reads to define a group | `3` |
| `-c` | Minimum T-to-C conversion locations per cluster | `1` |
| `-r` | Minimum fraction of T-to-C reads per cluster | `0.3` |

> **CLIP mode:** Set both `-c 0` and `-r 0` to run PARALUS in CLIP mode (no T-to-C signal expected). In this mode, T-to-C extraction is skipped and all groups are promoted to clusters. Both flags must be set to `0` together - setting only one to `0` will produce an error.

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

**Terminal output** - summary printed at the end of the CPSF1 run (`paralus_1.8.7.sh`, STAR, `-p 4`):

```
PARALUS finished successfully
  Binding sites (clusters): 1241
  Run time: 2m 28s
  Results in: /path/to/CPSF1_fromGEO_SRR13313335
```

**5-mer enrichment scatter plot** - each of the 1024 pentamers is plotted by its count in bound sequences (y-axis) vs. background sequences (x-axis). Pentamers enriched at binding sites fall above the diagonal; depleted ones fall below. The AATAAA polyadenylation signal is represented among the top enriched 5-mers.

![5-mer enrichment](images/FH-CPSF1_fPARCLIP_0.5Uul_RNaseI_5mer_enrichment.png)

**Metagene plot** - distribution of T-to-C crosslink sites across the 5'UTR, CDS, and 3'UTR of protein-coding genes, normalized by region length. A distinct spike at the 3' end of the 3'UTR reflects CPSF1's binding at the cleavage and polyadenylation site.

![Metagene UTR/CDS](images/FH-CPSF1_fPARCLIP_0.5Uul_RNaseI_ttoc_utr_cds.png)

---

## Repository structure

```
PARALUS/
├── paralus_1.8.7.sh         # Main pipeline script (includes Trierarch)
├── paralus                  # Launcher: runs the pipeline script as `paralus`
├── install_dependencies.sh  # Dependency installer (Linux and Apple Silicon Mac)
├── INSTALLATION.md          # Detailed dependencies installation instructions for all environments
├── PARALUS_manual.md        # Extended manual with full parameter reference and troubleshooting
├── images/                  # Workflow figure and example output plots (CPSF1) shown in this README
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
