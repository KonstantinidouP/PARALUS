# PARALUS - Reference Setup

This guide explains how to download the reference files PARALUS needs and build the genome and rRNA indexes. You only need to do this once per genome. The examples use the human genome (hg38). For installing PARALUS and its dependencies, see [INSTALLATION.md](INSTALLATION.md). For running PARALUS, see the main [README.md](README.md).

## Contents

- [1. Reference genome and genome index](#1-reference-genome-and-genome-index)
- [2. Gene annotation](#2-gene-annotation)
- [3. rRNA index (optional)](#3-rrna-index-optional)

---

## 1. Reference genome and genome index

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

---

## 2. Gene annotation

Trierarch requires a gene annotation GTF file and optionally a tRNA BED file and/or a miRNA GFF3 file. The GTF can be from GENCODE or Ensembl. All reference files are freely available from public databases.

We recommend downloading the latest comprehensive gene annotation GTF for your organism from **GENCODE** (https://www.gencodegenes.org) or **Ensembl** (https://www.ensembl.org), and optionally the tRNA annotation BED file from the **Genomic tRNA Database (GtRNAdb)** (https://gtrnadb.ucsc.edu) and/or the miRNA annotation GFF3 from **miRBase** (https://www.mirbase.org/download/). Download links and instructions for various organisms can be found in the PARALUS help page. Run the following command to view:

```bash
paralus --trierarch --help
```

### Build the annotation BED file with Trierarch

Trierarch is built into PARALUS and produces the annotation BED file from the GTF. Using canonical transcripts only (the `-C` option) is recommended for cleaner annotation:

```bash
paralus --trierarch \
  -G gencode.v49.annotation.gtf \
  -T hg38-tRNAs.bed \
  -C \
  -O gencode.v49.annotation_Trierarch_canonical
```

This produces `gencode.v49.annotation_Trierarch_canonical_Paralus.bed` to use with PARALUS via the `-A` flag. Alternatively, use `-G` directly in your PARALUS run and Trierarch builds the annotation before alignment.

---

## 3. rRNA index (optional)

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
