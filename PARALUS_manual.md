# PARALUS Extended Manual

This is the full manual for PARALUS and its companion annotation
builder, Trierarch (`paralus_1.7.3.sh`). For a guided first run, see
`README.md` instead - this document is for looking things up,
not for a first-time walkthrough.

---

## Table of contents

1. [Overview](#1-overview)
2. [Trierarch](#2-trierarch)
3. [PARALUS](#3-paralus)
4. [Pipeline walkthrough](#4-pipeline-walkthrough)
5. [Output files](#5-output-files)
6. [Advanced options](#6-advanced-options)
7. [Troubleshooting](#7-troubleshooting)

---

## 1\. Overview

PARALUS is a single Bash file that performs every step of PAR-CLIP/CLIP
analysis: adapter/UMI trimming, optional rRNA depletion, genome alignment
(STAR), T-to-C crosslink-site extraction (for PAR-CLIP data), a two-stage
group-to-cluster calling algorithm, genomic annotation of binding clusters,
5-mer sequence enrichment, and metagene analysis of crosslink site distribution across transcript regions -
all from a single command, with no dependencies beyond cutadapt, STAR,
samtools, bedtools, and (for the two plotting steps) python3 with
numpy/matplotlib.

Trierarch, its companion annotation builder, is merged into the same
file as a Bash function. It converts any GTF into the 7-column BED
PARALUS needs, resolving overlapping annotations by an explicit priority order.

Both tools are invoked from the same script:

```bash
./paralus_1.7.3.sh --trierarch -G genes.gtf -C -O prefix   # Trierarch standalone run to pre-build the annotation BED
./paralus_1.7.3.sh -i ... -x ... -A ... [options]         # a standard PARALUS run for PAR-CLIP / CLIP analysis
```

---

## 2\. Trierarch

PARALUS requires a non-overlapping genomic annotation BED file to assign each binding cluster to a gene, feature, and biotype. Its built-in annotation builder, Trierarch, builds this file from a standard GTF annotation, with optional tRNA and miRNA database integration for more complete annotation coverage. The annotation BED only needs to be generated once per genome assembly and can be reused across all standard PARALUS runs using that genome assembly. Trierarch can be run as a standalone step (using `--trierarch` mode), and the resulting BED file can be passed to PARALUS via `-A`. Alternatively, Trierarch can be triggered within a standard PARALUS run by providing a GTF file with the `-G` flag instead of a pre-built BED file.

### 2.1 Flags

|Flag|Required|Description|
|-|-|-|
|`-G genes.gtf`|yes|GTF annotation (GENCODE or Ensembl format). Recognises either the `gene_type` or `gene_biotype` attribute, whichever the file uses.|
|`-T tRNA_DB.bed`|no|tRNA coordinates, bed6 (chr, start, end, name, score, strand) - e.g. from GtRNAdb.|
|`-M miRNA_DB`|no|miRNA coordinates - either a miRBase GFF3 (auto-detected by column count) or a plain bed6 in the same shape as `-T`.|
|`-O output_prefix`|no|Output file prefix. Default: the GTF filename with `.gtf` stripped. If a full path is provided, output files are written to that location; otherwise they are written to the same directory as the GTF. The final output file is always named `<prefix>_Paralus.bed`.|
|`-C`|no|Restrict to canonical transcripts (Ensembl_canonical or MANE_Select tagged) for protein-coding genes. For non-coding genes, where canonical transcript designation is less consistently defined, all transcripts are included regardless of canonical tagging. Off by default, but strongly recommended when using recent GENCODE or Ensembl GTFs, which carry the required canonical transcript tags. Produces cleaner UTR/CDS boundaries for protein-coding genes.|
|`-U`|no|Keep genes with no `gene_name` attribute, falling back to `gene_id` (e.g. an ENSG... accession) as the display name. Off by default. Genes with no `gene_name` are excluded from the annotation, since an accession-only name usually indicates a novel/unannotated locus.|

### 2.2 Annotation construction and overlap resolution

Trierarch uses a two-level resolution strategy to ensure the final annotation BED contains no conflicting annotations.

Overlapping gene type annotations are resolved using the following fixed priority cascade (highest to lowest): `snoRNA > snRNA > rRNA > miRNA > tRNA > Mt_tRNA > scaRNA > lncRNA > protein_coding > everything else`. Where different gene types overlap, only the highest-priority one is retained at each genomic position during the final BED assembly.

Where multiple transcripts of the same gene are included, these are first collapsed into a unified, non-overlapping set of coordinates. Feature type conflicts across transcripts are resolved by priority: CDS takes precedence over UTR, and all exonic features take precedence over introns. For protein-coding genes, the `-C` flag can direct Trierarch to build annotation using the canonical transcript only.

> **Note:** Trierarch requires a GTF file as input. If a GFF3 file is provided instead, the script will stop with a clear error message and suggest a conversion command before re-running.

### 2.3 Downloading external databases (optional)

```bash
# tRNA (GtRNAdb) - human hg38 example
wget https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.bed

# miRNA (miRBase) - human example
wget https://www.mirbase.org/download/hsa.gff3
```

For organisms outside either resource's coverage, Trierarch runs fine without them, relying on whatever small-RNA biotypes the GTF itself annotates.

---

## 3\. PARALUS

### 3.1 Required input

|Flag|Description|
|-|-|
|`-i input_file`|Raw reads file: `.fastq.gz`, `.fastq`, or `.fasta`.|
|`-x genome_index`|Full path to the STAR genome index directory.|
|`-A annotation.bed` **or** `-G genes.gtf`|Either a pre-built Trierarch BED file (`-A`), or a raw GTF (`-G`) which triggers Trierarch to build the required BED file before the run starts.|
|`-f genome.fa`|Full path to the genome FASTA file. Required even if a STAR index already exists, as it is needed for T-to-C crosslink site extraction (`samtools calmd`) and for retrieving the binding site sequences in the output (`bedtools getfasta`).|

### 3.2 Preprocessing settings

|Flag|Default|Description|
|-|-|-|
|`-n`|`19`|Minimum read length after trimming.|
|`-a`|none|3' adapter sequence for cutadapt.|
|`-d`|`GTTCAGAGTTCTACAGTCCGACGATC`|5' adapter sequence. If the 5' adapter contains a barcode, include the full barcoded sequence here.|
|`-u`|`0`|Bases trimmed from the 5' end after adapter removal (UMI length).|
|`-p`|`8`|Number of threads. Do not exceed 12 - STAR itself recommends against more.|
|`-e`|none|Full path to a STAR index for an rRNA reference genome. When provided, reads are first mapped against this reference, the percentage of rRNA reads is reported in the stats file, and those reads are removed before genome alignment.|

### 3.3 Binding-site calling settings

|Flag|Default|Description|
|-|-|-|
|`-m`|`3`|Minimum reads to define a group.|
|`-c`|`1`|Minimum number of distinct T-to-C conversion positions required for a group to be called a cluster. Set both `-c` and `-r` to `0` for CLIP data with no T-to-C signal.|
|`-r`|`0.3`|Minimum fraction of reads within a group that must carry a T-to-C conversion for it to be called a cluster. Set both `-r` and `-c` to `0` for CLIP data with no T-to-C signal.|

### 3.4 Trierarch settings (used with `-G`)

When `-G` is used, the following Trierarch settings (described in Section 2.1) can also be passed directly in the same PARALUS command:

| Flag | Description |
|-|-|
| `-T` | tRNA BED file from GtRNAdb (optional) |
| `-M` | miRNA GFF3 or BED file (optional) |
| `-C` | Use canonical transcript only for protein-coding genes (recommended). See Section 2.1 for full details. |
| `-U` | Keep genes with no `gene_name` attribute, using their Ensembl accession ID as the gene name instead. Without this flag, these genes are excluded from the annotation BED and any clusters overlapping them will appear as unannotated. |

---

## 4\. Pipeline walkthrough

1. **Preprocessing.** Reads are collapsed by unique sequence (with read counts preserved in the FASTA headers for downstream depth weighting), rRNA-depleted (`-e`, optional),
and adapter/UMI-trimmed with cutadapt.
2. **Alignment.** STAR is used to align reads end-to-end without allowing spliced alignments (`--alignIntronMax 1`, `--alignEndsType EndToEnd`) against the reference genome. Intron/exon context is recovered later from the annotation.
3. **Crosslink-site extraction.** T-to-C mismatches — the hallmark of PAR-CLIP, arising when 4SU-incorporated sites are UV-crosslinked and result in T-to-C substitutions during cDNA synthesis — are identified and extracted at single-nucleotide resolution. This step is skipped in CLIP mode (`-c 0 -r 0`).
4. **Group calling.** Overlapping reads are merged into groups. Each group's read distribution is tested for normality. Groups with an irregular read distribution (e.g. broad or multi-modal) are set aside as "NonNormal" and handled separately in the next step.
5. **NonNormal group re-evaluation.** NonNormal groups (with broad or irregular read distributions) are re-evaluated on a representative subsample of their reads. This helps in high-coverage regions where elevated background can obscure a genuine binding signal: subsampling reduces the noise and allows the underlying signal to emerge. Groups with sparse, scattered reads are typically discarded at this stage. The subsampling is deterministic, ensuring reproducible results across runs.
6. **Cluster promotion.** Groups that clear the `-c`/`-r` thresholds are promoted to clusters. In CLIP mode (`-c 0 -r 0`) all groups are promoted to clusters.
7. **Annotation.** Clusters are intersected with the annotation BED file. Each cluster is assigned a single gene name, feature, and biotype. For details on how overlapping annotations are resolved, see Section 2.2 (Trierarch).
8. **5-mer enrichment.** Counts all 1024 possible pentamers (5-mers) within the bound cluster sequences and compares them against the counts of 5-mers within the flanking genomic sequences upstream and downstream of each cluster (used as background). Raw counts are saved to `*_5mers.txt` and enrichment ratios to `*_5mer_ratios.tsv`.
9. **Metagene analysis.** Distribution of crosslinked reads (or all aligned reads in CLIP mode) is visualized across two metagene profiles: the full mature mRNA of protein-coding genes (5'UTR → CDS → 3'UTR as a continuous unit) and independently across introns. Each mRNA and intron is normalized by its total crosslinked reads count before averaging, so highly-bound loci do not dominate the profile.

---

## 5\. Output files

PARALUS produces the following output files. Temporary intermediate files created during the run are removed upon completion.

|File|Contents|
|-|-|
|`<name>.Aligned.sortedByCoord.out.bam` (+`.bai`)|Full genome alignment.|
|`<name>.aligned_TtoC.bam` (+`.bai`)|Subset of aligned reads carrying a T-to-C conversion. Empty in CLIP mode.|
|`<name>.aligned_TtoC.bed`|T-to-C crosslink-supporting reads in BED format. Empty in CLIP mode.|
|`<name>.groups.csv`|Every candidate group, before cluster-threshold filtering (`-c`/`-r`).|
|`<name>.clusters.csv`|Final annotated PAR-CLIP / CLIP binding sites: genomic coordinates, annotated gene name, feature type (e.g. 3'UTR, CDS, intron), biotype (e.g. protein_coding, lncRNA), binding sequence, read counts, and T-to-C conversion metrics, per row.|
|`<name>.clusters.bed`|Bed version of the same clusters.|
|`<name>_run_<timestamp>.log`|Full stdout/stderr of the run, plus the exact invoked command.|
|`<name>_stats.txt`|Run parameters, total runtime, and summary counts (reads per biotype/feature, cluster counts, etc).|
|`<name>_5mers.txt`|Raw bound and background 5-mer counts.|
|`<name>_5mer_ratios.tsv`|Bound vs. background enrichment ratio per 5-mer, ranked from most to least enriched.|
|`<name>_5mer_enrichment.pdf`|Scatter plot illustrating 5-mer enrichment. Only generated if python3 with numpy and matplotlib are installed.|
|`<name>_ttoc_utr_cds.pdf`|Metagene plot of crosslinked read distribution across 5'UTR → CDS → 3'UTR. In CLIP mode, this plot shows the distribution of all aligned reads instead. Only generated if python3 with numpy and matplotlib are installed.|
|`<name>_ttoc_intron.pdf`|Metagene plot of crosslinked read distribution across introns (5' splice site → 3' splice site). In CLIP mode, this plot shows the distribution of all aligned reads instead. Only generated if python3 with numpy and matplotlib are installed.|
|`rRNA_aligned_<name>Aligned.sortedByCoord.out.bam`|Only if `-e` was used: reads that mapped to the rRNA reference. Retained for cases where rRNA binding is of biological interest.|

---

## 6\. Advanced options

### 6.1 rRNA pre-mapping (optional)

Build a small STAR index from the standard human/mouse rDNA repeat-unit
references, and pass it with `-e`:

```bash
wget -O human_rDNA_U13369.fa "https://www.ebi.ac.uk/ena/browser/api/fasta/U13369.1?download=true"
wget -O human_5S_X12811.fa "https://www.ebi.ac.uk/ena/browser/api/fasta/X12811.1?download=true"
cat human_rDNA_U13369.fa human_5S_X12811.fa > rRNA.fa
mkdir -p /path/to/rRNA_genome
STAR --runMode genomeGenerate \\
     --genomeDir /path/to/rRNA_genome \\
     --genomeFastaFiles rRNA.fa \\
     --genomeSAindexNbases 6 \\
     --runThreadN 8
```

`--genomeSAindexNbases 6` is required because the rRNA genome is very small (~43-45 kb) - refer to the STAR manual for details. For mouse, use accession `BK000964.3` and the same command.

---

## 7\. Troubleshooting

* **5-mer/metagene PDFs are missing** - these steps are skipped with a note in `_stats.txt` if `python3`/`matplotlib`/`numpy` aren't available, or if there is insufficient data for a given plot (e.g. too few clusters mapping to a particular region). Check `_stats.txt` for the specific reason.
* **`STAR --version` fails with `GLIBC_X.XX not found`** - a known issue
with bioconda's compiled STAR build on some systems. Install the
official precompiled binary from the STAR GitHub release page directly (for Linux systems) or compile from source (for Mac). For detailed instructions on appropriate STAR installation in the different environments see our extended installation guide (INSTALLATION.md).

**Error messages**

* **"illegal option" errors** - flags are case-sensitive: `-G`/`-g` and
`-C`/`-c` are different flags with different meanings.
* **"STAR genome index not found"** - `-x` must point at the index
*directory* (containing `SAindex`), not the genome FASTA.
* **"please provide a genome.fa using the -f option"** - `-f` is required
even when a STAR index already exists, since later steps
(`samtools calmd`, `bedtools getfasta`) need the FASTA directly.
* **"-A does not look like a proper annotation bed file"** - The `-A` flag expects a pre-built Trierarch BED file (`*_Paralus.bed`), not a raw GTF. To fix this, either run Trierarch standalone first with `--trierarch` and pass the resulting `*_Paralus.bed` to `-A`, or provide the raw GTF with `-G` and PARALUS will build the annotation BED before proceeding with the analysis.
* **"This looks like GFF3, not GTF"** - The `-G` flag expects a raw GTF file, not a GFF3. If needed, convert the GFF3 to GTF first.
* **"bedtools: command not found"** or similar - the `rnaseq` conda
environment isn't active in this terminal session; run
`conda activate rnaseq` again (it doesn't persist across terminal
windows/tabs).

* **"-c is set to 0 but -r is not. For CLIP mode (no T-to-C signal expected), both -c and -r must be set to 0. For PAR-CLIP mode, use the default values or set both above 0."** - `-c` and `-r` must either both be set to `0` (CLIP mode) or both be above `0` (PAR-CLIP mode). Setting only one of them to `0` is not allowed.
* **"-r is set to 0 but -c is not. For CLIP mode (no T-to-C signal expected), both -c and -r must be set to 0. For PAR-CLIP mode, use the default values or set both above 0."** - same as above, triggered when `-r` is `0` but `-c` is not.

