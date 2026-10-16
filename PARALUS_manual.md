# PARALUS Extended Manual

This is the full manual for PARALUS and its companion annotation builder, Trierarch (both in `paralus_1.8.7.sh`). For a guided first run, see `README.md` instead - this document is for looking things up, not for a first-time walkthrough.

---

## Table of contents

1. [Overview](#1-overview)
2. [Trierarch](#2-trierarch)
   - [2.1 Standalone run settings](#21-standalone-run-settings)
   - [2.2 Annotation construction and overlap resolution](#22-annotation-construction-and-overlap-resolution)
   - [2.3 Downloading external databases (optional)](#23-downloading-external-databases-optional)
3. [PARALUS](#3-paralus)
   - [3.1 Required input](#31-required-input)
   - [3.2 Preprocessing settings](#32-preprocessing-settings)
   - [3.3 Binding-site calling settings](#33-binding-site-calling-settings)
   - [3.4 Trierarch settings (used with `-G`)](#34-trierarch-settings-used-with--g)
4. [Pipeline walkthrough](#4-pipeline-walkthrough)
5. [Output files](#5-output-files)
6. [Advanced options](#6-advanced-options)
7. [Troubleshooting](#7-troubleshooting)

---

## 1\. Overview

PARALUS is a single Bash file that performs every step of PAR-CLIP/CLIP analysis: adapter/UMI trimming, optional rRNA depletion, genome alignment (STAR by default, or HISAT2 with `-H` - see Section 3.1 and 6.2), T-to-C crosslink-site extraction (for PAR-CLIP data), a two-stage group-to-cluster calling algorithm, genomic annotation of binding clusters, 5-mer sequence enrichment, and metagene analysis of crosslink site distribution across transcript regions - all from a single command, with no dependencies beyond cutadapt, an aligner (STAR, or HISAT2 for `-H`), samtools, bedtools, and (for the two plotting steps) python3 with numpy/matplotlib.

Trierarch, its companion annotation builder, is merged into the same file as a Bash function. It converts any GTF into the 7-column BED PARALUS needs, resolving overlapping annotations by an explicit priority order.

Both tools run with the same command, `paralus`. The installer (`install_paralus.sh`) creates a conda environment called `paralus` with all dependencies and the `paralus` command. Once the environment is activated (`conda activate paralus`), `paralus` works from any folder and runs the PARALUS script in your PARALUS folder (see [INSTALLATION.md](INSTALLATION.md)).

```bash
paralus --trierarch -G genes.gtf -C -O prefix   # Trierarch standalone run to pre-build the annotation BED
paralus -i ... -x ... -A ... [options]          # a standard PARALUS run for PAR-CLIP / CLIP analysis
```

---

## 2\. Trierarch

PARALUS requires a non-overlapping genomic annotation BED file to assign each binding cluster to a gene, feature, and biotype. Its built-in annotation builder, Trierarch, builds this file from a standard GTF annotation, with optional tRNA and miRNA database integration for more complete annotation coverage. The annotation BED only needs to be generated once per GTF (and set of Trierarch flags) and can be reused across all PARALUS runs that use that annotation. Trierarch can be run as a standalone step (using `--trierarch` mode), and the resulting BED file can be passed to PARALUS via `-A`. Alternatively, Trierarch can be triggered within a standard PARALUS run by providing a GTF file with the `-G` flag instead of a pre-built BED file. This is the simplest option for a first or single run, but each run with `-G` builds the annotation BED from the GTF again, even if it was built in an earlier run. To save time when analysing more samples from the same organism, build it once and reuse it with `-A`.

### 2.1 Standalone run settings

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
curl -O https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.bed

# miRNA (miRBase) - human example
curl -O https://www.mirbase.org/download/hsa.gff3
```

> **Note:** `curl` is preinstalled on macOS and most Linux systems. If yours doesn't have it, use `wget` with the same link instead.

For organisms outside either resource's coverage, Trierarch runs fine without them, relying on whatever small-RNA biotypes the GTF itself annotates.

---

## 3\. PARALUS

### 3.1 Required input

|Flag|Description|
|-|-|
|`-i input_file`|Raw reads file: `.fastq.gz`, `.fastq`, or `.fasta`, given as a file name or a full path. Output files are written to the input file's folder.|
|`-x genome_index`|Full path to the genome index: a STAR genomeDir by default, or a HISAT2 index basename (built with `hisat2-build`) if `-H` is set.|
|`-A annotation.bed` **or** `-G genes.gtf`|Either a pre-built Trierarch BED file (`-A`), or a raw GTF (`-G`) which triggers Trierarch to build the required BED file before the run starts.|
|`-f genome.fa`|Full path to the genome FASTA file. Required even if a genome index already exists, as it is needed for T-to-C crosslink site extraction (`samtools calmd`) and for retrieving the binding site sequences in the output (`bedtools getfasta`).|

### 3.2 Preprocessing settings

|Flag|Default|Description|
|-|-|-|
|`-n`|`19`|Minimum read length after trimming.|
|`-a`|none|3' adapter sequence for cutadapt.|
|`-d`|`GTTCAGAGTTCTACAGTCCGACGATC`|5' adapter sequence. If the 5' adapter contains a barcode, include the full barcoded sequence here.|
|`-u`|`0`|Bases trimmed from the 5' end after adapter removal (UMI length).|
|`-p`|`8`|Number of threads. The best value depends on the CPU cores and memory available on your machine.|
|`-e`|none|Full path to an rRNA reference genome index (STAR genomeDir, or HISAT2 basename if `-H` is set). When provided, reads are first mapped against this reference and those reads are removed before genome alignment.|
|`-H`|off (STAR is the default aligner)|Use HISAT2 instead of STAR for every alignment step. STAR needs its whole genome index in RAM to build and to align against (~30 GB for human/mouse), which most 16 GB Mac laptops can't do; HISAT2 uses a compressed genome index instead, needing about 8 GB of RAM to build the human hg38 index and about 5 GB to align against it. When set, `-x`/`-e` must be HISAT2 index basenames, not STAR genomeDir paths - see Section 6.2. STAR is recommended when memory allows, as it reports more clusters.|

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

1. **Preprocessing.** Identical reads are collapsed into one, removing PCR duplicates. The reads are then trimmed with cutadapt if a 3' adapter (`-a`) or UMI length (`-u`) is given, and rRNA-depleted if an rRNA index (`-e`) is given. Libraries without adapters or UMIs skip the trimming step.
2. **Alignment.** STAR (or HISAT2 with `-H`) is used to align reads end-to-end without allowing spliced alignments (`--alignIntronMax 1`/`--alignEndsType EndToEnd` for STAR, `--no-spliced-alignment`/`--no-softclip` for HISAT2) against the reference genome.
3. **Crosslink-site extraction.** T-to-C mismatches (the hallmark of PAR-CLIP, arising when 4SU-incorporated sites are UV-crosslinked and result in T-to-C substitutions during cDNA synthesis) are identified and extracted at single-nucleotide resolution. This step is skipped in CLIP mode (`-c 0 -r 0`).
4. **Group calling.** Overlapping reads are merged into groups. Each group's read distribution is tested for normality. Groups with an irregular read distribution (e.g. broad or multi-modal) are set aside as "NonNormal" and handled separately in the next step.
5. **NonNormal group re-evaluation.** NonNormal groups (with broad or irregular read distributions) are re-evaluated on a representative subsample of their reads. This helps in high-coverage regions where elevated background can obscure a genuine binding signal: subsampling reduces the noise and allows the underlying signal to emerge. Groups with sparse, scattered reads are typically discarded at this stage. The subsampling is deterministic, ensuring reproducible results across runs.
6. **Cluster promotion.** Groups that clear the `-c`/`-r` thresholds are promoted to clusters. In CLIP mode (`-c 0 -r 0`) all groups are promoted to clusters.
7. **Annotation.** Clusters are intersected with the annotation BED file. Each cluster is assigned a single gene name, feature, and biotype. For details on how overlapping annotations are resolved, see Section 2.2 (Trierarch).
8. **5-mer enrichment.** Counts all 1024 possible pentamers (5-mers) within the bound cluster sequences and compares them against the counts of 5-mers within the flanking genomic sequences upstream and downstream of each cluster (used as background). Raw counts are saved to `*_5mers.txt` and enrichment ratios to `*_5mer_ratios.tsv`.
9. **Metagene analysis.** Distribution of crosslinked reads (or all aligned reads in CLIP mode) is visualized across two metagene profiles: the full mature mRNA of protein-coding genes (5'UTR → CDS → 3'UTR as a continuous unit) and independently across introns. Each mRNA and intron is normalized by its total crosslinked reads count before averaging, so highly-bound loci do not dominate the profile.

---

## 5\. Output files

PARALUS produces the following output files. Temporary intermediate files created during the run are removed upon completion.

At the end of a run, PARALUS prints a short summary in the terminal: the number of binding sites (clusters), the run time and the folder with the results. If no clusters were found, it says so and points to the run log.

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

Build a small STAR index from the standard human/mouse rDNA repeat-unit references, and pass it with `-e`:

```bash
curl -o human_rDNA_U13369.fa "https://www.ebi.ac.uk/ena/browser/api/fasta/U13369.1?download=true"
curl -o human_5S_X12811.fa "https://www.ebi.ac.uk/ena/browser/api/fasta/X12811.1?download=true"
cat human_rDNA_U13369.fa human_5S_X12811.fa > rRNA.fa
mkdir -p /path/to/rRNA_genome
STAR --runMode genomeGenerate \
     --genomeDir /path/to/rRNA_genome \
     --genomeFastaFiles rRNA.fa \
     --genomeSAindexNbases 6 \
     --runThreadN 8
```

> **Note:** If you don't have `curl`, use `wget -O human_rDNA_U13369.fa "<link>"` instead (capital `-O` in wget).

`--genomeSAindexNbases 6` is required because the rRNA genome is very small (~43-45 kb) - refer to the STAR manual for details. For mouse, use accession `BK000964.3` and the same command.

### 6.2 HISAT2 mode (`-H`, for machines with limited RAM)

STAR needs its whole genome index in RAM, both to build it and to align against it (~30 GB for human/mouse). On a Mac laptop with 16 GB of RAM or less, that usually isn't practical. HISAT2 is offered as a second aligner for exactly this case: it uses a compressed genome index instead, needing about 8 GB of RAM to build the human hg38 index and about 5 GB to align against it.

Pass `-H` to switch every alignment step in PARALUS - the main genome alignment and, if used, the rRNA pre-mapping - from STAR to HISAT2. HISAT2 recovers the same biological signal (on CPSF1, both aligners showed the 3'UTR binding spike and enrichment of the AATAAA polyadenylation signal), but reports fewer clusters (845 vs. 1,241 with STAR), so STAR is recommended when memory allows. Internally, `--no-spliced-alignment --no-softclip` are used to match STAR's `--alignIntronMax 1 --alignEndsType EndToEnd` (no spliced alignment, no soft-clipping, on either aligner).

The only practical difference is what `-x` and `-e` expect: with `-H`, they take a HISAT2 index **basename** (built with `hisat2-build`, not a directory) instead of a STAR genomeDir:

```bash
mkdir -p /path/to/genome_index
hisat2-build /path/to/hg38.fa /path/to/genome_index/genome
# pass -x /path/to/genome_index/genome -H
```

The optional rRNA index (`-e`) is built the same way:

```bash
mkdir -p /path/to/rRNA_index
hisat2-build rRNA.fa /path/to/rRNA_index/rRNA
# pass -e /path/to/rRNA_index/rRNA -H
```

See [INSTALLATION.md](INSTALLATION.md) for installing HISAT2 via conda.

---

## 7\. Troubleshooting

* **5-mer/metagene PDFs are missing** - these steps are skipped with a note in `_stats.txt` if `python3`/`matplotlib`/`numpy` aren't available, or if there is insufficient data for a given plot (e.g. too few clusters mapping to a particular region). Check `_stats.txt` for the specific reason.
* **`STAR --version` fails with `GLIBC_X.XX not found`** - a known issue with bioconda's compiled STAR build on some systems. Install the official precompiled binary from the STAR GitHub release page directly (for Linux systems) or compile from source (for Mac). For detailed instructions on appropriate STAR installation in the different environments see our extended installation guide (INSTALLATION.md). The installer (`install_paralus.sh`) handles this for you: it installs the official STAR binary on Linux and compiles STAR on a Mac. If STAR remains impractical on your machine (e.g. limited RAM), consider HISAT2 mode (`-H`, Section 6.2) instead.
* **STAR runs without errors but aligns 0 reads (Apple Silicon Macs)** - this is a known problem with the STAR build from conda (bioconda) on Apple Silicon. Compile STAR from source instead (see [INSTALLATION.md](INSTALLATION.md)), or use the installer (`install_paralus.sh`), which compiles it for you.

**Error messages**

* **"paralus: command not found"** - the `paralus` command isn't available in this terminal session. If you installed with conda (the default installer setting), activate the environment with `conda activate paralus`. If you ran the installer with `--launcher-only`, run `source ~/.bashrc` (`source ~/.zshrc` on a Mac). If you have moved the PARALUS folder, run `bash install_paralus.sh` again from inside the PARALUS folder.
* **"error: input folder ... not found"** - the folder in the `-i` path doesn't exist, so PARALUS can't find the input file. Check the path for typos.
* **"PARALUS finished, but no binding sites (clusters) were found"** - check the run log named in the message for errors. If the run itself completed, look at the group count in `_stats.txt`: if groups were found but none became clusters, the data may have little T-to-C signal (for CLIP data, use `-c 0 -r 0`). If few reads aligned, check the adapter sequence (`-a`) and the genome index (`-x`).
* **"Could not find any combinations of options for the input"** - this combination of input file and options isn't supported. Check that the file ends in `.fastq.gz`, `.fastq` or `.fasta`, that `-u` (UMI length) is only used together with `-a` (3' adapter), and that `.fasta` input (already trimmed) is run without `-a` or `-u`.
* **"illegal option" errors** - flags are case-sensitive: `-G`/`-g` and `-C`/`-c` are different flags with different meanings.
* **"STAR genome index not found"** / **"HISAT2 genome index not found"** - `-x` must point at the index *directory* (containing `SAindex`) for STAR, or at the index *basename* (e.g. `/path/genome_index/genome`, matching files like `genome.1.ht2`) for HISAT2 (`-H`) - not the genome FASTA either way.
* **"please provide a genome.fa using the -f option"** - `-f` is required even when a genome index already exists, since later steps need the FASTA directly.
* **"-A does not look like a proper annotation bed file"** - The `-A` flag expects a pre-built Trierarch BED file (`*_Paralus.bed`), not a raw GTF. To fix this, either run Trierarch standalone first with `--trierarch` and pass the resulting `*_Paralus.bed` to `-A`, or provide the raw GTF with `-G` and PARALUS will build the annotation BED before proceeding with the analysis.
* **"This looks like GFF3, not GTF"** - The `-G` flag expects a raw GTF file, not a GFF3. If needed, convert the GFF3 to GTF first.
* **"error: missing dependencies: ..."** - one or more required tools aren't available in this terminal session. If you used the installer (`install_paralus.sh`), activate the environment with `conda activate paralus` (it doesn't persist across terminal windows/tabs). Otherwise, install the listed tools or load their modules.
* **"-c is set to 0 but -r is not. For CLIP mode (no T-to-C signal expected), both -c and -r must be set to 0. For PAR-CLIP mode, use the default values or set both above 0."** - `-c` and `-r` must either both be set to `0` (CLIP mode) or both be above `0` (PAR-CLIP mode). Setting only one of them to `0` is not allowed.
* **"-r is set to 0 but -c is not. For CLIP mode (no T-to-C signal expected), both -c and -r must be set to 0. For PAR-CLIP mode, use the default values or set both above 0."** - same as above, triggered when `-r` is `0` but `-c` is not.

