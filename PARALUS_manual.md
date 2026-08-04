# PARALUS Extended Manual

This is the full reference manual for PARALUS and its companion annotation
builder, Trierarch (`paralus\_1.7.0.sh`). For a guided first run, see
`PARALUS\_quickstart.md` instead - this document is for looking things up,
not for a first-time walkthrough.

\---

## Table of contents

1. [Overview](#1-overview)
2. [Trierarch reference](#2-trierarch-reference)
3. [PARALUS reference](#3-paralus-reference)
4. [Pipeline walkthrough](#4-pipeline-walkthrough)
5. [Output files](#5-output-files)
6. [Advanced recipes](#6-advanced-recipes)
7. [Troubleshooting](#7-troubleshooting)

\---

## 1\. Overview

PARALUS is a single Bash file that performs every step of PAR-CLIP/CLIP
analysis: adapter/UMI trimming, optional rRNA depletion, genome alignment
(STAR), T-to-C (or G-to-A) crosslink-site extraction, a two-stage
group-to-cluster calling algorithm, annotation against a gene model,
5-mer sequence enrichment, and transcript-region metagene profiling -
all from one command, with no dependencies beyond cutadapt, STAR,
samtools, bedtools, and (for the two plotting steps) python3 with
numpy/matplotlib.

Trierarch, its companion annotation builder, is merged into the same
file as a Bash function. It converts any GTF into the 7-column BED
PARALUS needs, resolving overlapping small-RNA/protein-coding
annotations by an explicit priority order rather than ad-hoc patching.

Both tools are invoked from the same script:

```bash
./paralus\_1.7.0.sh -i ... -x ... -A ... \[options]         # a normal PARALUS run
./paralus\_1.7.0.sh --trierarch -G genes.gtf -C -O prefix   # Trierarch standalone
```

\---

## 2\. Trierarch reference

### 2.1 Flags

|Flag|Required|Description|
|-|-|-|
|`-G genes.gtf`|yes|GTF annotation (GENCODE or Ensembl format). Recognises either the `gene\_type` or `gene\_biotype` attribute, whichever the file uses.|
|`-T tRNA\_DB.bed`|no|tRNA coordinates, bed6 (chr, start, end, name, score, strand) - e.g. from GtRNAdb.|
|`-M miRNA\_DB`|no|miRNA coordinates - either a miRBase GFF3 (auto-detected by column count) or a plain bed6 in the same shape as `-T`.|
|`-O output\_prefix`|no|Output file prefix. Default: the GTF's own filename with `.gtf` stripped. A bare name is placed next to the GTF; a path (containing `/`) is used as given.|
|`-C`|no|Restrict extraction to `Ensembl\_canonical`/`MANE\_Select`-tagged transcripts only. **Off by default.** Strongly recommended when you plan to use the metagene tools (Section 6.2), since without it a gene's 5'UTR/CDS/3'UTR boundaries are ambiguous (they depend on which isoform).|
|`-U`|no|Keep genes with no `gene\_name` attribute, falling back to `gene\_id` (e.g. an ENSG... accession) as the display name. **Off by default** - these rows are dropped instead, since an accession-only name usually indicates a novel/unannotated locus.|

### 2.2 What it builds, and in what order

Trierarch resolves overlapping annotations using a fixed priority cascade,
highest first: `snoRNA > snRNA > rRNA > miRNA > tRNA > Mt\_tRNA > scaRNA > lncRNA > protein\_coding > everything else`. Each class's bases are
subtracted from every lower-priority class before that class is added to
the final BED, so the output is conflict-free by construction - no
downstream patching is needed.

Within `protein\_coding`, 5'UTR/3'UTR are derived from the union of a
gene's exonic bases *minus* the union of its CDS bases (not the widest
start-codon-to-stop-codon span across all isoforms, which over-calls UTR
on multi-isoform genes). A remaining piece that falls inside a gene's
overall CDS envelope but isn't itself CDS (e.g. a cassette exon that's
coding in one isoform and skipped in another) is labelled `exon` rather
than guessed at.

### 2.3 GTF vs GFF3

Trierarch's attribute parsing assumes GTF's `key "value";` syntax. Given a
real GFF3 (`key=value;`), it fails loudly with a clear error and a
suggested `gffread`/AGAT conversion command, rather than silently
extracting nothing.

### 2.4 Downloading external databases

```bash
# tRNA (GtRNAdb) - human hg38 example
wget https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.bed

# miRNA (miRBase) - human example
wget https://www.mirbase.org/download/hsa.gff3
```

Both are eukaryote-specific resources (bacteria have no tRNA-database
equivalent conventionally used here, and no miRNA pathway at all) - for
organisms outside either resource's coverage, Trierarch runs fine without
them, relying on whatever small-RNA biotypes the GTF itself annotates.

\---

## 3\. PARALUS reference

### 3.1 Required flags

|Flag|Description|
|-|-|
|`-i input\_file`|Reads: `.fastq.gz`, `.fastq`, or `.fasta`.|
|`-x genome\_index`|STAR genome index directory.|
|`-A annotation.bed` **or** `-G genes.gtf`|Either a ready-made Trierarch bed (`-A`), or a raw GTF (`-G`), which triggers Trierarch automatically the first time and reuses the cached result on later runs with the same GTF/flags.|
|`-f genome.fa`|Genome FASTA - required even if a STAR index already exists (needed later for `samtools calmd` and `bedtools getfasta`).|

### 3.2 Optional preprocessing flags

|Flag|Default|Description|
|-|-|-|
|`-n`|`19`|Minimum read length after trimming.|
|`-m`|`3`|Minimum reads to define a group.|
|`-a`|none|3' adapter sequence for cutadapt.|
|`-d`|`GTTCAGAGTTCTACAGTCCGACGATC`|5' adapter sequence.|
|`-u`|`0`|Bases trimmed from the 5' end after adapter removal (UMI length). Omit if the 5' adapter itself is barcoded.|
|`-p`|`8`|Threads. Don't exceed \~12 - STAR itself recommends against more.|
|`-e`|none|STAR index for an rRNA reference, to pre-map and report/deplete rRNA reads before genome alignment (Section 6.1).|

### 3.3 Binding-site calling flags

|Flag|Default|Description|
|-|-|-|
|`-c`|`1`|Minimum T-to-C conversion locations to consider a group a cluster.|
|`-r`|`0.3`|Minimum fraction of reads at a group that must carry a T-to-C conversion for it to be called a cluster.|

### 3.4 Annotation-passthrough flags (used with `-G`)

Same letters Trierarch uses natively (Section 2.1) - passed straight
through when `-G` triggers it, no translation between standalone and
`-G` usage: `-T`, `-M`, `-C`, `-U`.

\---

## 4\. Pipeline walkthrough

1. **Preprocessing.** Reads are collapsed by unique sequence (with counts
preserved in the fasta headers), optionally rRNA-pre-mapped (`-e`),
then adapter/UMI-trimmed with cutadapt.
2. **Alignment.** STAR aligns in ungapped, end-to-end mode
(`--alignIntronMax 1`) - splicing is not modelled at the alignment
stage; intron/exon context is instead recovered later from the
annotation.
3. **Crosslink-site extraction.** `samtools calmd` regenerates the MD tag,
then a strand-aware awk pass identifies T-to-C (forward strand,
SAM flag 0/256) or G-to-A (reverse strand, SAM flag 16/272)
mismatches consistent with 4SU/6SG-induced conversions, handling
both single- and double-conversion reads.
4. **Group calling.** Overlapping reads are merged into groups; each
group's read-start/end spread is compared against a
standard-deviation-based normality threshold (68%, i.e. roughly
1 SD either side, the expected coverage for a normal distribution).
Groups that fail this test (satellite/repetitive loci, broad or
multi-modal read stacks) are set aside as "NonNormal" rather than
silently kept or dropped.
5. **NonNormal subsampling.** NonNormal groups are re-evaluated on a
subsample of their reads, drawn via a deterministic hash of each read
name (not `rand()`, which was found to be non-reproducible in this
environment across repeated runs) - verified reproducible and
well-distributed even when read names are sequential.
6. **Cluster promotion.** Groups (from both the normal and re-evaluated
NonNormal passes) become clusters once they clear the `-c`/`-r`
thresholds.
7. **Annotation.** Clusters are intersected against the annotation bed.
Each candidate hit's gene/feature/biotype is packed into one unit
before merging, so a cluster's reported annotation always describes
one real, internally-consistent feature - ranked by the same fixed
biotype/feature priority Trierarch itself uses, with overlap length
as a tiebreaker within a tied priority.
8. **5-mer enrichment.** Both an upstream and downstream real flanking
region are taken per cluster (each independently dropped if it would
cross a chromosome boundary or land on another real cluster); 5-mer
counts from both sides are averaged and compared against bound-site
counts, ranked by enrichment ratio with a pseudocount.
9. **Metagene analysis.** Crosslink locations are mapped onto a
transcript-ordered concatenation of each gene's 5'UTR/CDS/3'UTR
segments (from the canonical annotation), and independently onto each
intron's own splice-site-to-splice-site span, weighted by real
per-location read depth, normalized per gene/intron before averaging
so no single locus dominates.

\---

## 5\. Output files

Only files that survive to the end of a run are listed - everything else
is an intermediate PARALUS creates and deletes automatically.

|File|Contents|
|-|-|
|`<name>.Aligned.sortedByCoord.out.bam` (+`.bai`)|Full genome alignment.|
|`<name>.aligned\_TtoC.bam` (+`.bai`)|Reads carrying a real T-to-C/G-to-A conversion.|
|`<name>.aligned\_TtoC.bed`|Per-read bed of conversion-supporting reads (kept specifically so depth-weighted downstream analyses, e.g. the metagene tools, can be built from it later).|
|`<name>.groups.csv`|Every candidate group, before cluster-threshold filtering - useful for CLIP data with no conversion signal to filter on (Section 6.3).|
|`<name>.clusters.csv`|Final annotated binding sites: gene, feature, biotype, sequence, read/conversion support, per row.|
|`<name>.clusters.bed`|Bed version of the same clusters.|
|`<name>\_run\_<timestamp>.log`|Full stdout/stderr of the run, plus the exact invoked command.|
|`<name>\_stats.txt`|Run parameters, per-stage timing, and summary statistics (reads per biotype/feature, cluster counts, etc).|
|`<name>\_5mers.txt`|Raw background vs bound 5-mer counts.|
|`<name>\_5mer\_ratios.tsv`|Same data, ranked by enrichment ratio.|
|`<name>\_5mer\_enrichment.pdf`|Labelled scatter of the above (top/bottom enriched motifs).|
|`<name>\_ttoc\_utr\_cds.pdf`|5'UTR→CDS→3'UTR metagene, depth-weighted.|
|`<name>\_ttoc\_intron.pdf`|Intron metagene (5' splice site → 3' splice site), depth-weighted.|
|`rRNA\_aligned\_<name>Aligned.sortedByCoord.out.bam`|Only if `-e` was used: reads that mapped to the rRNA reference - kept deliberately, since binding on rRNA itself can be real biological signal, not just background to discard.|

\---

## 6\. Advanced recipes

### 6.1 rRNA pre-mapping

Build a small STAR index from the standard human/mouse rDNA repeat-unit
references, and pass it with `-e`:

```bash
wget -O human\_rDNA\_U13369.fa "https://www.ebi.ac.uk/ena/browser/api/fasta/U13369.1?download=true"
mkdir -p /path/to/rRNA\_genome
STAR --runMode genomeGenerate \\
     --genomeDir /path/to/rRNA\_genome \\
     --genomeFastaFiles human\_rDNA\_U13369.fa \\
     --genomeSAindexNbases 6 \\
     --runThreadN 8
```

`--genomeSAindexNbases 6` matters here specifically - these references are
tiny (43-45kb) compared to a real genome, and STAR's own guidance
(`min(14, log2(GenomeLength)/2 - 1)`) requires a much smaller value than
the default for genomes this size. Mouse equivalent: accession
`BK000964.3`, same command otherwise.

### 6.2 Canonical annotation for the metagene tools

```bash
./paralus\_1.7.0.sh --trierarch -G genes.gtf -C -O canonical\_prefix
```

Use the resulting `canonical\_prefix\_Paralus.bed` with `-A` in your real
run. Without `-C`, a gene can have multiple isoforms with different
UTR/CDS boundaries, so "this gene's 5'UTR" isn't a single well-defined
thing, and fewer genes will have all three regions (5'UTR, CDS, 3'UTR)
cleanly resolvable.

### 6.3 Plain CLIP (no 4SU/6SG labelling)

`clusters.csv` requires T-to-C support to promote a group into a cluster

* for CLIP data with no metabolic labelling and thus no conversion signal
at all, this filter can discard nearly everything real. Use
`groups.csv` instead for downstream analysis (e.g. building a metagene
input bed), since it has the same column layout without the
conversion-based filter applied.

### 6.4 Weighting downstream analyses by real read depth

```bash
# clusters, weighted by read count
awk -F "," 'NR>1{OFS="\\t"; print $1,$2,$3,$4,$11,$6}' name.clusters.csv > clusters\_weighted.bed

# groups, weighted by read count (for plain CLIP data)
awk -F "," 'NR>1{OFS="\\t"; print $1,$2,$3,$4,$12,$6}' name.groups.csv > groups\_weighted.bed
```

Column 5 (score) of a bed file feeding any tool that supports weighting
is read as a per-site confidence weight if numeric, or defaults to 1 if
`.` - so an unweighted bed (score column `.`) and a weighted one (real
counts) both work with the same downstream tooling.

\---

## 7\. Troubleshooting

* **"illegal option" errors** - flags are case-sensitive: `-G`/`-g` and
`-C`/`-c` are different flags with different meanings.
* **"STAR genome index not found"** - `-x` must point at the index
*directory* (containing `SAindex`), not the genome FASTA.
* **"please provide a genome.fa using the -f option"** - `-f` is required
even when a STAR index already exists, since later steps
(`samtools calmd`, `bedtools getfasta`) need the FASTA directly.
* **"-A does not look like a proper annotation bed file"** - a raw GTF
was passed to `-A`. Either use `-G` instead (Trierarch runs
automatically), or build the bed with `--trierarch` first and pass the
resulting `\*\_Paralus.bed` to `-A`.
* **This looks like GFF3, not GTF** - convert first:
`gffread input.gff3 -T -o output.gtf` (or AGAT's
`agat\_convert\_sp\_gff2gtf.pl` for messier/non-standard GFF3).
* **`STAR --version` fails with `GLIBC\_X.XX not found`** - a known issue
with bioconda's compiled STAR build on some systems. Install the
official precompiled binary directly instead:

```bash
  wget https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR\_2.7.11b.zip
  unzip STAR\_2.7.11b.zip
  sudo cp STAR\_2.7.11b/Linux\_x86\_64\_static/STAR /usr/local/bin/
  ```

* **"bedtools: command not found"** or similar - the `paralus` conda
environment isn't active in this terminal session; run
`conda activate paralus` again (it doesn't persist across terminal
windows/tabs).
* **5-mer/metagene PDFs are missing** - these steps are skipped
(silently, with a note in `\_stats.txt`) if `python3`/`matplotlib`
aren't available, or if a genuine zero-qualifying-data case is hit
(e.g. too few sites landing in a given region). Check `\_stats.txt`
for the specific reason.

