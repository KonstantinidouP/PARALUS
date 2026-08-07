#!/bin/bash
# PARAlus is a bedtools and awk based program to analyze PAR-CLIP or CLIP data
# created by Dimitrios Anastasakis
# PARalus 1.5.1 this version was with addition of --outSAMmultNmax 1 --outMultimapperOrder Random 
# PARalus 1.5.3 this version was with addition of rRNA mapper toget rRNA maping stats and also with a cluster.csv cleaner with priority in annotations as 3ss > 5ss > start_codon > stop_codon> 3utr> 5utr > CDS > exon > intron 
# PARalus 1.6:
#   - fixed a real strand bug in TtoC crosslink-site extraction: the two
#     single-mismatch reverse-strand (G->A) branches were gated on the
#     forward-strand SAM flag (0/256) instead of reverse (16/272), silently
#     dropping single-conversion reverse-strand reads entirely (no error,
#     just missing signal - hits minus-strand genes hardest)
#   - fixed cluster annotation resolution: gene/annotation/biotype used to
#     be resolved as three independent "|"-joined lists with no positional
#     correspondence, which could report a gene+annotation+biotype
#     combination that matched no real feature (e.g. GENE_A's annotation
#     stitched to GENE_B's biotype). Now each hit's gene:annotation:biotype
#     is packed as one unit that survives the merge intact, and a single
#     pass picks one whole winning hit
#   - added overlap-bp as a tiebreaker when two hits are tied on the same
#     biotype/annotation priority tier (previously arbitrary/first-in-list)
#   - the six trailing sed patches (e.g. 3utr,tRNA -> exon,tRNA) are removed:
#     they were compensating for exactly the annotation-bed overlap bug this
#     fix addresses at the source, and the new packed-unit format doesn't
#     match their old text pattern anyway
#   - -f (genome fasta) is now validated up front and independently of
#     whether the STAR index already exists, instead of only erroring when
#     BOTH are missing (previously could run for 10+ minutes of STAR
#     mapping before failing much later at samtools calmd with an empty
#     reference path)
#   - added an explicit STAR index check as its own clear error, split out
#     from the -f check above
#   - full run logging: stdout+stderr now tee'd to a timestamped log file,
#     and the exact invoked command is recorded in both that log and
#     _stats.txt
#   - cleanup rm calls now use -f (several would throw a cosmetic "cannot
#     remove: no such file" if a given run's code path never created that
#     intermediate file); removed one fully redundant duplicate rm at the
#     very end of the script
#   - removed dead variable references ($gtf, $tmapa, $genome - never
#     assigned by any option, leftover from an older version) from the
#     stats output; replaced with the annotation file path instead
#   - fixed lincRNA/lncRNA typo (the pipeline only ever produces "lncRNA" as
#     a biotype value, so the "lincRNA_XL_reads" stat was permanently zero)
#   - fixed swapped snoRNA_clusters/snRNA_clusters stat labels (each was
#     being printed under the other's name)
# PARalus 1.6.1:
#   - -A now requires Trierarch's bed output, not a raw gtf: a format guard
#     checks column count and numeric start/end before anything runs, and
#     fails immediately with a clear message instead of burning 10+ minutes
#     of STAR/samtools before crashing deep in a bedtools step
#   - alt/unplaced scaffold contigs (e.g. GL000008.2-style GenBank
#     accessions) are filtered out of every downstream bed file via a
#     chromosome whitelist built from whatever contigs the annotation bed
#     actually uses - organism/naming-convention agnostic, since it's driven
#     by the annotation rather than a hardcoded "^chr" pattern. Note this
#     will NOT remove a scaffold that genuinely carries an annotated gene
#     (e.g. GENCODE does annotate some genes on GL000008.2) - that's real
#     biology, not naming noise
#   - 5mer enrichment background redesigned: same-gene-constrained (nearest
#     available proxy to "same transcript" - the annotation is gene-level),
#     strand-aware (shift direction now follows transcript orientation, not
#     just genomic +/-), clamped to the gene's own boundaries, and excluded
#     entirely (from both bound and background, keeping them paired 1:1) if
#     the candidate background would overlap any real cluster - including
#     its own, if clamping pushed it that close. This also fixes the old
#     background/foreground sample-size mismatch near chromosome starts, and
#     the k-mer counting implementation is now a single fast pass instead of
#     1024 separate grep calls (counting convention unchanged - intentionally
#     non-overlapping, matching the original grep -o behaviour)
#   - -G added: point at a raw gtf and Trierarch runs automatically to build
#     the annotation bed (~1 min), cached by gtf/database modification time
#     on later runs. -t/-M/-C/-U forward tRNA_DB/miRNA_DB/canonical-only/
#     keep-unnamed to Trierarch (uppercase M/C/U since lowercase m/c/u are
#     already Paralus's own options). -A still works unchanged if you
#     already have a bed built
# PARalus 1.6.2:
#   - Trierarch merged into this same file as a function (run_trierarch),
#     so the whole pipeline is one script for GitHub instead of two. Every
#     variable it uses is declared local, and all of its own `exit` calls
#     became `return`, so it can never leak state into or be killed
#     independently of the rest of this script. Two ways to reach it:
#     standalone (./paralus.sh --trierarch -G genes.gtf -C -o output) or
#     automatically from a normal run via -G, same as before - -G no longer
#     needs to locate a separate trierarch.sh file at all
# PARalus 1.6.3:
#   - -M/-C/-U are now Trierarch's own native flags (its getopts previously
#     still used lowercase -m/-c/-u internally, requiring Paralus's -g
#     integration to translate case on every call - now identical letters
#     work standalone and via -g, no translation, one less thing to get
#     out of sync)
#   - -g renamed to -G throughout (both Paralus's own flag and Trierarch's
#     native one): pairs with -A (annotation bed) the same way -G pairs
#     with it conceptually - both uppercase, both "how the annotation is
#     obtained", vs the mostly-lowercase preprocessing flags below
# PARalus 1.6.4:
#   - -t/-o renamed to -T/-O (Trierarch's tRNA_DB and output-prefix flags):
#     no collision ever forced these two to be lowercase, unlike M/C/U/G,
#     but leaving them lowercase broke the "Trierarch's flags are all
#     uppercase, standalone or via Paralus" rule for no real reason
# PARalus 1.6.5:
#   - k-mer enrichment generalized from a hardcoded 5-mer to 3-mer, 4-mer,
#     and 5-mer, run automatically (perl's mer-list generator now builds any
#     k iteratively instead of 5 hardcoded nested loops)
#   - background redesigned to use BOTH an upstream and a downstream real
#     flanking region per cluster (strand-aware, same-gene-clamped, as
#     before), instead of one side only. A cluster is kept only if both
#     sides come back clean (full length, no overlap with any real cluster);
#     otherwise it's dropped from bound/upstream/downstream alike, keeping
#     all three paired 1:1:1. The two background counts are averaged per
#     k-mer before comparing against bound
#   - added a log-log scatter plot per k (avg background on x, bound on y)
#     via python3/matplotlib, written as <name>_kmer<k>_enrichment.pdf; the underlying
#     counts are also kept as <name>_kmer<k>.tsv. Falls back to a note in
#     _stats.txt (skips the plot, keeps the tsv) if python3/matplotlib isn't available
#     rather than failing the whole run over a missing plotting tool
# PARalus 1.6.6:
#   - fixed a real, pre-existing bug that made k-mer enrichment silently
#     produce empty/all-zero output on real data: "${name}".clusters.bed has
#     always been written with columns 4 and 5 swapped relative to normal
#     BED convention (column 4 is a literal "." placeholder, column 5 is the
#     actual cluster ID - inherited from how the CSV-to-bed conversion was
#     already written, predating the k-mer feature entirely). The k-mer
#     code assumed the standard convention (cluster ID in column 4), so
#     every cluster's ID looked identical ("."), collapsing the whole
#     genome's clusters down to a single row during deduplication. Fixed by
#     reading column 5 where the real cluster ID actually lives. My own
#     synthetic tests never caught this because they used hand-built
#     fixtures with the ID already (incorrectly) in column 4, matching my
#     wrong assumption rather than the real file
#   - matplotlib now degrades gracefully (skips the plot, notes it in
#     _stats.txt, keeps the tsv) if a k-mer's data ends up empty/all-zero,
#     instead of crashing with a traceback
# PARalus 1.6.7:
#   - the 1.6.6 fix was verified correct on re-inspection, but the same
#     "Unable to open file ..._candidate_upstream.bed" error recurred on a
#     real run, meaning something else is going on that couldn't be seen
#     from the error alone - the diagnostic file (_cluster_gene.tsv) was
#     being deleted immediately after use, so there was no way to tell where
#     the data actually disappeared. Added real row-count diagnostics to
#     _stats.txt at each stage (clusters, gene envelopes, clusters matched to
#     a gene, upstream/downstream candidates that fit within their gene
#     bounds) so the next run points at the actual failure point instead of
#     another guess
#   - both candidate bed files are now touched into existence even if empty,
#     so a genuine zero-candidate case degrades gracefully (empty k-mer
#     data, plot skipped, already handled by the 1.6.6 guard) instead of
#     bedtools hard-crashing with "Unable to open file"
# PARalus 1.6.8:
#   - reverted the k-mer section entirely back to the original, pre-1.6.5
#     logic (single flanking-shift background, unpaired, 1024 separate grep
#     calls, 5-mer only) - the 1.6.5-1.6.7 redesign (dual-sided background,
#     3/4/5-mer, gene-envelope pairing) kept surfacing real but hard-to-
#     reproduce bugs on real genome-scale data, so going back to the known-
#     stable version rather than continuing to debug blind
#   - added a small python3/matplotlib addition on top of the ORIGINAL
#     logic only: a log-log scatter plot of the existing background/bound
#     counts already written to <name>_5mers.txt, saved as
#     <name>_5mer_enrichment.pdf. Nothing about the background/counting
#     itself changed - same file, same columns, same numbers as before;
#     this just reads that output and plots it. Skipped silently if
#     matplotlib isn't available (core pipeline is unaffected either way)
# PARalus 1.6.9:
#   - replaced the inline plotting snippet with a small standalone-quality
#     tool (kmer_enrichment.py), embedded as a heredoc and written out at
#     runtime so it's still one file for GitHub. Adds: ratio-based ranking
#     (with a pseudocount, since real k-mer data is sparse enough that raw
#     count2/count1 constantly hits divide-by-zero/infinite ratios that
#     would otherwise dominate the "most enriched" list from pure sampling
#     noise), a top/bottom-10 enriched-motif summary appended to
#     _stats.txt, and a labeled PNG scatter (red=enriched, blue=de-enriched)
#     instead of the earlier unlabeled log-log PDF. Output is now
#     <name>_5mer_ratios.tsv (ranked) + <name>_5mer_enrichment.png

# PARalus 1.6.18 (branched from the known-good 1.6.9 baseline, NOT the same
# lineage as the abandoned 1.6.10-1.6.14 - those had a real reproducibility
# bug in the subsampling rewrite, not caught until a user compared repeated
# runs; fixed here from scratch, one loop at a time, each verified on real
# data before moving to the next):
#   - NonNormal-cluster subsampling: used to spawn one separate samtools
#     process per region (`samtools view -s fraction region`), replaced
#     with one pass - bedtools intersect pairs every read with its region's
#     fraction, one awk pass makes the keep/drop draw per read (a
#     deterministic hash of the read name, not awk's rand() - confirmed
#     non-reproducible across process invocations in this environment even
#     with an explicit seed), one final samtools call extracts the kept
#     reads. Verified reproducible via md5sum on two independent runs of
#     the exact assembled block, and the hash itself verified
#     well-distributed (5000 sequential test names landed 490-519 per
#     decile, expected ~500 - important since read names here are often
#     sequential, and a weak hash would bias the subsample toward
#     contiguous blocks rather than a genuine scatter)
#   - 5-mer counting: used to spawn a separate grep process per mer per
#     fasta file (2048 subprocess spawns total), replaced with one awk pass
#     per fasta file. An EARLIER attempt at this same fix (in the abandoned
#     1.6.10-1.6.14 lineage) was shipped as "verified byte-identical,"
#     which was wrong - that test only used 3-4 hand-picked mers and never
#     exposed either of two real bugs, both only caught because a user
#     insisted the counts really were wrong on real data rather than
#     accepting "already tested":
#       1) a single shared "skip 5 ahead after a match" pointer is not the
#          same as giving each of the 1024 mers its own independent grep
#          scan - two different mers can legitimately match at overlapping
#          starting positions (e.g. "AAAATTT" contains AAAAT, AAATT, and
#          AATTT, at positions 1/2/3 - grep finds all three independently;
#          a shared pointer finds only the first and skips past the other
#          two). Fixed by tracking a separate "next allowed position" per
#          mer instead of one shared pointer
#       2) that per-mer tracking has to reset at the start of every
#          sequence (fasta line/record) - grep's matching is inherently
#          line-scoped, so a mer's match near the end of one sequence must
#          not block a legitimate match near the start of an unrelated one
#     This time verified byte-identical against the real grep -o loop at
#     full 1024-mer scale, on a multi-sequence test file, not a handful of
#     hand-picked mers - the exact conditions the previous "verified" claim
#     was missing
#   - removed the last remaining perl dependency: this branch started from
#     the pre-perl-removal 1.6.9 baseline, so it never got the fix applied
#     earlier in the abandoned 1.6.10 lineage. The 1024-5mer list (5 nested
#     perl loops over A/T/C/G) is replaced with bash's own brace expansion
#     (`echo {A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}`), verified to
#     produce the identical set of 1024 5-mers - one less tool to install,
#     since bash is already a hard requirement anyway
#   - added timing lines to _stats.txt: total pipeline wall-clock time, plus
#     per-stage time for the NonNormal subsampling step, so the actual time
#     saved on real data is measurable rather than assumed

# PARalus 1.6.19:
#   - 5-mer background is now both an upstream (-200) and downstream (+200)
#     flanking region per cluster, not just one side - kept deliberately
#     simple (flat genomic shift, not strand-aware/gene-constrained like an
#     earlier, more complex attempt at this that kept surfacing
#     hard-to-reproduce bugs on real data and was reverted). Each side is
#     independently dropped if the shift would go past either end of the
#     chromosome (checked against the genome fasta's .fai index) or if it
#     lands on any real cluster (bedtools subtract -A against
#     clusters.bed) - a cluster can lose one side and keep the other. The
#     final background count per k-mer is the average of its upstream and
#     downstream counts, rounded to the nearest integer (round-half-up).
#     Verified all three conditions (normal case, chromosome-start edge,
#     chromosome-end edge, overlap with another real cluster) on a test set
#     before shipping, including one case not originally designed for
#     (a cluster's downstream shift landing on a different real cluster,
#     not just itself) that the same exclusion rule correctly caught

# PARalus 1.6.21:
#   - stopped deleting "${name}".aligned_TtoC.bed in cleanup. This is the
#     per-READ TtoC bed (one row per read showing a real conversion),
#     before it gets deduped into aligned_TtoC_loc.bed (one row per unique
#     location, no depth information). Standalone downstream tools (e.g.
#     the metagene scripts) that want to weight sites by read depth/
#     confidence - rather than counting every site equally regardless of
#     how many reads support it - need this per-read file to compute that;
#     it was otherwise gone by the time any such tool could use it. No
#     other behavior changed - every other cleanup deletion is untouched

# PARalus 1.6.22:
#   - k-mer enrichment: bumped --top from 10 to 32 (both the labeled top/
#     bottom summary in _stats.txt and the labeled points on the scatter)
#   - added the 5'UTR/CDS/3'UTR and intron metagene tools directly into the
#     pipeline, run automatically at the end against three site sources:
#       - clusters (from clusters.csv) - unweighted, every site counts
#         equally, matching how clusters are already reported elsewhere
#       - groups (from groups.csv, the pre-TtoC-filter version of the same
#         data) - also unweighted. Using groups rather than only clusters
#         matters for plain CLIP data with no 4SU/6SG labeling at all:
#         clusters.csv requires TtoC support to promote a group into a
#         cluster, which can throw out nearly everything real when there's
#         no conversion signal to filter on in the first place
#       - TtoC locations - weighted by real per-location read depth
#         (recovered from the per-read aligned_TtoC.bed that 1.6.21 stopped
#         deleting specifically for this), so a location supported by 7
#         reads counts more than one supported by 1
#     Output: <name>_clusters_utr_cds.pdf, <name>_groups_utr_cds.pdf,
#     <name>_ttoc_utr_cds.pdf, <name>_clusters_intron.pdf,
#     <name>_groups_intron.pdf, <name>_ttoc_intron.pdf. A genuine
#     zero-qualifying-data case (e.g. too few clusters landing in intron
#     regions) is skipped with a clear note in _stats.txt rather than
#     producing an empty or misleading plot - verified this exact case on
#     synthetic data before shipping, not just the cases with real output
#   - verified the whole block against synthetic clusters.csv/groups.csv/
#     aligned_TtoC.bed matching the real 17-column format, then re-verified
#     by extracting and running the exact assembled code from this
#     deployed file (not just the standalone metagene scripts in
#     isolation), confirming the heredoc embedding didn't corrupt anything

# PARalus 1.6.23:
#   - dropped the clusters/groups metagene outputs, keeping only the
#     depth-weighted TtoC one: <name>_ttoc_utr_cds.pdf and
#     <name>_ttoc_intron.pdf. Cuts the added runtime from 6 metagene runs
#     down to 2, and TtoC was the one actually wanted going forward

# PARalus 1.6.24:
#   - 5-mer enrichment plot output changed from .png to .pdf
#     (<name>_5mer_enrichment.pdf), matching the metagene outputs. No
#     change to the plot itself or the underlying ranked tsv

# PARalus 1.6.25:
#   - cleanup now also removes the rRNA pre-mapping intermediates (only
#     present if -e was used): pre_rRNA_dep_*_trimmed.fasta,
#     rRNA_depleted_*_trimmed.fasta, and STAR's own small log files for
#     that run (SJ.out.tab, Log.progress.out, Log.out, Log.final.out).
#     Deliberately NOT deleting rRNA_aligned_*Aligned.sortedByCoord.out.bam
#     - that's the reads that actually mapped to rRNA, real data worth
#     keeping if you want to look at binding on rRNA itself, not just an
#     intermediate
#   - also removes "${name}".aligned_TtoC_loc.bed (the deduped, one-row-
#     per-location TtoC bed) - it's fully consumed by the coverage steps
#     earlier in the pipeline, and now reconstructable on demand from
#     aligned_TtoC.bed (kept since 1.6.21) via the same merge command
#     already used for the metagene weighting, so keeping both was
#     redundant. Verified with dummy files matching the exact real
#     filenames before shipping - confirmed the rRNA BAM and
#     aligned_TtoC.bed survive, everything else is removed

# --------------------------------------------------------------------------
# run_trierarch: Trierarch's annotation-building logic, as a function so it
# lives in the same process as Paralus - no separate file to locate, no
# subprocess fork. Every variable it touches is declared local here so it
# can never leak into or collide with Paralus's own variables of the same
# name (file_dir, name, etc. are used by both scripts for different things).
# All of Trierarch's own `exit` calls became `return` for the same reason -
# an exit here would kill the whole combined script, not just this function.
#
# Two ways to reach this function:
#   1) standalone: ./paralus.sh --trierarch -G genes.gtf -C -o output
#   2) automatic, from inside a normal Paralus run via -G (see below)
# --------------------------------------------------------------------------
run_trierarch() {
    local OPTIND=1
    local gtf tRNA_DB miRNA_DB outprefix
    local canonical_only=0 keep_unnamed=0
    local gtfname file_dir canonical_filter gff3_signal ncol OPT cls classfile remainder priority_order
    local prog scriptdir

# Trierarch (formerly "masterchef") - builds the annotation .bed file used by Paralus
# In Athens, the trierarch was the citizen responsible for equipping and
# outfitting a trireme before it sailed. This script equips the bed file
# Paralus runs on.
#
# v3.0 changes vs masterchef_v2:
#   - gene_type OR gene_biotype attribute now recognised (works with GENCODE and
#     Ensembl/other gtfs, not just gencode.vNN)
#   - added -M for an external miRNA database (miRBase gff3 or plain bed6), on top
#     of the existing -t for tRNA database (GtRNAdb bed6)
#   - overlap resolution is now a general priority cascade (see PRIORITY below)
#     instead of a single ad-hoc "subtract protein_coding from smallRNA" step.
#     This mirrors the exact annotation/biotype priority Paralus itself already
#     assumes downstream in its cluster-cleaning step, so the fixups Paralus used
#     to do with sed (3utr,tRNA -> exon,tRNA etc.) are no longer needed because
#     the bed file is already conflict-free coming out of this script.

unset ALL_PROXY
prog=$0
scriptdir=$(dirname "${prog}")

if [ $# -eq 0 ]; then
    echo 'Trierarch v3.0 (formerly masterchef)
Please give me something to do
No input files provided
Please type "--help" for instructions'
    return 0
fi

if [ "$1" = "--help" ] || [ "$1" = "-h" ]; then
    cat << 'EOF'
#####################################################
Usage:
  ./paralus1_6_25.sh --trierarch -G genes.gtf [-T tRNA_DB.bed] [-M miRNA_DB.gff3] [-O output_prefix]

  Required fields:

    -G   any gtf annotation (GENCODE, Ensembl, etc - gene_type or gene_biotype
         attribute is used, whichever is present)

  Optional fields:

    -T   tRNA bed6 (e.g. from GtRNAdb): chr,start,end,name,score,strand
    -M   miRNA database - either a miRBase gff3, or a plain bed6 the same shape
         as the tRNA file
    -O   output prefix (default: gtf filename minus .gtf)
    -C   canonical/MANE only - restrict to transcripts tagged
         "Ensembl_canonical" or "MANE_Select" in the gtf. Recommended for
         messy/well-populated annotations (e.g. GENCODE): cuts out noise from
         rarely-expressed alternative transcripts (extra/rare exons, stop
         codons far downstream of the real one, etc.) that otherwise distort
         the UTR/CDS calls for a gene. Off by default since it needs the gtf
         to actually carry those tags (most recent GENCODE/Ensembl do).
    -U   keep genes with no gene_name (falls back to gene_id, e.g. an
         ENSG00000... accession as the name). Off by default - these rows
         are dropped instead, since an accession-only name is usually a
         novel/unannotated locus you don't want cluttering the bed.

  Output:
    <prefix>_Paralus.bed  - the 7-column bed (chr,start,end,gene,annotation,
                             strand,biotype) ready for `Paralus -A`

  Downloading tRNA_DB / miRNA_DB from the terminal:

    tRNA_DB (GtRNAdb, e.g. human hg38):
      wget https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.bed
      wget https://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc39/mm39-tRNAs.bed     # mouse mm39
      wget https://gtrnadb.ucsc.edu/genomes/eukaryota/Schi_pomb_972h/schiPomb_972H-tRNAs.bed  # S. pombe
      wget https://gtrnadb.ucsc.edu/genomes/eukaryota/Dmela6/dm6-tRNAs.bed      # fruit fly dm6
      wget https://gtrnadb.ucsc.edu/genomes/bacteria/Esch_coli_K_12_MG1655/eschColi_K_12_MG1655-tRNAs.bed  # E. coli
      # bacteria live under genomes/bacteria/, not genomes/eukaryota/ - easy to miss
      # if a direct .bed 404s, fall back to the tarball (works for every genome):
      #   wget https://gtrnadb.ucsc.edu/genomes/eukaryota/<Folder>/<build>-tRNAs.tar.gz
      #   tar -xzf <build>-tRNAs.tar.gz
      # to find a folder/build name not listed here, browse
      # https://gtrnadb.ucsc.edu/ (Browse by Phylogeny) and read it off the URL

    miRNA_DB (miRBase, e.g. human):
      wget https://www.mirbase.org/download/hsa.gff3   # human
      wget https://www.mirbase.org/download/mmu.gff3   # mouse
      wget https://www.mirbase.org/download/dme.gff3   # fruit fly
      # -> <code>.gff3 is your -M file, used as-is
      # unlike tRNAs, miRNAs are eukaryote-only and not even universal there:
      #  - no entry for fungi like S. pombe/S. cerevisiae - most fungi lack
      #    the canonical miRNA pathway (no Drosha/DGCR8) entirely
      #  - no entry for bacteria (e.g. E. coli) - no miRNA pathway in
      #    prokaryotes at all, full stop
      # full list of the ~30 organisms miRBase has genome coordinates for:
      # https://www.mirbase.org/download/ (species code, e.g. mmu=mouse, rno=rat)
#####################################################
EOF
    return 0
fi

outprefix=""
canonical_only=0
keep_unnamed=0
while getopts hG:T:M:O:CUL: OPT; do
    case "$OPT" in
        G) gtf=$OPTARG ;;
        T) tRNA_DB=$OPTARG ;;
        M) miRNA_DB=$OPTARG ;;
        O) outprefix=$OPTARG ;;
        C) canonical_only=1 ;;
        U) keep_unnamed=1 ;;
        h) return 0 ;;
    esac
done

if [ ! -f "${gtf}" ]; then
    echo "error:
please provide a gtf file using the -G option"
    return 0
else
    echo "GTF file provided"
fi

# --------------------------------------------------------------------------
# GTF vs GFF3 guard.
# GFF3 and GTF look similar (both 9-column, tab-separated) but use a
# completely different attribute syntax: GTF uses `key "value";` pairs
# (quoted), GFF3 uses `key=value;` (unquoted, hierarchical ID=/Parent=
# references instead of repeating gene/transcript info per line). This
# script's attribute parsing assumes GTF's flat quoted-pair style. Fed a
# real GFF3 it won't error - it'll just silently extract nothing/garbage,
# which is worse than failing loudly. Catch it here instead.
# --------------------------------------------------------------------------
gff3_signal=$(grep -v "^#" "${gtf}" | awk -F "\t" 'NF>=9 && $0!=""{print $9; exit}')
if echo "${gtf}" | grep -qi '\.gff3?$' || ! echo "${gff3_signal}" | grep -q '"'; then
    echo "error:
this looks like GFF3, not GTF (no quoted attribute values found, e.g. gene_id \"X\"; -
GFF3 uses gene_id=X; instead). This script's parsing assumes GTF's flat quoted-pair
attribute style and will silently extract nothing useful from a real GFF3.

Convert it to GTF first, then re-run with the converted file:
  gffread input.gff3 -T -o output.gtf
(gffread: conda install -c bioconda gffread. For messier/non-standard gff3, AGAT's
agat_convert_sp_gff2gtf.pl is more robust: conda install -c bioconda agat)"
    return 1
fi

file_dir=$(dirname "${gtf}")

if [ -n "${outprefix}" ]; then
    case "${outprefix}" in
        */*) gtfname="${outprefix}" ;;                    # already a path - use as given
        *)   gtfname="${file_dir}/${outprefix}" ;;        # bare name - put it next to the gtf
    esac
else
    gtfname="${gtf%.gtf}"
fi

echo "extracting all features into a master bed, then building the Paralus annotation bed"

if [ "${canonical_only}" -eq 1 ]; then
    echo "restricting to Ensembl_canonical / MANE_Select transcripts only (-C) - trims rare/lowly-expressed isoform noise"
fi

# --------------------------------------------------------------------------
# 1) pull exon/start_codon/stop_codon/CDS lines, normalise the attribute
#    column, and grab gene_name + whichever of gene_type/gene_biotype exists
# --------------------------------------------------------------------------
canonical_filter='cat'
if [ "${canonical_only}" -eq 1 ]; then
    canonical_filter='grep -i -w -E "tag \"Ensembl_canonical\"|tag \"MANE_Select\""'
fi

cat "$gtf" | awk -F "\t" 'tolower($3) == "exon" || tolower($3) == "start_codon" || tolower($3) == "stop_codon" || tolower($3) == "cds"' \
     | awk -v canonical="${canonical_only}" '
        canonical == 0 { print; next }
        /gene_type "protein_coding"/ || /gene_biotype "protein_coding"/ {
            if (/tag "Ensembl_canonical"/ || /tag "MANE_Select"/) print
            next
        }
        { print }
    ' \
    | grep -v "##" | sed 's/"//g' | sed 's/;//g' | sed 's/[[:space:]]/\t/g' | sed 's/\t\t/\t/g' | sed 's/"="/\t/g' \
    | awk -v keep_unnamed="${keep_unnamed}" '{
        OFS="\t";
        n=""; t="";
        for (I=1; I<NF; I++) {
            if ($I == "gene_name") n=$(I+1);
            if ($I == "gene_type" || $I == "gene_biotype") t=$(I+1);
        }
        if (n == "") {
            if (keep_unnamed == 1) { for (I=1; I<NF; I++) if ($I=="gene_id") n=$(I+1) }
            else next
        }
        if (t == "") t="unknown";
        print $1,$4-1,$5,n,tolower($3),$7,t
    }' \
    | sed 's/,/-/g' | bedtools sort | uniq \
    | awk -F "\t" '$5 !~ "UTR"' | awk -F "\t" 'toupper($4) !~ "Y_RNA"' | awk -F "\t" 'toupper($4) !~ "PRNA"' \
    | grep -v -i transcript | grep -v -i "antisense\|Metazoa_SRP\|misc_RNA\|sense_intronic\|-AS1\|mascRNA-menRNA\|pseudogene" \
    > "${gtfname}"_gtf.bed

echo "extracting gene names..."
awk -F "\t" '{print $4}' "${gtfname}"_gtf.bed | sort | uniq > "${gtfname}"_genes

echo "grouping annotation to create full gene body bed"
bedtools groupby -i <(sort -k4,4 -k2,2n "${gtfname}"_gtf.bed) -g 4 -c 1,2,3,4,5,6,7 -o distinct,min,max,distinct,distinct,distinct,distinct \
    | awk -F "\t" '{OFS="\t"; print $2,$3,$4,$5,$6,$7,$8}' > "${gtfname}"_genes_1.bed
uniq "${gtfname}"_genes_1.bed | awk -F "\t" '{OFS="\t"; print $1,$2,$3,$4,"gene",$6,$7}' | bedtools sort | uniq \
    > "${gtfname}"_genes.bed

echo "merging start/stop codons"
grep -i "start_codon\|stop_codon" "${gtfname}"_gtf.bed \
    | bedtools merge -s -c 4,5,6,7 -o distinct,distinct,distinct,distinct > "${gtfname}"_start_stop.bed

echo "removing small RNA gene bodies before building introns"
grep -v -w -i "miRNA\|snoRNA\|snRNA\|sRNA" "${gtfname}"_gtf.bed > "${gtfname}"_gtf_nomirs.bed
grep -v -w -i "miRNA\|snoRNA\|snRNA\|sRNA" "${gtfname}"_genes.bed > "${gtfname}"_genes_nomirs.bed
bedtools subtract -s -a "${gtfname}"_genes_nomirs.bed -b "${gtfname}"_gtf_nomirs.bed \
    | awk -F "\t" '{OFS="\t"; print $1,$2,$3,$4,"intron",$6,$7}' | bedtools sort > "${gtfname}"_introns.bed

echo "extracting CDS and splice sites"
awk -F "\t" '$5 == "cds"' "${gtfname}"_gtf.bed | bedtools sort \
    | bedtools merge -s -c 4,5,6,7 -o distinct,distinct,distinct,distinct \
    | awk -F "\t" '{OFS="\t"; print $1,$2,$3,$4,"CDS",$6,$7}' > "${gtfname}"_cds.bed

# --------------------------------------------------------------------------
# creating 3utr / 5utr
#
# The old approach took the widest possible span across ALL transcript
# isoforms (earliest start_codon to latest stop_codon) and called that the
# UTR. On real multi-isoform genes that span massively overshoots the true
# UTR and mostly overlaps CDS territory from other isoforms - which is why
# forcing "CDS beats UTR" on it wiped out UTRs almost everywhere. Building
# it that way was the actual bug, not the overlap-resolution step.
#
# Correct approach: take the union of all real exonic bases (across every
# isoform), subtract the union of all CDS bases (across every isoform).
# What's left is genuine non-coding exonic sequence, and it can never
# overlap CDS because it's built by removing CDS directly. Each remaining
# piece is then called 5utr/3utr/exon depending on where it sits relative
# to the gene's overall CDS span (its "envelope") and strand. A piece that
# falls inside the CDS envelope but isn't CDS itself (e.g. a cassette exon
# that's coding in one isoform, skipped in another) is ambiguous - it gets
# labelled "exon" rather than guessed at, since Paralus already treats
# "exon" as its own top-priority annotation.
# --------------------------------------------------------------------------
echo "creating 3utr / 5utr"
awk -F "\t" '$5 == "exon"' "${gtfname}"_gtf.bed | bedtools sort \
    | bedtools merge -s -c 4,5,6,7 -o distinct,distinct,distinct,distinct > "${gtfname}"_exon_merged.bed

bedtools subtract -s -a "${gtfname}"_exon_merged.bed -b "${gtfname}"_cds.bed | bedtools sort \
    > "${gtfname}"_noncoding_exonic.bed

sort -k4,4 -k2,2n "${gtfname}"_cds.bed \
    | bedtools groupby -g 4 -c 1,2,3,4,6,7 -o distinct,min,max,distinct,distinct,distinct \
    | awk -F "\t" '{OFS="\t"; print $2,$3,$4,$1,$5,$6}' > "${gtfname}"_cds_envelope.bed

awk -F "\t" '
    NR==FNR { env_start[$4]=$2; env_end[$4]=$3; next }
    {
        OFS="\t";
        gene=$4; strand=$6;
        if (!(gene in env_start)) { print $1,$2,$3,$4,"exon",$6,$7; next }
        es=env_start[gene]; ee=env_end[gene];
        if (strand == "+") {
            if ($3 <= es) label="5utr";
            else if ($2 >= ee) label="3utr";
            else label="exon";
        } else if (strand == "-") {
            if ($2 >= ee) label="5utr";
            else if ($3 <= es) label="3utr";
            else label="exon";
        } else label="exon";
        print $1,$2,$3,$4,label,$6,$7;
    }' "${gtfname}"_cds_envelope.bed "${gtfname}"_noncoding_exonic.bed \
    | awk -F "\t" '$3-$2 > 1' | bedtools sort > "${gtfname}"_3utr_5utr.bed

rm -f "${gtfname}"_exon_merged.bed "${gtfname}"_noncoding_exonic.bed "${gtfname}"_cds_envelope.bed

awk -F "\t" '{ if ($6=="+"){OFS="\t"; print $1,$2-1,$2+1,$4,"5ss",$6,$7"\n"$1,$3-1,$3+1,$4,"3ss",$6,$7}
               else if ($6=="-"){OFS="\t"; print $1,$2-1,$2+1,$4,"3ss",$6,$7"\n"$1,$3-1,$3+1,$4,"5ss",$6,$7} }' \
    "${gtfname}"_introns.bed > "${gtfname}"_5ss_3ss.bed

merge_class () {
    local infile outfile
    infile="$1"; outfile="$2"
    if [ -s "${infile}" ]; then
        bedtools sort -i "${infile}" \
            | bedtools merge -s -c 4,5,6,7 -o distinct,distinct,distinct,distinct > "${outfile}"
    else
        : > "${outfile}"
    fi
}

echo "extracting small RNA gene bodies from the gtf itself"
grep -w -i "miRNA\|snoRNA\|snRNA\|sRNA" "${gtfname}"_gtf.bed > "${gtfname}"_mirs_raw.bed
merge_class "${gtfname}"_mirs_raw.bed "${gtfname}"_gtf_mirs.bed
rm -f "${gtfname}"_mirs_raw.bed

# --------------------------------------------------------------------------
# 2) external small-RNA databases: tRNA_DB (existing) + miRNA_DB (new)
#    both accepted as bed6 (chr,start,end,name,score,strand); a miRBase gff3
#    is auto-detected and converted.
# --------------------------------------------------------------------------
if [ -n "${tRNA_DB}" ] && [ -f "${tRNA_DB}" ]; then
    echo "preparing tRNA_DB"
    awk -F "\t" '{OFS="\t"; print $1,$2,$3,$4,"exon",$6,"tRNA"}' "${tRNA_DB}" > "${gtfname}"_tRNA_ext.bed
else
    : > "${gtfname}"_tRNA_ext.bed
fi

if [ -n "${miRNA_DB}" ] && [ -f "${miRNA_DB}" ]; then
    echo "preparing miRNA_DB"
    ncol=$(grep -v "^#" "${miRNA_DB}" | awk -F "\t" 'NF>0{print NF; exit}')
    if [ "${ncol}" -ge 8 ]; then
        # looks like a gff3 (miRBase): chr,source,type,start,end,score,strand,frame,attrs
        # Use only the primary_transcript (hairpin/precursor) row, not the mature
        # -3p/-5p "miRNA" sub-features - those sit inside the same interval and,
        # if kept, get merged into one locus with 3 names stitched together
        # ("hairpin|hsa-miR-x-3p|hsa-miR-x-5p"). One row per locus is what we want.
        awk -F "\t" '$3 == "miRNA_primary_transcript"' "${miRNA_DB}" \
            | sed 's/;/\t/g' \
            | awk -F "\t" '{
                OFS="\t"; name="";
                for (I=9; I<=NF; I++) { if ($I ~ /^Name=/) { name=$I; sub("Name=","",name) } }
                if (name=="") name="miRNA_"NR;
                print $1,$4-1,$5,name,"exon",$7,"miRNA"
              }' > "${gtfname}"_miRNA_ext.bed
    else
        # plain bed6, same convention as tRNA_DB
        awk -F "\t" '{OFS="\t"; print $1,$2,$3,$4,"exon",$6,"miRNA"}' "${miRNA_DB}" > "${gtfname}"_miRNA_ext.bed
    fi
else
    : > "${gtfname}"_miRNA_ext.bed
fi

# --------------------------------------------------------------------------
# 3) build one bed per biotype class, ranked by priority (highest wins the
#    overlapping bases). This is the same order Paralus already expects in
#    its own cluster-cleaning step, just enforced here instead:
#      snoRNA > snRNA > rRNA > miRNA > tRNA > Mt_tRNA > scaRNA > lncRNA >
#      protein_coding > everything else
# --------------------------------------------------------------------------
echo "building per-biotype class beds"

grep -w -i snoRNA "${gtfname}"_gtf_mirs.bed > "${gtfname}"_class_snoRNA.bed
grep -w -i snRNA  "${gtfname}"_gtf_mirs.bed > "${gtfname}"_class_snRNA.bed
grep -w -i -E "rRNA|Mt_rRNA" "${gtfname}"_gtf.bed > "${gtfname}"_rRNA_raw.bed
merge_class "${gtfname}"_rRNA_raw.bed "${gtfname}"_class_rRNA.bed
if [ -n "${miRNA_DB}" ] && [ -f "${miRNA_DB}" ]; then
    # an external miRNA_DB was given - use it exclusively for the miRNA class.
    # Mixing it with the gtf's own miRNA gene rows just produces one merged
    # locus carrying both names stitched together, which is worse than either
    # source alone.
    grep -w -i miRNA "${gtfname}"_miRNA_ext.bed > "${gtfname}"_miRNA_raw.bed
else
    grep -w -i miRNA "${gtfname}"_gtf_mirs.bed > "${gtfname}"_miRNA_raw.bed
fi
merge_class "${gtfname}"_miRNA_raw.bed "${gtfname}"_class_miRNA.bed
cat "${gtfname}"_tRNA_ext.bed | grep -w -i -v Mt_tRNA | bedtools sort > "${gtfname}"_class_tRNA.bed
grep -w -i Mt_tRNA "${gtfname}"_gtf_mirs.bed > "${gtfname}"_class_Mt_tRNA.bed
grep -w -i scaRNA "${gtfname}"_gtf_mirs.bed > "${gtfname}"_class_scaRNA.bed
grep -w -i -E "lncRNA|lincRNA" "${gtfname}"_gtf.bed > "${gtfname}"_lncRNA_raw.bed
merge_class "${gtfname}"_lncRNA_raw.bed "${gtfname}"_class_lncRNA.bed
rm -f "${gtfname}"_rRNA_raw.bed "${gtfname}"_miRNA_raw.bed "${gtfname}"_lncRNA_raw.bed

echo "building the protein_coding class (utrs, cds, introns, splice sites, start/stop)"
cat "${gtfname}"_3utr_5utr.bed "${gtfname}"_5ss_3ss.bed "${gtfname}"_cds.bed "${gtfname}"_start_stop.bed \
    "${gtfname}"_introns.bed | grep -w -i protein_coding | bedtools sort | uniq > "${gtfname}"_class_protein_coding.bed

echo "building the catch-all class for anything not covered above"
cat "${gtfname}"_gtf.bed \
    | grep -v -w -i "snoRNA\|snRNA\|rRNA\|Mt_rRNA\|miRNA\|tRNA\|Mt_tRNA\|scaRNA\|lncRNA\|lincRNA\|protein_coding\|CDS\|start_codon\|stop_codon" \
    | bedtools sort | uniq > "${gtfname}"_class_other.bed

# --------------------------------------------------------------------------
# 4) priority cascade: subtract everything already claimed by a higher
#    priority class from every lower priority class, strand-aware, then
#    accumulate. Order = priority list above.
# --------------------------------------------------------------------------
echo "resolving overlaps by priority (non-coding wins)"
priority_order="snoRNA snRNA rRNA miRNA tRNA Mt_tRNA scaRNA lncRNA protein_coding other"

: > "${gtfname}"_Paralus.bed
: > "${gtfname}"_claimed.bed

for cls in ${priority_order}; do
    classfile="${gtfname}_class_${cls}.bed"
    if [ ! -s "${classfile}" ]; then
        continue
    fi
    if [ -s "${gtfname}"_claimed.bed ]; then
        remainder=$(mktemp)
        bedtools subtract -s -a "${classfile}" -b "${gtfname}"_claimed.bed > "${remainder}"
    else
        remainder=$(mktemp)
        cp "${classfile}" "${remainder}"
    fi
    awk -F "\t" '$3-$2 > 0' "${remainder}" >> "${gtfname}"_Paralus.bed
    cat "${remainder}" >> "${gtfname}"_claimed.bed
    bedtools sort -i "${gtfname}"_claimed.bed | bedtools merge -s -c 4,5,6,7 -o distinct,distinct,distinct,distinct \
        > "${gtfname}"_claimed_merged.bed
    mv "${gtfname}"_claimed_merged.bed "${gtfname}"_claimed.bed
    rm -f "${remainder}"
done

bedtools sort -i "${gtfname}"_Paralus.bed | uniq | sed 's/,/|/g' > "${gtfname}"_Paralus_sorted.bed
mv "${gtfname}"_Paralus_sorted.bed "${gtfname}"_Paralus.bed

echo "final Paralus annotation bed written to: ${gtfname}_Paralus.bed"

# --------------------------------------------------------------------------
# cleanup intermediates
# --------------------------------------------------------------------------
rm -f "${gtfname}"_gtf.bed "${gtfname}"_genes.bed "${gtfname}"_genes_1.bed "${gtfname}"_3utr_5utr.bed \
      "${gtfname}"_start_stop.bed "${gtfname}"_gtf_nomirs.bed "${gtfname}"_introns.bed \
      "${gtfname}"_gtf_mirs.bed "${gtfname}"_cds.bed "${gtfname}"_5ss_3ss.bed "${gtfname}"_genes \
      "${gtfname}"_genes_nomirs.bed "${gtfname}"_tRNA_ext.bed \
      "${gtfname}"_miRNA_ext.bed "${gtfname}"_claimed.bed \
      "${gtfname}"_class_snoRNA.bed "${gtfname}"_class_snRNA.bed "${gtfname}"_class_rRNA.bed \
      "${gtfname}"_class_miRNA.bed "${gtfname}"_class_tRNA.bed "${gtfname}"_class_Mt_tRNA.bed \
      "${gtfname}"_class_scaRNA.bed "${gtfname}"_class_lncRNA.bed "${gtfname}"_class_protein_coding.bed \
      "${gtfname}"_class_other.bed

echo
#################################
#####         DONE        #######
#################################
}

unset ALL_PROXY

prog=$0
scriptdir=$(dirname "${prog}")

if [ "$1" = "--trierarch" ]; then
    shift
    run_trierarch "$@"
    exit $?
fi

if [ $# -eq 0 ]; then
  echo 'PARalus v1.6.25
Please give me something to do 
No input files provided
Please type "--help" for instructions'
  exit 0
fi

if [ "$1" = "--help" ] || [ "$1" = "-h" ]; then
  cat << 'EOF'
  
          ||
       ##################
       ##################   
       ###Paralus1.6.25###
       ##################
       ##################
      ##################
           ||         
   &  __@_/  ||       &
   \ \  ||    ||       //
    \__\_/\___@___@___@___@____//   
      \   /  /  /  /  /
  ~~~~~~~~\____/___/___/___/____))))~~~~~~~~~~
   ~ ~~   /  /  /  /   ~~   ~
   ~~    /  /  /  /  ~    ~~~
 
 Usage: ./Paralus -x genome_index_base -i input_file -A annotation.bed
        ./Paralus -x genome_index_base -i input_file -G genes.gtf
        ./Paralus --trierarch -G genes.gtf -C -O output_prefix    (Trierarch standalone, no Paralus run)
 
 Required fields:
  
  -i  Your input file (fastqz/fastq/fasta/bam/SRA) 
  -x  Genome index folder
    		May build one like this:
  			make dir /path/Genome_index
  			STAR --runMode genomeGenerate --genomeDir /path/Genome_index --genomeFastaFiles rRNA_genome.fasta

  Annotation - provide EITHER -A (a ready-made bed) OR -G (a raw gtf, Trierarch
  runs automatically - takes about a minute, cached after the first run). Or
  run Trierarch on its own with --trierarch, e.g. to pre-build a bed once and
  reuse it across many samples with -A - see "./Paralus --trierarch --help"
  for Trierarch's own flags:

  -A  A 7-column bed file already produced by Trierarch. Use this if you've
      already run Trierarch yourself, or want to reuse the same annotation
      across many samples without re-checking it every time.
  -G  A raw gtf (GENCODE or Ensembl format, gene_type or gene_biotype either
      recognised). Trierarch is run automatically to build the annotation bed
      the first time; on later runs with the same gtf and flags below, the
      cached *_Paralus.bed next to the gtf is reused instead of rebuilding it.
      Flags below are the same letters Trierarch uses natively (see
      "./Paralus --trierarch --help") and are passed straight through when
      -G triggers it - no translation between standalone and -G usage:
        -T  tRNA database bed6 (e.g. from GtRNAdb)
        -M  miRNA database, miRBase gff3 or bed6
        -C  restrict to Ensembl_canonical/MANE_Select transcripts
        -U  keep genes with no gene_name (fallback to gene_id)
      (uppercase M/C/U since lowercase m/c/u are already Paralus's own
      options below - Trierarch itself uses M/C/U natively for exactly
      this reason, whether run standalone or auto-triggered via -g)

  *Make sure you have loaded samtools bedtools and cutadapt acording to your enviroment
  
  
  
  
 Optional parameters:
  
  Preprocess:
  
  -n  Minimum length allowed after trimming. Default: 19
  -m  Minimum number of aligned reads to consider. Default: 3
  -a  Your 3' adapter. Default: No adapter
  -d  Your 5' adapter. Default: GTTCAGAGTTCTACAGTCCGACGATC
  -u  Number of wildcards anchored to the 5' end of the sequence. Default: 0
     (do not add this if 5' adapter is barcoded)
  -p  Threads. Default: 8 dont use more than 12!! else star wont run
  -e  rRNA genome to pre map sequnces and provide % of rRNA reads 
  		May build one like this:
  		make dir /path/rRNA_genome
  		STAR --runMode genomeGenerate --genomeDir /path/rRNA_genome --genomeFastaFiles rRNA_genome.fasta
  
  PAR-CLIP binding site identification - parameters:
  
  -c  Conversion location count. How many TtoC location to be considered cluster.(1-2) Default: 1
  -r  A group will be considered a cluster if TtoC containing reads are equal or more than this fraction: 0.3 (30%) 
  -f  If genome.fa (<Genome index basename>.fa) does not exist in the subfolder of genome index,
     please provide a the genome.fa using the -f option or copy the genome.fa into the genome index folder
  
#######################
EOF
  exit 0
fi

minlen=19
minRead=3
five_adapt=GTTCAGAGTTCTACAGTCCGACGATC
umis=0
f=fq
threads=8
convloc=1
fraction=0.3

while getopts hi:x:A:p:r:v:n:m:a:d:u:s:e:c:f:L:G:T:M:CU OPT; do
  case "$OPT" in
    i)
      input=$OPTARG
      ;;
    x)
		index=$OPTARG
      ;;
    A)
      annotation=$OPTARG
      ;;
    G)
      raw_gtf=$OPTARG
      ;;
    T)
      trierarch_tRNA_DB=$OPTARG
      ;;
    M)
      trierarch_miRNA_DB=$OPTARG
      ;;
    C)
      trierarch_canonical=1
      ;;
    U)
      trierarch_keep_unnamed=1
      ;;
    p)
      threads=$OPTARG
      ;;
    r)
      fraction=$OPTARG
      ;;
    n)
      minlen=$OPTARG
      ;;
    m)
      minRead=$OPTARG
      ;;
    a)
      three_adapt=$OPTARG
      ;;
    d)
      five_adapt=$OPTARG
      ;;
    u)
      umis=$OPTARG
      ;;

    c)
    	convloc=$OPTARG
      ;;
    f)
    	genomefasta=$OPTARG
      ;;
     e)
    	rRNA_genome=$OPTARG
      ;;

    h) exit 0;;
  esac
  
done
input=$(basename "${input}")
name=${input%.fast*}

file_dir=$(dirname "${input}")


ext="${input##*.}"


genomedir=$(dirname "${index}")
rRNA_genomedir=$(dirname "${rRNA_genome}")

# --------------------------------------------------------------------------
# -G integration: if a raw gtf was given instead of a ready-made -A bed,
# run Trierarch automatically to build one. Resolve paths to absolute here,
# before any `cd`, so a relative -g/-t/-M keeps working regardless of which
# directory Paralus itself changes into below.
# --------------------------------------------------------------------------
if [ -n "${raw_gtf}" ]; then
    raw_gtf=$(cd "$(dirname "${raw_gtf}")" && pwd)/$(basename "${raw_gtf}")
fi
if [ -n "${trierarch_tRNA_DB}" ]; then
    trierarch_tRNA_DB=$(cd "$(dirname "${trierarch_tRNA_DB}")" && pwd)/$(basename "${trierarch_tRNA_DB}")
fi
if [ -n "${trierarch_miRNA_DB}" ]; then
    trierarch_miRNA_DB=$(cd "$(dirname "${trierarch_miRNA_DB}")" && pwd)/$(basename "${trierarch_miRNA_DB}")
fi

if [ -z "${annotation}" ] && [ -n "${raw_gtf}" ]; then
    expected_annotation="${raw_gtf%.gtf}_Paralus.bed"

    need_rebuild=1
    if [ -f "${expected_annotation}" ]; then
        need_rebuild=0
        [ "${raw_gtf}" -nt "${expected_annotation}" ] && need_rebuild=1
        [ -n "${trierarch_tRNA_DB}" ] && [ "${trierarch_tRNA_DB}" -nt "${expected_annotation}" ] && need_rebuild=1
        [ -n "${trierarch_miRNA_DB}" ] && [ "${trierarch_miRNA_DB}" -nt "${expected_annotation}" ] && need_rebuild=1
    fi

    if [ "${need_rebuild}" -eq 1 ]; then
        echo "-G given: building annotation from ${raw_gtf} with Trierarch (about a minute)..."
        trierarch_call=(-G "${raw_gtf}")
        [ -n "${trierarch_tRNA_DB}" ] && trierarch_call+=(-T "${trierarch_tRNA_DB}")
        [ -n "${trierarch_miRNA_DB}" ] && trierarch_call+=(-M "${trierarch_miRNA_DB}")
        [ "${trierarch_canonical}" = "1" ] && trierarch_call+=(-C)
        [ "${trierarch_keep_unnamed}" = "1" ] && trierarch_call+=(-U)
        run_trierarch "${trierarch_call[@]}"
        if [ ! -f "${expected_annotation}" ]; then
            echo "error:
Trierarch ran but ${expected_annotation} was not produced - check the Trierarch output
above for errors."
            exit 0
        fi
    else
        echo "-G given: reusing cached annotation ${expected_annotation} (gtf/databases unchanged since it was last built)"
    fi

    annotation="${expected_annotation}"
fi

cd $file_dir

# --------------------------------------------------------------------------
# sanity-check the annotation file BEFORE running anything: it needs to be
# Trierarch's bed output, not the raw GTF. A GTF has 9 columns with a
# free-text attribute field and looks superficially bed-like, so this can
# otherwise run for 10+ minutes of STAR/samtools before crashing deep in a
# bedtools step with a cryptic "differing number of BED fields" error.
# --------------------------------------------------------------------------
annot_bad=0
if [ -f "${annotation}" ]; then
    annot_ncol=$(awk -F "\t" 'NF>0 && $0 !~ /^#/{print NF; exit}' "${annotation}")
    if [ "${annot_ncol}" = "9" ]; then annot_bad=1; fi
    if ! awk -F "\t" 'NF>0 && $0 !~ /^#/{ if ($2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/) exit 0; else exit 1 }' "${annotation}"; then
        annot_bad=1
    fi
fi

if [ ! -f "${index}/SAindex" ]; then
echo "error:
STAR genome index not found at ${index}/SAindex - build one first (see --help) or
check the -x path"
 exit 0

elif [ -z "${genomefasta}" ] || [ ! -f "${genomefasta}" ]; then
echo "error:
please provide a genome.fa using the -f option (needed for samtools calmd and
bedtools getfasta later in the pipeline, even if a STAR index already exists)"
 exit 0

elif [ ! -f "${annotation}" ]; then
echo "error:
please provide an annotation bed with -A, or a raw gtf with -G to have Trierarch
build one automatically"
 exit 0

elif [ "${annot_bad}" = "1" ]; then
echo "error:
-A does not look like a proper annotation bed file (9 tab-separated columns,
or non-numeric start/end) - this usually means a raw .gtf was passed to -A.
Either pass the gtf with -G instead (Trierarch then runs automatically), or
run Trierarch yourself first:
  ./paralus1_6_25.sh --trierarch -G your.gtf -C -O output_prefix
and pass the resulting *_Paralus.bed here with -A."
 exit 0

elif [ ! -f "${input}" ]; then
echo "error:
please provide an input file (fastqz/fastq/fasta/"
 exit 0
 elif [[ ! -f "${rRNA_genome}/SAindex" ]]; then
echo "Warning:
rRNA_genome was not provided"
else
echo "Input files provided"
fi

# --------------------------------------------------------------------------
# chromosome whitelist: alt/unplaced scaffolds (e.g. GL000008.2-style
# GenBank accessions in a GENCODE primary_assembly genome) essentially never
# carry real gene annotations, so they show up in the annotation bed almost
# never - meaning they're also not analytically useful here regardless.
# Rather than hardcode a naming pattern like "^chr" (which would silently
# break on any genome that doesn't use that prefix, e.g. plain Ensembl-style
# "1","2","X"), build the whitelist from whatever contigs actually appear in
# the annotation file itself - organism/naming-convention agnostic, and it's
# also what eliminates the "inconsistent naming convention" bedtools
# warnings, since those alt scaffolds never enter any bed file from here on.
# --------------------------------------------------------------------------
awk -F "\t" '{print $1}' "${annotation}" | sort -u > "${name}"_valid_chroms.txt

if [ -f "${name}"_stats.txt ] ; then rm -f "${name}"_stats.txt ; fi 
if [ -f "${name}".log ] ; then rm -f "${name}".log ; fi

# --------------------------------------------------------------------------
# full run log: capture everything (stdout+stderr) from here on, plus the
# exact command that was invoked - so a run can always be traced back to
# precisely which flags/inputs produced it.
# --------------------------------------------------------------------------
runlog="${name}"_run_$(date +%Y%m%d_%H%M%S).log
exec > >(tee -a "${runlog}") 2>&1
echo "COMMAND: $0 $*"
echo "RUN DATE: $(date)"
pipeline_start_ts=$(date +%s)

echo "Your parameters:" >> "${name}"_stats.txt
echo "COMMAND : $0 $*" >> "${name}"_stats.txt
echo "Input file : "${input}"" >> "${name}"_stats.txt
echo "genome file : "${index}"" >> "${name}"_stats.txt
echo "annotation file : "${annotation}"" >> "${name}"_stats.txt
echo "threads : "${threads}"" >> "${name}"_stats.txt
echo "Minimun TtoC fraction : "${fraction}"" >> "${name}"_stats.txt
echo "Minimun conversion locations : "${convloc}"" >> "${name}"_stats.txt
echo "Minimun reads for group definition : "${minRead}"" >> "${name}"_stats.txt
echo "Minumun length allowed by cutadapt : "${minlen}"" >> "${name}"_stats.txt
echo "Your 3' adapter : "${three_adapt}"" >> "${name}"_stats.txt
echo "Your 5' adapter : "${five_adapt}"" >> "${name}"_stats.txt
echo "bases trimed from the 5' end of the read ( wildcards after 5' adapter) : "${umis}"" >> "${name}"_stats.txt
num_of_lines=$(wc -l < "${input}")
echo "number of lines in fastq: ${num_of_lines}" >> "${name}"_stats.txt
num_of_seq=$(echo $(( num_of_lines / 4 )))
echo "number of sequences : ${num_of_seq}" >> "${name}"_stats.txt

##########################################################################################################################
## preprocessing #in theory adapters should be removed prior 
## collapsing but error rates of illumina are < 0.1% making it insignificant #
##########################################################################################################################


if [[ ! -f ${rRNA_genome}/SAindex && $ext == "gz" && ${#three_adapt} -gt 0 && $umis -gt 0 ]]
	then
		echo "You have a fastq.gz file and have UMIs and adapters and did not provide a rRNA genome"
		gzip -dc "${name}".fastq.gz | awk 'NR%4==2' | sort \
		| uniq -c | sed 's/[[:space:]][[:space:]]//g' | awk '{OFS= "\n"; print ">"NR"-"$1,$2}' > ${name}_collapsed.fasta
		cutadapt -u "$umis" -m "$minlen" -a "$three_adapt" -g "$five_adapt;min_overlap=10" \
		 -o "${name}"_trimmed.fasta ${name}_collapsed.fasta >> "${name}"_stats.txt	

  	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --runThreadN "$threads" --genomeDir "${index}" \
	--readFilesType Fastx --readFilesIn "${name}"_trimmed.fasta --outFileNamePrefix "${name}." --outSAMtype BAM SortedByCoordinate \
	--outSAMmultNmax 1 --outMultimapperOrder Random --outSAMattributes NH MD HI AS 

elif [[ ! -f ${rRNA_genome}/SAindex && $ext == "gz" && ${#three_adapt} -gt 0 && $umis -eq 0 ]]
	then
		echo "You have a fastq.gz file, provided adapters and did not provide UMIs and rRNA genome"
		gzip -dc "${name}".fastq.gz | awk 'NR%4==2' | sort \
		| uniq -c | sed 's/[[:space:]][[:space:]]//g' | awk '{OFS= "\n"; print ">"NR"-"$1,$2}' > ${name}_collapsed.fasta
		cutadapt -m "$minlen" -a "$three_adapt" -g "$five_adapt;min_overlap=10" \
		 -o "${name}"_trimmed.fasta ${name}_collapsed.fasta >> "${name}"_stats.txt	

  	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --runThreadN "$threads" --genomeDir "${index}" \
	--readFilesType Fastx --readFilesIn "${name}"_trimmed.fasta --outFileNamePrefix "${name}." --outSAMtype BAM SortedByCoordinate \
	--outSAMmultNmax 1 --outMultimapperOrder Random --outSAMattributes NH MD HI AS 


elif [[ ! -f ${rRNA_genome}/SAindex && $ext == "fastq" && ${#three_adapt} -gt 0 && $umis -gt 0 ]]
	then
		echo "You have a fastq file and have UMIs and adapters and did not provide a rRNA genome"
		awk 'NR%4==2' "${name}".fastq | sort | uniq -c | sed 's/[[:space:]][[:space:]]//g' \
		| awk '{OFS= "\n"; print ">"NR"-"$1,$2}' > ${name}_collapsed.fasta
		cutadapt -u "$umis" -m "$minlen" -a "$three_adapt" -g "$five_adapt;min_overlap=10" \
		-o "${name}"_trimmed.fasta ${name}_collapsed.fasta >> "${name}"_stats.txt	

  	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --runThreadN "$threads" --genomeDir "${index}" \
	--readFilesType Fastx --readFilesIn "${name}"_trimmed.fasta --outFileNamePrefix "${name}." --outSAMtype BAM SortedByCoordinate \
	--outSAMmultNmax 1 --outMultimapperOrder Random --outSAMattributes NH MD HI AS 

elif [[ ! -f ${rRNA_genome}/SAindex && $ext == "fasta" && ${#three_adapt} -eq 0 && $umis -eq 0 ]]
	then
		echo "You have a fasta file and did not provide UMIs and adapters and a rRNA genome"
  	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --runThreadN "$threads" --genomeDir "${index}" \
	--readFilesType Fastx --readFilesIn "${name}"_trimmed.fasta --outFileNamePrefix "${name}." --outSAMtype BAM SortedByCoordinate \
	--outSAMmultNmax 1 --outMultimapperOrder Random --outSAMattributes NH MD HI AS 

elif [[ -f ${rRNA_genome}/SAindex && $ext == "gz" && ${#three_adapt} -gt 0 && $umis -gt 0 ]]
	then
				echo "You have a fastq.gz file and have UMIs and adapters and provided a rRNA genome"
		gzip -dc "${name}".fastq.gz | awk 'NR%4==2' | sort | uniq -c \
		| sed 's/[[:space:]][[:space:]]//g' | awk '{OFS= "\n"; print ">"NR"-"$1,$2}' > ${name}_collapsed.fasta
		cutadapt -u "$umis" -m "$minlen" -a "$three_adapt" -g "$five_adapt;min_overlap=10" \
		-o pre_rRNA_dep_"${name}"_trimmed.fasta ${name}_collapsed.fasta >> "${name}"_stats.txt

	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --genomeDir "${rRNA_genome}" --readFilesIn pre_rRNA_dep_"${name}"_trimmed.fasta \
	--readFilesType Fastx --outFileNamePrefix rRNA_aligned_"${name}" --outSAMtype BAM SortedByCoordinate --outSAMattributes NH MD HI AS \
	--runThreadN "$threads" --outReadsUnmapped Fastx --limitBAMsortRAM 2000000000
mv rRNA_aligned_"${name}"Unmapped.out.mate1 rRNA_depleted_"${name}"_trimmed.fasta 2>/dev/null || mv rRNA_aligned_"${name}"Unmapped.out.mate0 rRNA_depleted_"${name}"_trimmed.fasta 

  	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --runThreadN "$threads" --genomeDir "${index}" \
	--readFilesType Fastx --readFilesIn rRNA_depleted_"${name}"_trimmed.fasta --outFileNamePrefix "${name}." --outSAMtype BAM SortedByCoordinate \
	--outSAMmultNmax 1 --outMultimapperOrder Random --outSAMattributes NH MD HI AS 
		echo "done!" >> "${name}".log
elif [[ -f ${rRNA_genome}/SAindex && $ext == "fastq" && ${#three_adapt} -gt 0 && $umis -gt 0 ]]
	then
			echo "You have a fastq file and have UMIs and adapters and provided a rRNA genome"
			
		awk 'NR%4==2' "${name}".fastq | sort | uniq -c | sed 's/[[:space:]][[:space:]]//g' \
		| awk '{OFS= "\n"; print ">"NR"-"$1,$2}' > ${name}_collapsed.fasta
		cutadapt -u "$umis" -m "$minlen" -a "$three_adapt" -g "$five_adapt;min_overlap=10" \
		-o pre_rRNA_dep_"${name}"_trimmed.fasta ${name}_collapsed.fasta >> "${name}"_stats.txt	

	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --genomeDir "${rRNA_genome}" \
	--readFilesType Fastx --readFilesIn pre_rRNA_dep_"${name}"_trimmed.fasta --outFileNamePrefix rRNA_aligned_"${name}" \
	--outSAMtype BAM SortedByCoordinate --outSAMattributes NH MD HI AS --runThreadN "$threads" --outReadsUnmapped Fastx --limitBAMsortRAM 2000000000

mv rRNA_aligned_"${name}"Unmapped.out.mate1 rRNA_depleted_"${name}"_trimmed.fasta 2>/dev/null || mv rRNA_aligned_"${name}"Unmapped.out.mate0 rRNA_depleted_"${name}"_trimmed.fasta

  	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --runThreadN "$threads" --genomeDir "${index}" \
	--readFilesType Fastx --readFilesIn rRNA_depleted_"${name}"_trimmed.fasta --outFileNamePrefix "${name}." \
	--outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1 --outMultimapperOrder Random --outSAMattributes NH MD HI AS 
		echo "done!" >> "${name}".log
elif [[ -f ${rRNA_genome}/SAindex && $ext == "fasta" && ${#three_adapt} -eq 0 && $umis -eq 0 ]]
	then
		echo "You have a fasta file and no UMIs and adapters and provided a rRNA genome"
		
	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --genomeDir "${rRNA_genome}" --readFilesIn "${input}" \
	--readFilesType Fastx --outFileNamePrefix rRNA_aligned_"${name}" --outSAMtype BAM SortedByCoordinate \
	--outSAMattributes NH MD HI AS --runThreadN "$threads" --outReadsUnmapped Fastx --limitBAMsortRAM 2000000000
mv rRNA_aligned_"${name}"Unmapped.out.mate1 rRNA_depleted_"${name}"_trimmed.fasta 2>/dev/null || mv rRNA_aligned_"${name}"Unmapped.out.mate0 rRNA_depleted_"${name}"_trimmed.fasta 


	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --runThreadN "$threads" --genomeDir "${index}" \
	--readFilesType Fastx --readFilesIn rRNA_depleted_"${name}"_trimmed.fasta --outFileNamePrefix "${name}." \
	--outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1 --outMultimapperOrder Random --outSAMattributes NH MD HI AS 
		echo "done!" >> "${name}".log
elif [[ ! -f ${rRNA_genome}/SAindex && $ext == "gz" && ${#three_adapt} -eq 0 && $umis -eq 0 ]]
	then
		echo "You have a fastq.gz file and not provie UMIs, adapters and a rRNA genome"
	
		gzip -dc "${name}".fastq.gz | awk 'NR%4==2' | sort \
		| uniq -c | sed 's/[[:space:]][[:space:]]//g' | awk '{OFS= "\n"; print ">"NR"-"$1,$2}' > ${name}_collapsed.fasta
		echo "aligning sequences to genome" >> "${name}".log
  	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --runThreadN "$threads" --genomeDir "${index}" \
	--readFilesType Fastx --readFilesIn ${name}_collapsed.fasta --outFileNamePrefix "${name}." --outSAMtype BAM SortedByCoordinate \
	--outSAMmultNmax 1 --outMultimapperOrder Random --outSAMattributes NH MD HI AS 
		echo "done!" >> "${name}".log

elif [[ -f ${rRNA_genome}/SAindex && $ext == "gz" && ${#three_adapt} -eq 0 && $umis -eq 0 ]]
	then
	echo "You have a fastq.gz file and no UMIs and adapters and provided a rRNA genome"
		gzip -dc "${name}".fastq.gz | awk 'NR%4==2' | sort \
		| uniq -c | sed 's/[[:space:]][[:space:]]//g' | awk '{OFS= "\n"; print ">"NR"-"$1,$2}' > ${name}_collapsed.fasta
		echo "aligning sequences to rRNA_genome" >> "${name}".log
	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --genomeDir "${rRNA_genome}" \
	--readFilesType Fastx --readFilesIn ${name}_collapsed.fasta --outFileNamePrefix rRNA_aligned_"${name}" \
	--outSAMtype BAM SortedByCoordinate --outSAMattributes NH MD HI AS --runThreadN "$threads" --outReadsUnmapped Fastx --limitBAMsortRAM 2000000000
mv rRNA_aligned_"${name}"Unmapped.out.mate1 rRNA_depleted_"${name}"_trimmed.fasta 2>/dev/null || mv rRNA_aligned_"${name}"Unmapped.out.mate0 rRNA_depleted_"${name}"_trimmed.fasta
		echo "aligning sequences to genome" >> "${name}".log
  	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --runThreadN "$threads" --genomeDir "${index}" \
	--readFilesType Fastx --readFilesIn rRNA_depleted_"${name}"_trimmed.fasta --outFileNamePrefix "${name}." \
	--outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1 --outMultimapperOrder Random --outSAMattributes NH MD HI AS 

elif [[ -f ${rRNA_genome}/SAindex && $ext == "gz" && ${#three_adapt} -gt 0 && $umis -eq 0 ]]
	then
				echo "You have a fastq.gz file and have UMIs and adapters and provided a rRNA genome"
		gzip -dc "${name}".fastq.gz | awk 'NR%4==2' | sort | uniq -c \
		| sed 's/[[:space:]][[:space:]]//g' | awk '{OFS= "\n"; print ">"NR"-"$1,$2}' > ${name}_collapsed.fasta
		cutadapt -m "$minlen" -a "$three_adapt" -g "$five_adapt;min_overlap=10" \
		-o pre_rRNA_dep_"${name}"_trimmed.fasta ${name}_collapsed.fasta >> "${name}"_stats.txt

	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --genomeDir "${rRNA_genome}" --readFilesIn pre_rRNA_dep_"${name}"_trimmed.fasta \
	--readFilesType Fastx --outFileNamePrefix rRNA_aligned_"${name}" --outSAMtype BAM SortedByCoordinate --outSAMattributes NH MD HI AS \
	--runThreadN "$threads" --outReadsUnmapped Fastx --limitBAMsortRAM 2000000000
mv rRNA_aligned_"${name}"Unmapped.out.mate1 rRNA_depleted_"${name}"_trimmed.fasta 2>/dev/null || mv rRNA_aligned_"${name}"Unmapped.out.mate0 rRNA_depleted_"${name}"_trimmed.fasta 

  	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --runThreadN "$threads" --genomeDir "${index}" \
	--readFilesType Fastx --readFilesIn rRNA_depleted_"${name}"_trimmed.fasta --outFileNamePrefix "${name}." --outSAMtype BAM SortedByCoordinate \
	--outSAMmultNmax 1 --outMultimapperOrder Random --outSAMattributes NH MD HI AS 
		echo "done!" >> "${name}".log

echo "Could not find any combinations of options for the input"


fi

##############
# index bam ##
##############

samtools index "${name}".Aligned.sortedByCoord.out.bam


###############################
## converting bamtobed    ##
## good for bedtools coverage #
###############################


bedtools bamtobed -i "${name}".Aligned.sortedByCoord.out.bam \
	| awk -F "\t" 'NR==FNR{chr[$1]=1; next} ($1 in chr)' "${name}"_valid_chroms.txt - > "${name}".aligned_sorted.bed

###############################
## get TtoC bam and index it ##
###############################

samtools calmd -e "${name}".Aligned.sortedByCoord.out.bam "${genomefasta}" | awk -F "\t" '(($1 ~"^@") || ($2 == "0" || $2 == "256") && ($10 ~ "C" && $13 ~ "T")) || (($1 ~"^@") || ($2 == "16" || $2 == "272") && ($10 ~ "G" && $13 ~ "A"))' | samtools view -Sbh --threads "$threads" - | samtools sort --threads "$threads" - > "${name}".aligned_TtoC.bam


samtools index "${name}".aligned_TtoC.bam

########################################################################################
### make TtoC bed to run coverage faster (bedtools coverage with bams would collapse ###
########################################################################################


bedtools bamtobed -i "${name}".aligned_TtoC.bam \
	| awk -F "\t" 'NR==FNR{chr[$1]=1; next} ($1 in chr)' "${name}"_valid_chroms.txt - > "${name}".aligned_TtoC.bed


##############################
#### extract TtoC locationg ##
##############################


samtools view --threads "$threads" "${name}".aligned_TtoC.bam \
	| awk -F "\t" 'NR==FNR{chr[$1]=1; next} ($3 in chr)' "${name}"_valid_chroms.txt - \
	| awk '{OFS="\t" } { print $1,$2,$3,$4,$10,$13}' | sed -e 's/MD:Z://g' -e 's/\^//g' -e 's/=*//g' | awk 'BEGIN{FS=OFS="\t"} {gsub("A|G|C|T", "\t&\t", $6)} {gsub("A|G|C|T", "&\t", $5)}1' | awk 'NF >5' | awk '{OFS="\t"} \
  { if ((NF == 8) && ($2 == 0 || $2 == 256) && ($5~"C" && $7~"T")) \
    { print $3,$4+$6-1,$4+$6,$1,".","+"} \
  else if ((NF == 8) && ($2 == 16 || $2 == 272) && ($5~"G" && $7~"A")) \
    { print $3,$4+$6-1,$4+$6,$1,".","-"}
  else if ((NF == 11) && ($2 == 0 || $2 == 256) && ($5~"C" && $8~"T") && ($6~"C" && $10~"T")) \
    { print $3,$4+$7-1,$4+$7,$1,".","+","\n"$3,$4+$7+$9,$4+$7+$9+1,$1,".","+"}
  else if ((NF == 11) && ($2 == 0 || $2 == 256) && ($5~"C" && $8~"T") && ($6!~"C" || $10!~"T")) \
    { print $3,$4+$7-1,$4+$7,$1,".","+"}
  else if ((NF == 11) && ($2 == 0 || $2 == 256) && ($5!~"C" || $8!~"T") && ($6~"C" || $10~"T")) \
    { print $3,$4+$7+$9,$4+$7+$9+1,$1,".","+"}
  else if ((NF == 11) && ($2 == 16 || $2 == 272) && ($5~"G" && $8~"A") && ($6~"G" && $10~"A")) \
    { print $3,$4+$7-1,$4+$7,$1,".","-","\n"$3,$4+$7+$9,$4+$7+$9+1,$1,".","-"}
  else if ((NF == 11) && ($2 == 16 || $2 == 272) && ($5~"G" && $8~"A") && ($6!~"G" || $10!~"A")) \
    { print $3,$4+$7-1,$4+$7,$1,".","-"}
  else if ((NF == 11) && ($2 == 16 || $2 == 272) && ($5!~"G" || $8~"A") && ($6~"G" || $10~"A")) \
    { print $3,$4+$7+$9,$4+$7+$9+1,$1,".","-"}}'  | bedtools sort | bedtools merge -d -1 -s -i - -c 6 -o distinct \
| awk '{OFS="\t"; print $1,$2,$3,"TtoC",".",$4}' | bedtools sort > "${name}".aligned_TtoC_loc.bed


######################
## merging bedfile ##
#####################


## the ($2 != $4 || $3 != $5 ) means that group of reads that have exact same start and end are discarded 
## these are artifacts
#########################################################################################################

bedtools merge -s -i "${name}".aligned_sorted.bed -c 2,3,4,5,6,2 -o mean,mean,count,distinct,distinct,collapse \
	| awk -v var=${minRead} '($6 >= var ) && ($3-$2 > 1) && ($2 != $4 || $3 != $5 )' \
	| sed 's/-[0-9]/1/g' > "${name}".aligned_merged.bed


####################################
##### first round of statistics ####
####################################


cat "${name}".aligned_merged.bed | awk '{OFS="\t"} { gsub(",","\t", $9)}1' | awk '{OFS="\t"} {sum=0; var=0 ; for (i=9; i<=NF; i++) {sum=sum+$i} mean=sum/(NF-8); for (i=9; i<=NF; i++) {var+=($i-mean)^2/(NF-9)} ; sd=sqrt(var) ; print $1,int($4-sd),int($5+sd),$6,"Cln"NR,$8,sd,$5-$4,$2,$3 }' | sed 's/-[0-9]/1/g' | bedtools sort | bedtools merge -s -c 4,5,6,7,8,9,10 -o sum,distinct,distinct,mean,mean,min,max -delim "" > "${name}".aligned_merged_stats.bed 
  
######################################
## coverage on mean or max locations #
######################################

awk '{OFS="\t"} {print $1,$9,$10,$4,$5,$6,$7,$8}' "${name}".aligned_merged_stats.bed \
	| bedtools coverage -counts -s -a - -b "${name}".aligned_TtoC_loc.bed \
	| bedtools coverage -counts -s -F 0.66 -a - -b "${name}".aligned_TtoC.bed \
	| bedtools coverage -counts -s -F 0.66 -a - -b "${name}".aligned_sorted.bed > "${name}".MinMax_Cov.bed 


awk '{OFS="\t"} { print $1,$2,$3,$4,$5,$6,$7,$8}' "${name}".aligned_merged_stats.bed \
	| bedtools coverage -counts -s -a - -b "${name}".aligned_TtoC_loc.bed \
	| bedtools coverage -counts -s -F 0.66 -a - -b "${name}".aligned_TtoC.bed \
	| bedtools coverage -counts -s -F 0.66 -a - -b "${name}".aligned_sorted.bed > "${name}".MeanStartMeanEnd_Cov.bed



##############################################################################
# Filter by distribution SDV of reads set > 68% for normal at sd-1 - sd +1g 
# Peaks that either will not follow normal distribution or the SDV is as
##############################################################################


paste "${name}".MinMax_Cov.bed "${name}".MeanStartMeanEnd_Cov.bed | awk -F "\t" '$11 > 0' > "${name}"_pasted.txt


awk -F "\t" '(($22/$11 >= 0.68) && ($7 < $8))' "${name}"_pasted.txt | awk -F "\t" '{OFS="\t"} \
  { if (($2 < $13) && ($3 > $14)) \
    { print $1,$13,$14,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 > $13) && ($3 > $14)) \
    { print $1,$2,$14,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 > $13) && ($3 < $14)) \
    { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 < $13) && ($3 < $14)) \
    { print $1,$13,$3,$4,$5,$6,$7,$8,$9,$10,$11}
  else if (($2 == $13) && ($3 > $14)) \
    { print $1,$13,$14,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 > $13) && ($3 == $14)) \
    { print $1,$2,$14,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 == $13) && ($3 < $14)) \
    { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 < $13) && ($3 == $14)) \
    { print $1,$13,$3,$4,$5,$6,$7,$8,$9,$10,$11} \
   else { print $1,$13,$3,$4,$5,$6,$7,$8,$9,$10,$11} }' > "${name}".groups_normal.bed
    
awk '{OFS="\t"} {if (($22/$11 < 0.68) || ($7 > $8)) print $1,$2,$3,sqrt($4)/$4,$5,$6,$7,$8,$9,$10,$11}' ${name}_pasted.txt > "${name}".groups_NonNormal.bed 


rm -f ${name}_pasted.txt


#############################################
# subsample bam file for NonNormal clusters #
#############################################

number_of_NonNormal=`wc -l < "${name}".groups_NonNormal.bed`


samtools view -H "${name}".Aligned.sortedByCoord.out.bam > "${name}".aligned_sqrt.sam


bedtools intersect -abam "${name}".Aligned.sortedByCoord.out.bam -b "${name}".groups_NonNormal.bed > "${name}".aligned_NonNormal.bam

samtools index "${name}".aligned_NonNormal.bam


# --------------------------------------------------------------------------
# Subsampling used to loop over groups_NonNormal.bed, spawning one separate
# samtools process per region (slow at scale, and reaching for GNU parallel
# to fix that would add a dependency that isn't always available/installable
# without sudo). Replaced with one pass instead: bedtools intersect pairs
# every read with its region's subsampling fraction in a single call, one
# awk pass makes the keep/drop draw per read, and one final samtools call
# extracts exactly the kept read names - no loop, no extra dependency
# (bedtools/samtools/awk are already required either way).
#
# The keep/drop decision is a deterministic hash of the read name, NOT
# awk's rand(): tested directly, awk's own rand()/srand() turned out to be
# non-reproducible across separate process invocations in this environment
# (mawk), even with an explicit seed - so two runs on identical input could
# silently subsample different reads. A hash avoids depending on any RNG
# implementation at all. A naive polynomial hash isn't good enough either -
# read names in this pipeline are frequently sequential, and a weak hash
# puts sequential names at near-identical values, biasing the "random"
# subsample toward whole contiguous blocks instead of a genuine scatter.
# The extra multiplicative mixing step below fixes that - verified uniform
# across 5000 sequential test names (490-519 per decile, expected ~500) and
# reproducible across repeated runs on the same input (confirmed by md5sum
# on the exact assembled block, not just an isolated snippet).
# --------------------------------------------------------------------------
stage_start_ts=$(date +%s)
bedtools intersect -abam "${name}".aligned_NonNormal.bam -b "${name}".groups_NonNormal.bed -bed -wa -wb \
	| awk -F "\t" '
function hash(s,   h, i, n, c) {
	n = length(s)
	h = 5381
	for (i = 1; i <= n; i++) {
		c = index("qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM0123456789_-.: ", substr(s,i,1))
		h = (h * 33 + c) % 4294967291
	}
	h = (h * 2654435761) % 4294967291
	return h
}
{ if ((hash($4) % 100000) / 100000.0 <= $16) print $4 }' \
	| sort -u > "${name}"_sqrt_kept_names.txt

samtools view -N "${name}"_sqrt_kept_names.txt "${name}".aligned_NonNormal.bam >> "${name}".aligned_sqrt.sam

rm -f "${name}"_sqrt_kept_names.txt
echo "timing: NonNormal subsampling took $(( $(date +%s) - stage_start_ts ))s (${number_of_NonNormal} regions)" >> "${name}"_stats.txt


samtools view -Shb "${name}".aligned_sqrt.sam --threads ${threads} \
	| samtools sort - --threads ${threads} | bedtools bamtobed -i - \
	| bedtools merge -s -i - -c 2,3,4,5,6,2 -o mean,mean,count,distinct,distinct,collapse \
	| awk -v var=${minRead} '($6 >= var ) && ($3-$2 > 1)' | sed 's/-[0-9]/1/g' > "${name}".aligned_merged_sqrt.bed



######################################################
# second round on subsamples sam round of statistics #
######################################################

cat "${name}".aligned_merged_sqrt.bed | awk '{OFS="\t"} { gsub(",","\t", $9)}1' \
	| awk '{OFS="\t"} {sum=0; var=0 ; for (i=9; i<=NF; i++) {sum=sum+$i} mean=sum/(NF-8); for (i=9; i<=NF; i++) {var+=($i-mean)^2/(NF-9)} ; sd=sqrt(var) ; print $1,int($4-sd),int($5+sd),$6,"Cla"NR,$8,sd,$5-$4,$2,$3 }' \
	| sed 's/-[0-9]/1/g' | bedtools sort | bedtools merge -s -c 4,5,6,7,8,9,10 -o sum,distinct,distinct,mean,mean,min,max -delim "" > "${name}".aligned_merged_stats_sqrt.bed 
  
awk '{OFS="\t"} {print $1,$9,$10,$4,$5,$6,$7,$8}' "${name}".aligned_merged_stats_sqrt.bed \
	| bedtools coverage -counts -s -a - -b "${name}".aligned_TtoC_loc.bed \
	| bedtools coverage -counts -s -F 0.66 -a - -b "${name}".aligned_TtoC.bed | bedtools coverage -counts -s -F 0.66 -a - -b "${name}".aligned_sorted.bed > "${name}".MinMax_Cov_sqrt.bed 

awk '{OFS="\t"} { print $1,$2,$3,$4,$5,$6,$7,$8}' "${name}".aligned_merged_stats_sqrt.bed \
	| bedtools coverage -counts -s -a - -b "${name}".aligned_TtoC_loc.bed \
	| bedtools coverage -counts -s -F 0.66 -a - -b "${name}".aligned_TtoC.bed \
	| bedtools coverage -counts -s -F 0.66 -a - -b "${name}".aligned_sorted.bed > "${name}".MeanStartMeanEnd_Cov_sqrt.bed

paste "${name}".MinMax_Cov_sqrt.bed "${name}".MeanStartMeanEnd_Cov_sqrt.bed | awk -F "\t" '$11 > 0' > "${name}"_pasted_sqrt.txt


awk -F "\t" '(($22/$11 >= 0.68) && ($7 < $8))' "${name}"_pasted_sqrt.txt |awk -F "\t" '{OFS="\t"} \
   { if (($2 < $13) && ($3 > $14)) \
    { print $1,$13,$14,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 > $13) && ($3 > $14)) \
    { print $1,$2,$14,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 > $13) && ($3 < $14)) \
    { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 < $13) && ($3 < $14)) \
    { print $1,$13,$3,$4,$5,$6,$7,$8,$9,$10,$11}
  else if (($2 == $13) && ($3 > $14)) \
    { print $1,$13,$14,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 > $13) && ($3 == $14)) \
    { print $1,$2,$14,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 == $13) && ($3 < $14)) \
    { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 < $13) && ($3 == $14)) \
    { print $1,$13,$3,$4,$5,$6,$7,$8,$9,$10,$11} \
   else { print $1,$13,$3,$4,$5,$6,$7,$8,$9,$10,$11} }' > "${name}".groups_normal_sqrt.bed
   
#################################################################################################################
# Merge normal stage 1 and normal stage 2 and filter satellites by tukeys factor is set by nt from peak to peak #
#################################################################################################################

factor=50

cat "${name}".groups_normal.bed "${name}".groups_normal_sqrt.bed | bedtools sort | awk -F "\t" -v var=$factor '{OFS="\t"; print $1,$2,$3,$5,".",$6,$7,$8,$9,$10,$11}' > "${name}"_filt.groups
 
 
 
###########################################
# annotate peaks using bedtools intersect #
###########################################

# --------------------------------------------------------------------------
# v1.6: this whole block was rebuilt. The old approach (-wa -wb -loj, then
# bedtools merge -o distinct independently on the gene/annotation/biotype
# columns) resolved each of those three columns separately, with no
# guarantee they stayed matched to the same underlying feature. A cluster
# hitting GENE_A (exon, protein_coding) and GENE_B (intron, tRNA) could end
# up reported as GENE_B + exon + tRNA - a combination that matches neither
# real feature. It also never captured overlap length, so ties between two
# same-tier hits were broken arbitrarily (whichever bedtools happened to
# list first), not by which one the cluster actually mostly sits on.
#
# Fixed by: capturing overlap-bp per hit (-wo), packing gene:annotation:
# biotype:overlap_bp into ONE field per hit so it survives the merge intact,
# then resolving in a single pass that picks one whole winning hit - ranked
# by the same biotype/annotation priority tiers as before, with overlap-bp
# only as a tiebreaker within a tied tier. (-wo and -loj can't be combined
# in this bedtools version - -loj's no-hit behavior is replicated by hand
# via a separate -v pass for clusters with zero overlaps.)
# --------------------------------------------------------------------------

bedtools intersect -s -wo -a "${name}"_filt.groups -b "${annotation}" > "${name}"_hits.tmp
bedtools intersect -s -v -a "${name}"_filt.groups -b "${annotation}" \
	| awk -F "\t" '{OFS="\t"; print $0,".",".",".",".",".",".",".",0}' > "${name}"_nohits.tmp

cat "${name}"_hits.tmp "${name}"_nohits.tmp | bedtools sort \
	| awk -F "\t" '{OFS="\t"; unit=$15":"$16":"$18":"$19; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,($10/$11),unit}' \
	| bedtools sort \
	| bedtools merge -s -c 4,5,6,7,8,9,10,11,12,13 \
	-o distinct,distinct,distinct,max,max,max,max,max,max,distinct -delim "|" > "${name}".groups.tsv

rm -f "${name}"_hits.tmp "${name}"_nohits.tmp

cat "${name}".groups.tsv | awk -F "\t" '{OFS="\t"; print $1,$2,$3,$4,$5,$6 }' \
	| bedtools getfasta -s -fi ${genomefasta} -bed - -tab > "${name}".groups.fa


echo "chromosome,start,end,cluster_id,dot,strand,gene_name,annotation,gene_type,sequence,crosslinked_reads,total_reads,fraction,TtoC_locations,avr_readlength,read_distance_SD,chr:start-end(strand)" > "${name}".groups.csv

paste "${name}".groups.tsv "${name}".groups.fa | awk -F "\t" 'BEGIN{
  OFS=",";
  priority_annot["exon"]=1; priority_annot["3ss"]=2; priority_annot["5ss"]=3;
  priority_annot["start_codon"]=4; priority_annot["stop_codon"]=5;
  priority_annot["3utr"]=6; priority_annot["5utr"]=7;
  priority_annot["CDS"]=8; priority_annot["intron"]=9;
  priority_biotype["snoRNA"]=1; priority_biotype["snRNA"]=2; priority_biotype["rRNA"]=3;
  priority_biotype["miRNA"]=4; priority_biotype["tRNA"]=5; priority_biotype["Mt_tRNA"]=6;
  priority_biotype["scaRNA"]=7; priority_biotype["protein_coding"]=8; priority_biotype["lncRNA"]=9;
}
{
  # $1-6 core, $7-11 numeric data, $12 fraction, $13 packed unit-list, $14 location, $15 sequence
  n = split($13, units, "|");
  best_bt=999; best_at=999; best_ov=-1;
  best_gene="NA"; best_annot="No_Annotation"; best_biotype="NA";
  for (u=1; u<=n; u++) {
    split(units[u], f, ":");
    gene=f[1]; annot=f[2]; biotype=f[3]; ov=f[4]+0;
    if (gene == ".") continue;
    bt = (biotype in priority_biotype) ? priority_biotype[biotype] : 99;
    at = (annot in priority_annot) ? priority_annot[annot] : 99;
    if (bt < best_bt || (bt==best_bt && at<best_at) || (bt==best_bt && at==best_at && ov>best_ov)) {
      best_bt=bt; best_at=at; best_ov=ov; best_gene=gene; best_annot=annot; best_biotype=biotype;
    }
  }
  print $1,$2,$3,$4,$5,$6,best_gene,best_annot,best_biotype,$15,$10,$11,$12,$9,$8,$7,$14
}' >> "${name}".groups.csv

cat "${name}".groups.csv | awk -F "," -v convloc="${convloc}" -v frac="${fraction}" '{ if (NR==1 || ($14 >= convloc && $13 >=frac)) print $0 } ' > "${name}".clusters.csv


grep -v chromosome "${name}".clusters.csv \
	| awk -F "," '{OFS="\t"; print $1,$2,$3,$5,$4,$6 }' \
	| bedtools sort > "${name}".clusters.bed


echo "Binding statistics" >> "${name}"_stats.txt

groups=$(wc -l < "${name}".groups.csv)
echo "groups: ${groups}" >> "${name}"_stats.txt
clusters=$(wc -l < "${name}".clusters.csv)
echo "clusters: ${clusters}" >> "${name}"_stats.txt

Protein_coding_5UTR_XL_reads=$(grep 5utr < "${name}".clusters.csv | awk -F',' '{sum+=$11;}END{print sum;}')
echo "Protein_coding-5UTR_XL_reads: ${Protein_coding_5UTR_XL_reads}" >> "${name}"_stats.txt
Protein_coding_CDS_XL_reads=$(grep CDS < "${name}".clusters.csv | awk -F',' '{sum+=$11;}END{print sum;}')
echo "Protein_coding-CDS_XL_reads: ${Protein_coding_CDS_XL_reads}" >> "${name}"_stats.txt
Protein_coding_intron_XL_reads=$(grep intron < "${name}".clusters.csv | awk -F',' '{sum+=$11;}END{print sum;}')
echo "Protein_coding-introns_XL_reads: ${Protein_coding_intron_XL_reads}" >> "${name}"_stats.txt
Protein_coding_3UTR_XL_reads=$(grep 3utr < "${name}".clusters.csv | awk -F',' '{sum+=$11;}END{print sum;}')
echo "Protein_coding-3UTR_XL_reads: ${Protein_coding_3UTR_XL_reads}" >> "${name}"_stats.txt
lncRNA_XL_reads=$(grep lncRNA < "${name}".clusters.csv | awk -F',' '{sum+=$11;}END{print sum;}')
echo "lncRNA_XL_reads: ${lncRNA_XL_reads}" >> "${name}"_stats.txt
snoRNA_XL_reads=$(grep snoRNA < "${name}".clusters.csv | awk -F',' '{sum+=$11;}END{print sum;}')
echo "snoRNA_XL_reads: ${snoRNA_XL_reads}" >> "${name}"_stats.txt
snRNA_XL_reads=$(grep snRNA < "${name}".clusters.csv | awk -F',' '{sum+=$11;}END{print sum;}')
echo "snRNA_XL_reads: ${snRNA_XL_reads}" >> "${name}"_stats.txt
rRNA_XL_reads=$(grep rRNA < "${name}".clusters.csv | awk -F',' '{sum+=$11;}END{print sum;}')
echo "rRNA_XL_reads: ${rRNA_XL_reads}" >> "${name}"_stats.txt

Protein_coding_5UTR_clusters=$(grep 5utr < "${name}".clusters.csv | wc -l)
echo "Protein_coding-5UTR_clusters: ${Protein_coding_5UTR_clusters}" >> "${name}"_stats.txt
Protein_coding_CDS_clusters=$(grep CDS < "${name}".clusters.csv | wc -l)
echo "Protein_coding-CDS_clusters: ${Protein_coding_CDS_clusters}" >> "${name}"_stats.txt
Protein_coding_intron_clusters=$(grep intron < "${name}".clusters.csv | wc -l)
echo "Protein_coding-introns_clusters: ${Protein_coding_intron_clusters}" >> "${name}"_stats.txt
Protein_coding_3UTR_clusters=$(grep 3utr < "${name}".clusters.csv | wc -l)
echo "Protein_coding-3UTR_clusters: ${Protein_coding_3UTR_clusters}" >> "${name}"_stats.txt
lncRNA_clusters=$(grep lncRNA < "${name}".clusters.csv | wc -l)
echo "lncRNA_clusters: ${lncRNA_clusters}" >> "${name}"_stats.txt
snoRNA_clusters=$(grep snoRNA < "${name}".clusters.csv | wc -l)
echo "snoRNA_clusters: ${snoRNA_clusters}" >> "${name}"_stats.txt
snRNA_clusters=$(grep snRNA < "${name}".clusters.csv | wc -l)
echo "snRNA_clusters: ${snRNA_clusters}" >> "${name}"_stats.txt
rRNA_clusters=$(grep rRNA < "${name}".clusters.csv | wc -l)
echo "rRNA_clusters: ${rRNA_clusters}" >> "${name}"_stats.txt



###################
## cleaning files #
###################

rm -f "${name}"_collapsed.fasta
#rm "${name}"_trimmed.fasta
rm -f "${name}".SJ.out.tab
rm -f "${name}".Log.progress.out
rm -f "${name}".Log.final.out
rm -f "${name}".Log.out
rm -f "${name}".aligned_sorted.bed
rm -f "${name}".aligned_merged.bed
rm -f "${name}".aligned_merged_stats.bed
rm -f "${name}".MinMax_Cov.bed
rm -f "${name}".MeanStartMeanEnd_Cov.bed
rm -f "${name}".groups_normal.bed
rm -f "${name}".aligned_sqrt.sam
rm -f "${name}".aligned_merged_sqrt.bed
rm -f "${name}".aligned_merged_stats_sqrt.bed
rm -f "${name}".MinMax_Cov_sqrt.bed
rm -f "${name}".MeanStartMeanEnd_Cov_sqrt.bed
rm -f "${name}"_pasted_sqrt.txt
rm -f "${name}".groups_normal_sqrt.bed
rm -f "${name}"_filt.groups
rm -f "${name}".log
rm -f "${name}".groups.tsv
rm -f "${name}".groups.fa
rm -f "${name}"_valid_chroms.txt
rm -f "${name}".aligned_NonNormal.bam
rm -f "${name}".aligned_NonNormal.bam.bai
rm -f "${name}".groups_NonNormal.bed 

# --------------------------------------------------------------------------
# rRNA pre-mapping intermediates (only created if -e was used): the
# trimmed-but-not-yet-rRNA-checked fasta, the rRNA-depleted fasta fed into
# the real genome alignment, and STAR's own small log files for that run -
# all fully consumed by this point, not needed again. Deliberately NOT
# deleting rRNA_aligned_"${name}"Aligned.sortedByCoord.out.bam - that's the
# actual reads that mapped TO rRNA, which is real, useful data if you want
# to look at binding on rRNA itself, not just an intermediate.
# --------------------------------------------------------------------------
rm -f pre_rRNA_dep_"${name}"_trimmed.fasta
rm -f rRNA_depleted_"${name}"_trimmed.fasta
rm -f rRNA_aligned_"${name}"SJ.out.tab
rm -f rRNA_aligned_"${name}"Log.progress.out
rm -f rRNA_aligned_"${name}"Log.out
rm -f rRNA_aligned_"${name}"Log.final.out

# --------------------------------------------------------------------------
# aligned_TtoC_loc.bed (deduped, single row per unique crosslink location) -
# fully consumed by the coverage steps earlier in the pipeline, and now
# reconstructable on demand from aligned_TtoC.bed (the per-read version,
# kept since 1.6.21) via the same merge command used for the metagene
# weighting: no need to keep both.
# --------------------------------------------------------------------------
rm -f "${name}".aligned_TtoC_loc.bed


###########################
# get the 5mer enrichment #
###########################

# --------------------------------------------------------------------------
# Background is now both an upstream (-200) and downstream (+200) flanking
# region per cluster, not just one side - kept deliberately simple (flat
# genomic shift, not strand-aware/gene-constrained like an earlier, more
# complex attempt at this that kept surfacing hard-to-reproduce bugs):
#   - a side is dropped if the shift would go past either end of the
#     chromosome (checked against the genome fasta's own .fai index)
#   - a side is dropped if the shifted region lands on any real cluster
#     (bedtools subtract -A against clusters.bed)
# Each side is dropped independently - a cluster can lose one side and
# still keep the other. The final background count per k-mer is the
# average of its upstream and downstream counts, rounded to the nearest
# integer (standard round-half-up: int(avg+0.5)).
# --------------------------------------------------------------------------

if [ ! -f "${genomefasta}.fai" ]; then samtools faidx "${genomefasta}"; fi

awk -F "\t" '{OFS="\t"; print $1,$2-200,$3-200,$4,$5,$6}' "${name}".clusters.bed \
	| awk -F "\t" '$2 >= 0' > "${name}"_bg_upstream_candidate.bed

awk -F "\t" '{OFS="\t"; print $1,$2+200,$3+200,$4,$5,$6}' "${name}".clusters.bed > "${name}"_bg_downstream_raw.bed
awk -F "\t" 'NR==FNR{len[$1]=$2; next} { if ($3 <= len[$1]) print }' "${genomefasta}.fai" "${name}"_bg_downstream_raw.bed \
	> "${name}"_bg_downstream_candidate.bed
rm -f "${name}"_bg_downstream_raw.bed

bedtools subtract -s -A -a "${name}"_bg_upstream_candidate.bed -b "${name}".clusters.bed > background_upstream_"${name}".clusters.bed
bedtools subtract -s -A -a "${name}"_bg_downstream_candidate.bed -b "${name}".clusters.bed > background_downstream_"${name}".clusters.bed
rm -f "${name}"_bg_upstream_candidate.bed "${name}"_bg_downstream_candidate.bed

up_n=$(wc -l < background_upstream_"${name}".clusters.bed)
down_n=$(wc -l < background_downstream_"${name}".clusters.bed)
echo "5mer enrichment: ${up_n} upstream / ${down_n} downstream background regions (each excludes chromosome-edge and other-cluster overlaps independently)" >> "${name}"_stats.txt

bedtools getfasta -s -fi "${genomefasta}" -bed "${name}".clusters.bed -fo bound_"${name}".fasta
bedtools getfasta -s -fi "${genomefasta}" -bed background_upstream_"${name}".clusters.bed -fo background_upstream_"${name}".fasta
bedtools getfasta -s -fi "${genomefasta}" -bed background_downstream_"${name}".clusters.bed -fo background_downstream_"${name}".fasta

rm -f background_upstream_"${name}".clusters.bed background_downstream_"${name}".clusters.bed

##############################################
#create all combinations of 5mers and count
##############################################

echo {A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G} | tr ' ' '\n' > "${name}"_5mers.txt


#count each mer in the fasta files

grep -v ">" background_upstream_"${name}".fasta | grep -v "@" > filt_background_upstream_"${name}".fasta
grep -v ">" background_downstream_"${name}".fasta | grep -v "@" > filt_background_downstream_"${name}".fasta
grep -v ">" bound_"${name}".fasta | grep -v "@" > filt_bound_"${name}".fasta

rm -f background_upstream_"${name}".fasta background_downstream_"${name}".fasta bound_"${name}".fasta



# --------------------------------------------------------------------------
# This used to loop over all 1024 mers, spawning a separate grep process per
# mer per fasta file (2048 subprocess spawns total). Replaced with one awk
# pass per fasta file - but this needs care to actually match grep -o's
# semantics, which an earlier, insufficiently-tested attempt at this same
# fix got wrong in two ways (caught by comparing real output, not assumed
# from a 3-4-mer toy test):
#
#   1) grep -o gives each of the 1024 mers its OWN independent scan of the
#      sequence - a single shared "skip 5 ahead after any match" pointer is
#      NOT equivalent, since two different mers can legitimately match at
#      overlapping starting positions (e.g. "AAAATTT" contains AAAAT at
#      position 1 AND AAATT at position 2 AND AATTT at position 3 - grep
#      finds all three independently; a shared pointer finds only the
#      first and skips past the other two entirely). Fixed by tracking a
#      separate "next allowed position" per mer, not one shared pointer.
#
#   2) that per-mer tracking must reset at the start of every sequence
#      (fasta line) - grep's matching is inherently line-scoped (a match
#      can't span a newline), so a mer's match near the end of one
#      sequence must not block a legitimate match near the start of the
#      next, unrelated sequence.
#
# Verified byte-identical against the real grep -o loop at full 1024-mer
# scale on a multi-sequence test file (not just a handful of hand-picked
# mers) before shipping this.
# --------------------------------------------------------------------------
count_kmers_verified () {
	local fasta="$1"
	awk -v merlist="${name}"_5mers.txt 'BEGIN{
		while ((getline m < merlist) > 0) counts[m]=0
		close(merlist)
	}
	{
		delete next_allowed
		seq=toupper($0); L=length(seq);
		for (i=1; i<=L-5+1; i++) {
			w=substr(seq,i,5);
			if (w in counts) { if (i >= next_allowed[w]) { counts[w]++; next_allowed[w]=i+5 } }
		}
	}
	END{ while ((getline m < merlist) > 0) print m"\t"counts[m]; close(merlist) }' "${fasta}"
}

count_kmers_verified filt_background_upstream_"${name}".fasta | sort > background_upstream_"${name}"_5mers_out1.txt
count_kmers_verified filt_background_downstream_"${name}".fasta | sort > background_downstream_"${name}"_5mers_out1.txt
count_kmers_verified filt_bound_"${name}".fasta | sort > bound_"${name}"_5mers_out1.txt

rm -f filt_background_upstream_"${name}".fasta filt_background_downstream_"${name}".fasta filt_bound_"${name}".fasta


paste background_upstream_"${name}"_5mers_out1.txt background_downstream_"${name}"_5mers_out1.txt bound_"${name}"_5mers_out1.txt \
	| awk -F "\t" '{OFS="\t"; avg=($2+$4)/2; print $1,int(avg+0.5),$6}' > "${name}"_5mers.txt


rm -f background_upstream_"${name}"_5mers_out1.txt background_downstream_"${name}"_5mers_out1.txt bound_"${name}"_5mers_out1.txt


##########################################################
# k-mer enrichment ranking + plot (background vs bound)  #
##########################################################

if command -v python3 >/dev/null 2>&1 && python3 -c "import matplotlib, numpy" >/dev/null 2>&1; then

cat > "${name}"_kmer_enrichment.py << 'PYEOF'
#!/usr/bin/env python3
"""
kmer_enrichment.py - minimal k-mer count comparison tool for PARALUS
Input: whitespace-delimited file, 3 columns: kmer  count1  count2
Output: TSV to stdout: kmer  count1  count2  ratio(c2/c1), sorted by ratio.
Also prints top/bottom N to stderr as a quick summary.
"""
import sys
import argparse
import numpy as np


def load(path):
    kmers, c1, c2 = [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            k, a, b = line.split()
            kmers.append(k)
            c1.append(int(a))
            c2.append(int(b))
    return np.array(kmers), np.array(c1), np.array(c2)


def main():
    ap = argparse.ArgumentParser(description="k-mer count enrichment comparison")
    ap.add_argument("infile")
    ap.add_argument("--top", type=int, default=32)
    ap.add_argument("--plot", metavar="PNG")
    ap.add_argument("--pseudocount", type=float, default=1.0,
                     help="added to both counts before the ratio - without this, any k-mer "
                          "with zero background produces an infinite/undefined ratio, which "
                          "real sparse k-mer data hits constantly. A single random hit against "
                          "zero background isn't real signal; the pseudocount stops low-count "
                          "k-mers from dominating the top list purely from sampling noise")
    args = ap.parse_args()

    kmers, c1, c2 = load(args.infile)
    if len(kmers) == 0:
        print("no data in input file", file=sys.stderr)
        return

    p = args.pseudocount
    ratio = (c2 + p) / (c1 + p)
    order = np.argsort(ratio)

    for i in order:
        print(f"{kmers[i]}\t{c1[i]}\t{c2[i]}\t{ratio[i]:.4f}")

    n = min(args.top, len(kmers))
    print(f"\n# Top {n} de-enriched in col2:", file=sys.stderr)
    for i in order[:n]:
        print(f"{kmers[i]}\tc1={c1[i]}\tc2={c2[i]}\tratio={ratio[i]:.3f}", file=sys.stderr)
    print(f"\n# Top {n} enriched in col2:", file=sys.stderr)
    for i in order[::-1][:n]:
        print(f"{kmers[i]}\tc1={c1[i]}\tc2={c2[i]}\tratio={ratio[i]:.3f}", file=sys.stderr)

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        top_idx = order[::-1][:n]
        bot_idx = order[:n]
        fig, ax = plt.subplots(figsize=(7, 7))
        ax.scatter(c1, c2, s=8, alpha=0.35, color="#2c7fb8")
        lims = [0, max(c1.max(), c2.max()) * 1.02] if max(c1.max(), c2.max()) > 0 else [0, 1]
        ax.plot(lims, lims, "r--", linewidth=1, label="y=x")
        for i in top_idx:
            ax.scatter(c1[i], c2[i], s=30, color="red", zorder=5)
            ax.annotate(kmers[i], (c1[i], c2[i]), textcoords="offset points", xytext=(5, 5), fontsize=8, color="red")
        for i in bot_idx:
            ax.scatter(c1[i], c2[i], s=30, color="blue", zorder=5)
            ax.annotate(kmers[i], (c1[i], c2[i]), textcoords="offset points", xytext=(5, -10), fontsize=8, color="blue")
        ax.set_xlabel("background count")
        ax.set_ylabel("bound count")
        ax.set_title(f"k-mer counts (top {n} enriched=red, de-enriched=blue)")
        ax.legend()
        plt.tight_layout()
        plt.savefig(args.plot, dpi=150)


if __name__ == "__main__":
    main()
PYEOF

python3 "${name}"_kmer_enrichment.py "${name}"_5mers.txt --plot "${name}"_5mer_enrichment.pdf --top 32 \
	> "${name}"_5mer_ratios.tsv 2>> "${name}"_stats.txt

rm -f "${name}"_kmer_enrichment.py
fi

##########################################################
# metagene analysis: 5'UTR/CDS/3'UTR and intron, run      #
# against clusters (unweighted), groups (unweighted), and #
# TtoC locations (weighted by real per-location read      #
# depth, recovered from the per-read aligned_TtoC.bed that #
# 1.6.21 stopped deleting specifically for this purpose)   #
##########################################################

if command -v python3 >/dev/null 2>&1 && python3 -c "import matplotlib, numpy" >/dev/null 2>&1; then

# --------------------------------------------------------------------------
# TtoC locations, weighted by real read depth per merged location (a
# location seen in 7 reads counts more than one seen in a single read) -
# recovered from the per-read aligned_TtoC.bed that 1.6.21 stopped
# deleting specifically for this.
# --------------------------------------------------------------------------
sort -k1,1 -k2,2n "${name}".aligned_TtoC.bed | bedtools merge -d -1 -s -i - -c 6,6 -o distinct,count \
	| awk '{OFS="\t"; print $1,$2,$3,"TtoC",$5,$4}' > "${name}"_metagene_ttoc.bed

cat > "${name}"_metagene_utr_cds.py << 'METAGENE_UTRCDS_PYEOF'
#!/usr/bin/env python3
"""
metagene_utr_cds.py - 5'UTR -> CDS -> 3'UTR composite metagene for PARALUS

Uses Trierarch's CANONICAL annotation bed (built with -C) specifically
because a gene's 5'UTR/CDS/3'UTR boundaries are only well-defined for ONE
transcript - without -C, different isoforms of the same gene can have
different UTR/CDS boundaries, and "this gene's 5'UTR is 340nt" stops being
a single well-defined number.

Algorithm:
  1. Parse the annotation bed, keep protein_coding 5utr/CDS/3utr rows.
  2. Per gene, concatenate each region's exonic pieces in TRANSCRIPT order
     (ascending genomic coordinate for + strand, descending for -), giving
     each gene a total length per region and a genomic->relative-position
     mapping. Only genes with all three regions present are used.
  3. Map each input site (a bed file of crosslink sites/clusters) onto this
     mapping via bedtools intersect, converting its genomic position to a
     0-1 relative position within whichever region it falls in.
  4. Bin into a fixed number of bins per region. Each gene's profile is
     normalized to sum to 1 across all its bins (so one highly-bound gene
     can't dominate the aggregate), then averaged across all qualifying
     genes (with a minimum total-site threshold to reduce noise from
     genes with only 1-2 sites).
  5. Plot as one continuous track (5'UTR | CDS | 3'UTR) with vertical
     lines at the two region boundaries.

Optional confidence weighting: column 5 (the BED "score" field) of the
sites file, if numeric, is used as a per-site weight instead of counting
every site equally - so a cluster backed by 200 reads counts more than one
barely passing threshold with 2. A "." or missing column 5 defaults to
weight 1 (fully backward compatible with a plain 6-column sites bed).

Usage:
    metagene_utr_cds.py annotation.bed sites.bed --plot out.png
"""
import sys
import argparse
import subprocess
import tempfile
import os
from collections import defaultdict

import numpy as np


def load_annotation(path):
    """Returns {gene: {"strand": s, "5utr": [(s,e),...], "CDS": [...], "3utr": [...]}}"""
    genes = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            chrom, start, end, gene, annot, strand, biotype = parts[:7]
            if biotype != "protein_coding":
                continue
            if annot not in ("5utr", "CDS", "3utr"):
                continue
            g = genes.setdefault(gene, {"strand": strand, "chrom": chrom,
                                         "5utr": [], "CDS": [], "3utr": []})
            g[annot].append((int(start), int(end)))
    return genes


def build_gene_mapping(genes, min_region_len=1):
    """
    For each gene with all three regions present, build:
      - per-region transcript-order segment list with cumulative mapping
      - per-region total length
    Returns {gene: {"strand":.., "chrom":.., region: {"segments": [...], "total": L}}}
    """
    mapped = {}
    for gene, g in genes.items():
        if not g["5utr"] or not g["CDS"] or not g["3utr"]:
            continue
        strand = g["strand"]
        entry = {"strand": strand, "chrom": g["chrom"]}
        ok = True
        for region in ("5utr", "CDS", "3utr"):
            segs = g[region]
            if strand == "+":
                ordered = sorted(segs, key=lambda x: x[0])
            else:
                ordered = sorted(segs, key=lambda x: -x[1])
            segmap = []
            cum = 0
            for (s, e) in ordered:
                length = e - s
                segmap.append((s, e, cum, cum + length))
                cum += length
            if cum < min_region_len:
                ok = False
                break
            entry[region] = {"segments": segmap, "total": cum}
        if ok:
            mapped[gene] = entry
    return mapped


def genomic_to_relative(pos, segmap):
    """pos: 1-based-ish genomic coordinate (bed start of the site).
    segmap: list of (genomic_start, genomic_end, cum_start, cum_end).
    Returns cumulative position within the region, or None if not found."""
    for (s, e, cs, ce) in segmap["segments"]:
        if s <= pos < e:
            offset = pos - s
            return cs + offset
    return None


def run_bedtools_intersect(sites_path, region_bed_path):
    """Runs bedtools intersect -wo, returns list of (site_chrom, site_start, site_end,
    site_name, ann_gene, ann_annot, ann_start_genomic)."""
    cmd = ["bedtools", "intersect", "-wo", "-a", sites_path, "-b", region_bed_path]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True)
    rows = []
    for line in out.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 13:
            continue
        rows.append(parts)
    return rows


def main():
    ap = argparse.ArgumentParser(description="5'UTR/CDS/3'UTR composite metagene")
    ap.add_argument("annotation", help="Trierarch canonical (-C) annotation bed")
    ap.add_argument("sites", help="bed file of crosslink sites or clusters")
    ap.add_argument("--plot", metavar="PNG", required=True)
    ap.add_argument("--bins-5utr", type=int, default=50)
    ap.add_argument("--bins-cds", type=int, default=100)
    ap.add_argument("--bins-3utr", type=int, default=50)
    ap.add_argument("--min-sites-per-gene", type=int, default=3,
                     help="genes with fewer than this many total sites are excluded, "
                          "to reduce noise from genes with only 1-2 sites")
    args = ap.parse_args()

    genes = load_annotation(args.annotation)
    mapped = build_gene_mapping(genes)
    if not mapped:
        print("no genes with all three regions (5utr/CDS/3utr) found - nothing to plot", file=sys.stderr)
        return

    # Build a per-region bed of just the annotation rows we're using, restricted
    # to genes that made it into `mapped`, so bedtools intersect only ever
    # matches a well-defined region.
    tmp_region_bed = tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
    for gene, entry in mapped.items():
        for region in ("5utr", "CDS", "3utr"):
            for (s, e, cs, ce) in entry[region]["segments"]:
                tmp_region_bed.write(f"{entry['chrom']}\t{s}\t{e}\t{gene}\t{region}\t{entry['strand']}\n")
    tmp_region_bed.close()

    rows = run_bedtools_intersect(args.sites, tmp_region_bed.name)
    os.unlink(tmp_region_bed.name)

    # per-gene, per-region bin counts
    bins_per_region = {"5utr": args.bins_5utr, "CDS": args.bins_cds, "3utr": args.bins_3utr}
    gene_region_bins = defaultdict(lambda: {r: np.zeros(bins_per_region[r]) for r in ("5utr", "CDS", "3utr")})
    gene_total_weight = defaultdict(float)
    gene_site_count = defaultdict(int)

    n_cols_site = 6  # chrom,start,end,name,score,strand for sites.bed
    for parts in rows:
        site_start = int(parts[1])
        ann_gene = parts[n_cols_site + 3]
        ann_region = parts[n_cols_site + 4]
        if ann_gene not in mapped:
            continue
        entry = mapped[ann_gene]
        rel = genomic_to_relative(site_start, entry[ann_region])
        if rel is None:
            continue
        total = entry[ann_region]["total"]
        frac = rel / total if total > 0 else 0.0
        frac = min(max(frac, 0.0), 0.999999)
        b = int(frac * bins_per_region[ann_region])
        try:
            weight = float(parts[4]) if parts[4] != "." else 1.0
        except (ValueError, IndexError):
            weight = 1.0
        gene_region_bins[ann_gene][ann_region][b] += weight
        gene_total_weight[ann_gene] += weight
        gene_site_count[ann_gene] += 1

    qualifying = [g for g in gene_region_bins if gene_site_count[g] >= args.min_sites_per_gene]
    if not qualifying:
        print("no genes passed the min-sites-per-gene threshold - nothing to plot", file=sys.stderr)
        return

    agg = {r: np.zeros(bins_per_region[r]) for r in ("5utr", "CDS", "3utr")}
    for g in qualifying:
        total_w = gene_total_weight[g]
        for r in ("5utr", "CDS", "3utr"):
            agg[r] += gene_region_bins[g][r] / total_w
    for r in ("5utr", "CDS", "3utr"):
        agg[r] /= len(qualifying)

    print(f"{len(mapped)} genes had all three regions; {len(qualifying)} passed the "
          f"min-{args.min_sites_per_gene}-sites threshold and were used", file=sys.stderr)

    combined = np.concatenate([agg["5utr"], agg["CDS"], agg["3utr"]])
    b5, bc, b3 = bins_per_region["5utr"], bins_per_region["CDS"], bins_per_region["3utr"]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(range(len(combined)), combined, color="#2c7fb8", linewidth=1.5)
    ax.axvline(b5, color="black", linestyle="--", linewidth=1)
    ax.axvline(b5 + bc, color="black", linestyle="--", linewidth=1)
    ax.set_xticks([b5 / 2, b5 + bc / 2, b5 + bc + b3 / 2])
    ax.set_xticklabels(["5'UTR", "CDS", "3'UTR"])
    ax.set_ylabel(f"avg. fraction of a gene's (weighted) signal per bin (n={len(qualifying)} genes)")
    ax.set_title("5'UTR -> CDS -> 3'UTR metagene")
    fig.tight_layout()
    fig.savefig(args.plot, dpi=150)


if __name__ == "__main__":
    main()
METAGENE_UTRCDS_PYEOF

cat > "${name}"_metagene_intron.py << 'METAGENE_INTRON_PYEOF'
#!/usr/bin/env python3
"""
metagene_intron.py - intron metagene for PARALUS

Each intron is treated independently (not concatenated per gene, since
introns are usually looked at relative to their own splice sites, not as
part of one continuous transcript coordinate). Position 0 = the intron's
own 5' splice site (in transcript orientation), position 1 = its 3' splice
site. Each qualifying intron's site distribution is normalized to sum to 1
(so one heavily-bound intron can't dominate), then averaged across all
qualifying introns genome-wide.

Optional confidence weighting: column 5 (the BED "score" field) of the
sites file, if numeric, is used as a per-site weight instead of counting
every site equally - so a cluster backed by 200 reads counts more than one
barely passing threshold with 2. A "." or missing column 5 defaults to
weight 1 (fully backward compatible with a plain 6-column sites bed).

Usage:
    metagene_intron.py annotation.bed sites.bed --plot out.png
"""
import sys
import argparse
import subprocess
import tempfile
import os
from collections import defaultdict

import numpy as np


def load_introns(path):
    """Returns list of (chrom, start, end, gene, strand)."""
    introns = []
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            chrom, start, end, gene, annot, strand, biotype = parts[:7]
            if annot != "intron":
                continue
            introns.append((chrom, int(start), int(end), gene, strand))
    return introns


def run_bedtools_intersect(sites_path, region_bed_path):
    cmd = ["bedtools", "intersect", "-wo", "-a", sites_path, "-b", region_bed_path]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True)
    rows = []
    for line in out.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 10:
            continue
        rows.append(parts)
    return rows


def main():
    ap = argparse.ArgumentParser(description="Intron metagene (pooled, per-intron normalized)")
    ap.add_argument("annotation", help="Trierarch annotation bed (canonical recommended)")
    ap.add_argument("sites", help="bed file of crosslink sites or clusters")
    ap.add_argument("--plot", metavar="PNG", required=True)
    ap.add_argument("--bins", type=int, default=100)
    ap.add_argument("--min-sites-per-intron", type=int, default=3,
                     help="introns with fewer than this many sites are excluded, "
                          "to reduce noise from introns with only 1-2 sites")
    ap.add_argument("--min-intron-len", type=int, default=50,
                     help="skip introns shorter than this - very short introns give "
                          "unstable relative-position estimates")
    args = ap.parse_args()

    introns = load_introns(args.annotation)
    introns = [iv for iv in introns if (iv[2] - iv[1]) >= args.min_intron_len]
    if not introns:
        print("no introns long enough to use - nothing to plot", file=sys.stderr)
        return

    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
    intron_id_lookup = {}
    for idx, (chrom, s, e, gene, strand) in enumerate(introns):
        iid = f"{gene}_intron{idx}"
        intron_id_lookup[iid] = (chrom, s, e, strand)
        tmp.write(f"{chrom}\t{s}\t{e}\t{iid}\t.\t{strand}\n")
    tmp.close()

    rows = run_bedtools_intersect(args.sites, tmp.name)
    os.unlink(tmp.name)

    intron_bins = defaultdict(lambda: np.zeros(args.bins))
    intron_total_weight = defaultdict(float)
    intron_site_count = defaultdict(int)

    n_cols_site = 6
    for parts in rows:
        site_start = int(parts[1])
        iid = parts[n_cols_site + 3]
        chrom, s, e, strand = intron_id_lookup[iid]
        length = e - s
        if strand == "+":
            offset = site_start - s
        else:
            offset = e - site_start
        if offset < 0 or offset >= length:
            continue
        frac = min(max(offset / length, 0.0), 0.999999)
        b = int(frac * args.bins)
        try:
            weight = float(parts[4]) if parts[4] != "." else 1.0
        except (ValueError, IndexError):
            weight = 1.0
        intron_bins[iid][b] += weight
        intron_total_weight[iid] += weight
        intron_site_count[iid] += 1

    qualifying = [i for i in intron_bins if intron_site_count[i] >= args.min_sites_per_intron]
    if not qualifying:
        print("no introns passed the min-sites-per-intron threshold - nothing to plot", file=sys.stderr)
        return

    agg = np.zeros(args.bins)
    for i in qualifying:
        agg += intron_bins[i] / intron_total_weight[i]
    agg /= len(qualifying)

    print(f"{len(introns)} introns long enough to use; {len(qualifying)} passed the "
          f"min-{args.min_sites_per_intron}-sites threshold and were used", file=sys.stderr)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(range(args.bins), agg, color="#e6550d", linewidth=1.5)
    ax.axvline(0, color="black", linestyle="--", linewidth=1)
    ax.axvline(args.bins - 1, color="black", linestyle="--", linewidth=1)
    ax.set_xticks([0, args.bins / 2, args.bins - 1])
    ax.set_xticklabels(["5' splice site", "intron (scaled)", "3' splice site"])
    ax.set_ylabel(f"avg. fraction of an intron's (weighted) signal per bin (n={len(qualifying)} introns)")
    ax.set_title("Intron metagene")
    fig.tight_layout()
    fig.savefig(args.plot, dpi=150)


if __name__ == "__main__":
    main()
METAGENE_INTRON_PYEOF

python3 "${name}"_metagene_utr_cds.py "${annotation}" "${name}"_metagene_ttoc.bed \
	--plot "${name}"_ttoc_utr_cds.pdf --min-sites-per-gene 3 2>> "${name}"_stats.txt
python3 "${name}"_metagene_intron.py "${annotation}" "${name}"_metagene_ttoc.bed \
	--plot "${name}"_ttoc_intron.pdf --min-sites-per-intron 3 2>> "${name}"_stats.txt

rm -f "${name}"_metagene_utr_cds.py "${name}"_metagene_intron.py
rm -f "${name}"_metagene_ttoc.bed

fi

echo "timing: total pipeline wall-clock time $(( $(date +%s) - pipeline_start_ts ))s" >> "${name}"_stats.txt

#############
## THE END ##
#############
