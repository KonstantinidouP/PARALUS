#!/bin/bash
# =============================================================================
# download_references.sh
# Helper script to download reference files required by PARALUS
#
# This script downloads:
#   1. GENCODE GTF annotation file (human hg38, latest release)
#   2. tRNA BED file from GtRNAdb (human hg38)
#
# IMPORTANT: Links in this script may become outdated as databases release
# new versions. If a download fails, please visit the source databases
# directly to find the current download links:
#
#   GENCODE GTF : https://www.gencodegenes.org/human/
#   tRNA BED    : https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/
#
# Usage:
#   bash download_references.sh
#   bash download_references.sh -o /path/to/output/folder
#
# =============================================================================

# --- Default output directory (current folder) ---
outdir="."

# --- Parse options ---
while getopts o: OPT; do
    case "$OPT" in
        o) outdir=$OPTARG ;;
    esac
done

mkdir -p "$outdir"
cd "$outdir" || exit 1

echo ""
echo "=============================================="
echo "  PARALUS Reference File Downloader"
echo "=============================================="
echo ""
echo "Files will be saved to: $(pwd)"
echo ""

# =============================================================================
# 1. GENCODE GTF
# =============================================================================
echo "----------------------------------------------"
echo " Step 1: Downloading GENCODE GTF annotation"
echo "----------------------------------------------"
echo ""
echo " Source  : GENCODE (https://www.gencodegenes.org)"
echo " Release : v49 (hg38 / GRCh38.p14, freeze 02.2025)"
echo ""
echo " NOTE: If this download fails, visit https://www.gencodegenes.org/human/"
echo " and download the latest 'Comprehensive gene annotation' GTF file."
echo ""

GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
GTF_FILE="gencode.v49.annotation.gtf.gz"

if [ -f "gencode.v49.annotation.gtf" ]; then
    echo " gencode.v49.annotation.gtf already exists, skipping download."
elif [ -f "$GTF_FILE" ]; then
    echo " $GTF_FILE already downloaded, decompressing..."
    gunzip "$GTF_FILE"
else
    echo " Downloading $GTF_FILE ..."
    wget -c "$GTF_URL"
    if [ $? -ne 0 ]; then
        echo ""
        echo " ERROR: GTF download failed."
        echo " Please visit https://www.gencodegenes.org/human/ to download manually."
        echo ""
    else
        echo " Decompressing..."
        gunzip "$GTF_FILE"
        echo " Done: gencode.v49.annotation.gtf"
    fi
fi

echo ""

# =============================================================================
# 2. tRNA BED file from GtRNAdb
# =============================================================================
echo "----------------------------------------------"
echo " Step 2: Downloading tRNA BED file"
echo "----------------------------------------------"
echo ""
echo " Source  : Genomic tRNA Database / GtRNAdb (https://gtrnadb.ucsc.edu)"
echo " Genome  : Human hg38"
echo ""
echo " NOTE: If this download fails, visit:"
echo " https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/"
echo " and click 'Download tRNAscan-SE Results' to get the tar.gz archive."
echo " Extract it and look for the file ending in -tRNAs.bed"
echo ""

TRNA_URL="https://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz"
TRNA_ARCHIVE="hg38-tRNAs.tar.gz"
TRNA_BED="hg38-tRNAs.bed"

if [ -f "$TRNA_BED" ]; then
    echo " $TRNA_BED already exists, skipping download."
else
    echo " Downloading $TRNA_ARCHIVE ..."
    wget -c "$TRNA_URL"
    if [ $? -ne 0 ]; then
        echo ""
        echo " ERROR: tRNA download failed."
        echo " Please visit https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/ to download manually."
        echo ""
    else
        echo " Extracting archive..."
        tar -xzf "$TRNA_ARCHIVE"
        if [ -f "$TRNA_BED" ]; then
            echo " Done: $TRNA_BED"
            rm "$TRNA_ARCHIVE"
        else
            echo ""
            echo " WARNING: Could not find $TRNA_BED after extraction."
            echo " Contents of archive:"
            tar -tzf "$TRNA_ARCHIVE" 2>/dev/null || ls .
            echo " Please identify the BED file manually."
            echo ""
        fi
    fi
fi

echo ""

# =============================================================================
# Summary
# =============================================================================
echo "=============================================="
echo " Download summary"
echo "=============================================="
echo ""

for f in "gencode.v49.annotation.gtf" "hg38-tRNAs.bed"; do
    if [ -f "$f" ]; then
        echo " [OK]  $f"
    else
        echo " [MISSING]  $f — please download manually (see notes above)"
    fi
done

echo ""
echo " You can now run masterchef:"
echo ""
echo "   bash masterchef_v2.sh -g gencode.v49.annotation.gtf -t hg38-tRNAs.bed"
echo ""
echo "=============================================="
