#!/bin/bash

# Script to build liftover chain file from two reference genome fasta files
# Usage: ./build_liftover_chain.sh <source_genome.fa> <target_genome.fa> <output_prefix>

set -e  # Exit on any error

# Check if correct number of arguments provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <source_genome.fa> <target_genome.fa> <output_prefix>"
    echo "Example: $0 hg38.fa hg19.fa hg38ToHg19"
    exit 1
fi

# Input parameters
SOURCE_GENOME="$1"
TARGET_GENOME="$2"
OUTPUT_PREFIX="$3"

# Check if input files exist
if [ ! -f "$SOURCE_GENOME" ]; then
    echo "Error: Source genome file '$SOURCE_GENOME' not found!"
    exit 1
fi

if [ ! -f "$TARGET_GENOME" ]; then
    echo "Error: Target genome file '$TARGET_GENOME' not found!"
    exit 1
fi

# Create working directory
WORK_DIR="${OUTPUT_PREFIX}_liftover_work"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

echo "Starting liftover chain construction..."
echo "Source genome: $SOURCE_GENOME"
echo "Target genome: $TARGET_GENOME"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Working directory: $WORK_DIR"

# Step 1: Create 2bit files from fasta files
echo "Step 1: Converting fasta files to 2bit format..."
faToTwoBit "../$SOURCE_GENOME" source.2bit
faToTwoBit "../$TARGET_GENOME" target.2bit

# Step 2: Create chromosome size files
echo "Step 2: Creating chromosome size files..."
twoBitInfo source.2bit source.chrom.sizes
twoBitInfo target.2bit target.chrom.sizes

# Step 3: Split source genome into chunks for parallel processing
echo "Step 3: Splitting source genome into chunks..."
mkdir -p source_chunks
faSplit size "../$SOURCE_GENOME" 10000000 source_chunks/chunk_

# Step 4: Run lastz alignment for each chunk
echo "Step 4: Running lastz alignments..."
mkdir -p psl_files

# Get list of chunk files
for chunk_file in source_chunks/chunk_*.fa; do
    if [ -f "$chunk_file" ]; then
        chunk_name=$(basename "$chunk_file" .fa)
        echo "Processing $chunk_name..."
        
        # Run lastz alignment
        lastz "$chunk_file" "../$TARGET_GENOME" \
            --format=psl \
            --chain \
            --gapped \
            --step=19 \
            --maxgap=50 \
            --ydrop=3400 \
            --gappedthresh=6000 \
            --scores=HoxD55.q \
            > "psl_files/${chunk_name}.psl"
    fi
done

# Step 5: Combine all PSL files
echo "Step 5: Combining PSL files..."
cat psl_files/*.psl > combined.psl

# Step 6: Convert PSL to chain format
echo "Step 6: Converting PSL to chain format..."
axtChain -psl -verbose=0 combined.psl source.2bit target.2bit stdout | \
    chainAntiRepeat source.2bit target.2bit stdin "${OUTPUT_PREFIX}.chain"

# Step 7: Sort the chain file
echo "Step 7: Sorting chain file..."
chainSort "${OUTPUT_PREFIX}.chain" "${OUTPUT_PREFIX}.sorted.chain"

# Step 8: Create net file (optional but recommended)
echo "Step 8: Creating net file..."
chainNet "${OUTPUT_PREFIX}.sorted.chain" source.chrom.sizes target.chrom.sizes \
    "${OUTPUT_PREFIX}.net" /dev/null

# Step 9: Extract best chains
echo "Step 9: Extracting best chains..."
netChainSubset "${OUTPUT_PREFIX}.net" "${OUTPUT_PREFIX}.sorted.chain" \
    "${OUTPUT_PREFIX}.over.chain"

# Step 10: Move final files to parent directory
echo "Step 10: Moving final files..."
mv "${OUTPUT_PREFIX}.over.chain" "../${OUTPUT_PREFIX}.over.chain"
mv "${OUTPUT_PREFIX}.net" "../${OUTPUT_PREFIX}.net"
mv source.chrom.sizes "../${OUTPUT_PREFIX}.source.chrom.sizes"
mv target.chrom.sizes "../${OUTPUT_PREFIX}.target.chrom.sizes"

# Clean up intermediate files (optional)
echo "Cleaning up intermediate files..."
cd ..
rm -rf "$WORK_DIR"

echo "Liftover chain construction completed!"
echo "Output files:"
echo "  - ${OUTPUT_PREFIX}.over.chain (main liftover chain file)"
echo "  - ${OUTPUT_PREFIX}.net (alignment net file)"
echo "  - ${OUTPUT_PREFIX}.source.chrom.sizes (source genome chromosome sizes)"
echo "  - ${OUTPUT_PREFIX}.target.chrom.sizes (target genome chromosome sizes)"
echo ""
echo "Usage example:"
echo "  liftOver input.bed ${OUTPUT_PREFIX}.over.chain output.bed unmapped.bed"