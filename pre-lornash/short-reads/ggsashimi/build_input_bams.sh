#!/bin/bash
set -eo pipefail

# Check if exactly two arguments (cell lines) are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <cell_line_1> <cell_line_2>"
    exit 1
fi

# Assign input arguments to variables
CELL_LINE1=$1
CELL_LINE2=$2

VIS_DIR="/scratch/asabe/projects/c2l2/vis"

# Define the output file
OUTPUT_FILE="${VIS_DIR}/input_bams.tsv"
# Remove the output file if it exists
if [ -f "$OUTPUT_FILE" ]; then
    rm "$OUTPUT_FILE"
fi

# Function to write bam file entries to the output file
write_entries () {
    local CELL_LINE=$1
    echo -e "${CELL_LINE}_1\t/scratch/asabe/projects/c2l2/data/short-reads/bam/${CELL_LINE}_1/${CELL_LINE}_1.Aligned.sortedByCoord.out.bam\t${CELL_LINE}" >> $OUTPUT_FILE
    echo -e "${CELL_LINE}_2\t/scratch/asabe/projects/c2l2/data/short-reads/bam/${CELL_LINE}_2/${CELL_LINE}_2.Aligned.sortedByCoord.out.bam\t${CELL_LINE}" >> $OUTPUT_FILE
}

# Write entries for the two cell lines
write_entries $CELL_LINE1
write_entries $CELL_LINE2

# Inform the user of completion
# echo "BAM entries for $CELL_LINE1 and $CELL_LINE2 added to $OUTPUT_FILE."
