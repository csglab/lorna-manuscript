#!/bin/bash
set -eo pipefail

# Define paths and variables
INPUT_CSV="res.csv"  # Replace with your actual CSV file path
VIS_DIR="/scratch/asabe/projects/c2l2/vis"
LOG_FILE="${VIS_DIR}/output_log.csv"

# Prepare the log file if it doesn't exist
if [ ! -f "$LOG_FILE" ]; then
    echo "event_id,cell_line,cell_line_2,event_region,output_png" > "$LOG_FILE"
fi

# Function to extract event coordinates and cell lines, then run the build_input_bams.sh and run_ggsashimi.sh scripts
process_event() {
    local event_id=$1
    local pos_mutation=$2
    local tissue=$3
    local cell_line=$4
    local cell_line_2=$5
    
    # Extract the chromosome and coordinates from event_id
    local chr=$(echo "$event_id" | awk -F'[:;]' '{print $3}')
    local start=$(echo "$event_id" | awk -F'[:;]' '{print $4}' | awk -F'-' '{print $1}')
    local end=$(echo "$event_id" | awk -F'[:;]' '{print $5}' | awk -F'-' '{print $2}')

    local numeric_part=$(echo "$event_id" | grep -oP '\d+-\d+:\d+-\d+' | sed 's/:/_/')

    # Adjust start and end by 300 bases
    local adj_start=$((start - 300))
    local adj_end=$((end + 300))

    # Form the event region string
    local event_region="${chr}:${adj_start}-${adj_end}"

    # Run the build_input_bams.sh script for the two cell lines
    echo $cell_line $cell_line_2
    ./build_input_bams.sh "$cell_line" "$cell_line_2"

    # Create the output filename for the sashimi plot
    # local output_filename="${VIS_DIR}/${chr}_${adj_start}-${adj_end}.png"
    local output_filename="${VIS_DIR}/${tissue}_${cell_line}_${cell_line_2}_${numeric_part}.png"

    # Run the run_ggsashimi.sh script with the event and cell lines
    echo $event_region $output_filename
    ./run_ggsashimi.sh -c "$event_region" -o "$output_filename"

    # Log the event and output PNG to the CSV log
    echo "$event_id,$cell_line,$cell_line_2,$event_region,$output_filename" >> "$LOG_FILE"
}

# Read the CSV file and process each event
while IFS=',' read -r event_id tissue cell_line cell_line_2; do

    # Check if the required fields (cell_line and cell_line_2) are present
    if [[ -n "$cell_line" && -n "$cell_line_2" ]]; then
        # Call the function to process the event
        process_event "$event_id" "$pos_mutation" "$tissue" "$cell_line" "$cell_line_2"
    else
        echo "Skipping event $event_id due to missing cell line information."
    fi
done < "$INPUT_CSV"

echo "Processing completed. Check $LOG_FILE for details."
