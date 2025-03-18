#!/bin/bash
set -eo pipefail

# Set project directory and change to it
project_dir="/large_storage/goodarzilab/saberi/mpaqt/benchmarks"
cd ${project_dir}

data_dir="${project_dir}/DATA/ILLUMINA_DATA"
bus_output_dir="${project_dir}/processed/bus"
mkdir -p ${bus_output_dir}

# Kallisto index file
kallisto_index="/home/saberi/projects/mpaqt/rebuttal/data/spike-ins/gencode.v47.sequins.v2.4.transcripts.kallisto.index"

# Use 64 threads
num_threads=64

# # Activate conda environment
# source /opt/conda/etc/profile.d/conda.sh
# conda activate kallisto

# Iterate over each sample directory
for sample_dir in ${data_dir}/*; do
    sample_name=$(basename ${sample_dir})
    
    # Define the paired fastq files based on the sample name
    fastq1="${sample_dir}/${sample_name}.sra_1.fastq"
    fastq2="${sample_dir}/${sample_name}.sra_2.fastq"
    
    # Set the output directory for the current sample
    cell_output_dir="${bus_output_dir}/${sample_name}"
    
    # Skip the sample if output files already exist
    if [ -f "${cell_output_dir}/output.bus" ] && [ -f "${cell_output_dir}/matrix.ec" ]; then
        echo "Skipping ${sample_name}, results already exist."
        continue
    fi
    
    echo "Running kallisto bus for ${sample_name}"
    mkdir -p ${cell_output_dir}
    
    kallisto bus \
        --index ${kallisto_index} \
        --threads ${num_threads} \
        --num \
        --paired \
        --technology bulk \
        --output-dir ${cell_output_dir} \
        ${fastq1} \
        ${fastq2}
done

# Deactivate conda environment
# conda deactivate