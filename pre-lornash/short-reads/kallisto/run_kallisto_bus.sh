#!/bin/bash
set -eo pipefail

docker pull quay.io/biocontainers/kallisto:0.51.0--h6de1650_0

project_dir='/scratch/asabe/projects/lornash'
cd ${project_dir}

## Short-reads
fastq_dir='data/short-reads/fastq'
kallisto_bus_dir="data/short-reads/bus"
mkdir -p ${kallisto_bus_dir}

kallisto_index='data/references/c2l2.index.kallisto'
num_threads=32

run_settings="--rm --volume ${project_dir}:/csglab/mpaqt/projects --workdir /csglab/mpaqt/projects --env TMPDIR=/csglab/mpaqt/projects"

### Preparing short-reads data
fastq_pairs_tsv="${fastq_dir}/cell_line_fastq_pairs.tsv"
> ${fastq_pairs_tsv}

for r1_file in ${fastq_dir}/*_R1_001.fastq.gz; do
    base_name=$(basename $r1_file | cut -d'_' -f1-2)
    r2_file="${fastq_dir}/${base_name}_$(basename $r1_file | cut -d'_' -f3-4)_R2_001.fastq.gz"
    echo -e "${base_name}\t${r1_file}\t${r2_file}" >> ${fastq_pairs_tsv}
done

while IFS=$'\t' read -r cell_line fastq1 fastq2; do

    cell_output_dir="${kallisto_bus_dir}/${cell_line}"
    
    # Check if the Kallisto BUS results already exist
    if [ -f "${cell_output_dir}/output.bus" ] && [ -f "${cell_output_dir}/matrix.ec" ]; then
        echo "Skipping ${cell_line}, results already exist."
        continue
    fi

    echo "Kallisto BUS: ${cell_line}"
    mkdir -p ${cell_output_dir}

    docker run \
        ${run_settings} \
        quay.io/biocontainers/kallisto:0.51.0--h6de1650_0 \
        kallisto bus \
            --index ${kallisto_index} \
            --threads ${num_threads} \
            --num \
            --paired \
            --technology bulk \
            --output-dir ${cell_output_dir} \
            ${fastq1} \
            ${fastq2}

done < ${fastq_pairs_tsv}