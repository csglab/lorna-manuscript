#!/bin/bash
set -eo pipefail

num_threads=128
ulimit -n 100000

project_dir='/scratch/asabe/projects/lornash'
cd ${project_dir}

# List of cell lines and replicates
cell_lines=(A2058 A375tr A549 BxPC3 C42B Colo320 H1299 H358 H23 HCC44 HEK293 HEPG2 HS707A HCT116 K562 LNCaP LS174Ttr MCF7 MDA-231u MDA-453 MP2 PANC1tr PC3 SW480tr SW620 ZR-75-1)
replicates=(1 2)

transcriptome_fasta="data/references/c2l2.annotation.sorted.gtf"
kallisto_index="data/references/c2l2.index.kallisto"

# Create Kallisto index
if [ ! -f "${kallisto_index}" ]; then
    echo "Creating Kallisto index..."
    kallisto index --index ${kallisto_index} ${transcriptome_fasta}
fi

# Quantify isoform abundance for each cell line and replicate
for cell_line in "${cell_lines[@]}"
do
    for replicate in "${replicates[@]}"
    do
        sample_prefix="${cell_line}_${replicate}"
        sample_dir="data/short-reads"
        output_dir="data/isoform-abundance/${sample_prefix}"

        mkdir -p ${output_dir}

        if [ -f ${sample_dir}/fastq/${sample_prefix}_1.fastq.gz ] && [ -f ${sample_dir}/fastq/${sample_prefix}_2.fastq.gz ]
        then
            echo "Kallisto: ${sample_prefix}"
            kallisto quant \
                --index ${kallisto_index} \
                --output-dir ${output_dir} \
                --plaintext \
                --rf-stranded \
                --threads ${num_threads} \
                ${sample_dir}/fastq/${sample_prefix}_1.fastq.gz \
                ${sample_dir}/fastq/${sample_prefix}_2.fastq.gz

            mv ${output_dir}/abundance.tsv ${output_dir}/${sample_prefix}.rf-stranded.c2l2.abundance.tsv
            mv ${output_dir}/run_info.json ${output_dir}/${sample_prefix}.rf-stranded.c2l2.run_info.json
        fi
    done
done