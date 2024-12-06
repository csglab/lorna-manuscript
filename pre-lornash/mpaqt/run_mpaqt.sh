#!/bin/bash
set -eo pipefail

docker pull csglab/mpaqt:0.3.1
docker pull quay.io/biocontainers/kallisto:0.51.1--heb0cbe2_0

### Directories
# project_dir='/scratch/asabe/projects/lornash'
project_dir='/scratch/asabe/amir/MPAQT/projects'
project='human.c2l2.revio.gencode_v47.GRCh38_p14'

## Short-reads
fastq_dir="${project}/data/short-reads/fastq"
kallisto_bus_dir="${project}/data/short-reads/bus"
mkdir -p ${kallisto_bus_dir}

## Long-reads
lr_quant_file="${project}/data/long-reads/quant/${project}.counts.csv"
lr_counts_dir="${project}/data/long-reads/counts"
mkdir -p ${lr_counts_dir}

## MPAQT
logs_dir="${project}/data/logs"

### Settings
cd ${project_dir}
index_run_settings="--rm --volume ${project_dir}:/csglab/mpaqt/projects --workdir /csglab/mpaqt/projects --env TMPDIR=/csglab/mpaqt/projects"
run_settings="--rm --volume ${project_dir}:/csglab/mpaqt/projects --workdir /csglab/mpaqt/projects --user $(id -u):$(id -g) --env TMPDIR=/csglab/mpaqt/projects"

num_threads=32

### Create the index
docker run \
    ${index_run_settings} \
    csglab/mpaqt:0.3.1 index \
    --ref_txome ${project}/data/references/${project}.extended_annotations.transcripts.fasta \
    --ref_annot ${project}/data/references/${project}.extended_annotations.sorted.biotyped.gtf \
    --output ${project}/data/references/${project}.index \
    --num_threads ${num_threads}

# ## Index
# mpaqt_index="${project}/data/references/${project}.index.mpaqt"
# kallisto_index="${project}/data/references/${project}.index.kallisto"

# ### Create the project
# docker run \
#     ${run_settings} \
#     csglab/mpaqt:0.3.1 create project \
#         --index ${mpaqt_index} \
#         ${project}

# ### Preparing short-reads data
# fastq_pairs_tsv="${fastq_dir}/cell_line_fastq_pairs.tsv"
# > ${fastq_pairs_tsv}

# for r1_file in ${fastq_dir}/*_R1_001.fastq.gz; do
#     base_name=$(basename $r1_file | cut -d'_' -f1-2)
#     r2_file="${fastq_dir}/${base_name}_$(basename $r1_file | cut -d'_' -f3-4)_R2_001.fastq.gz"
#     echo -e "${base_name}\t${r1_file}\t${r2_file}" >> ${fastq_pairs_tsv}
# done

# while IFS=$'\t' read -r cell_line fastq1 fastq2; do

#     cell_output_dir="${kallisto_bus_dir}/${cell_line}"
    
#     # Check if the Kallisto BUS results already exist
#     if [ -f "${cell_output_dir}/output.bus" ] && [ -f "${cell_output_dir}/matrix.ec" ]; then
#         echo "Skipping ${cell_line}, results already exist."
#         continue
#     fi

#     echo "Kallisto BUS: ${cell_line}"
#     mkdir -p ${cell_output_dir}

#     docker run \
#         ${run_settings} \
#         quay.io/biocontainers/kallisto:0.51.1--heb0cbe2_0 \
#         kallisto bus \
#             --index ${kallisto_index} \
#             --threads ${num_threads} \
#             --num \
#             --paired \
#             --technology bulk \
#             --output-dir ${cell_output_dir} \
#             ${fastq1} \
#             ${fastq2}

# done < ${fastq_pairs_tsv}

# ### Reformatting long-reads data
# header=$(head -n 1 ${lr_quant_file})
# cell_lines=($(echo $header | cut -f3- -d,))

# for ((i=0; i<${#cell_lines[@]}; i++)); do

#     echo "Processing Long-reads Counts: ${cell_lines[$i]}"
#     cell_line=${cell_lines[$i]}
#     mkdir -p ${lr_counts_dir}/${cell_line}
#     output_file="${lr_counts_dir}/${cell_line}/${cell_line}.counts.tsv"
#     cut -d, -f1,$((i+3)) ${lr_quant_file} > ${output_file}
#     sed -i '1s/\t[^\t]*$/\tcount/' ${output_file}

# done

# ### Getting cell-line samples names
# lr_samples=($(ls -d ${lr_counts_dir}/*/ | xargs -n 1 basename))
# sr_samples=($(ls -d ${kallisto_bus_dir}/*/ | xargs -n 1 basename))

# samples=()

# for cell_line in "${sr_samples[@]}"; do
#     if [[ " ${lr_samples[@]} " =~ " ${cell_line} " ]]; then
#         samples+=("$cell_line")
#     fi
# done

# ### Creating sample repositories
# docker run \
#     ${run_settings} \
#     csglab/mpaqt:0.3.1 create sample \
#         --project ${project} \
#         ${samples[@]}

# ### Preparing samples (one-threaded)
# # for sample in ${samples[@]}; do

# #     docker run \
# #     ${run_settings} \
# #     csglab/mpaqt:0.3.1 prepare short-read \
# #         --project ${project} \
# #         --sample ${sample} \
# #         --bus ${kallisto_bus_dir}/${sample}/output.bus \
# #         --matrix_ec ${kallisto_bus_dir}/${sample}/matrix.ec

# #     docker run \
# #     ${run_settings} \
# #     csglab/mpaqt:0.3.1 prepare long-read \
# #         --project ${project} \
# #         --sample ${sample} \
# #         --counts ${lr_counts_dir}/${sample}/${sample}.counts.tsv

# # done

# ### Preparing samples (parallel)
# export run_settings
# export project
# export kallisto_bus_dir
# export lr_counts_dir

# mpaqt_prepare() {
#     sample=$1

#     # Run the short-read prepare command
#     docker run \
#     ${run_settings} \
#     csglab/mpaqt:0.3.1 prepare short-read \
#         --project ${project} \
#         --sample ${sample} \
#         --bus ${kallisto_bus_dir}/${sample}/output.bus \
#         --matrix_ec ${kallisto_bus_dir}/${sample}/matrix.ec

#     # Run the long-read prepare command
#     docker run \
#     ${run_settings} \
#     csglab/mpaqt:0.3.1 prepare long-read \
#         --project ${project} \
#         --sample ${sample} \
#         --counts ${lr_counts_dir}/${sample}/${sample}.counts.tsv
# }

# export -f mpaqt_prepare

# parallel --tag -j ${num_threads} \
#     mpaqt_prepare ::: "${samples[@]}" \
#     >${logs_dir}/mpaqt_prepare.log 2>&1

# ### One-threaded MPAQT Quant
# # for sample in ${samples[@]}; do
# #     docker run \
# #     ${run_settings} \
# #     csglab/mpaqt:0.3.1 quant \
# #         --project ${project} \
# #         --sample ${sample}
# # done

# ### Parallel MPAQT Quant
# export run_settings
# export project

# mpaqt_quant() {
#     sample=$1
#     docker run \
#     ${run_settings} \
#     csglab/mpaqt:0.3.1 quant \
#         --project ${project} \
#         --sample ${sample}
# }

# export -f mpaqt_quant
# mkdir -p ${logs_dir}

# parallel --tag -j ${num_threads} \
#     mpaqt_quant ::: "${samples[@]}" \
#     >${logs_dir}/mpaqt_quant.log 2>&1
