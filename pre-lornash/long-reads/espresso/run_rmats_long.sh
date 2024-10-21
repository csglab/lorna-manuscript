#!/bin/bash
set -eo pipefail

cell_line=$1

espresso_dir="espresso/${cell_line}"

ref_annot_gtf='references/gencode.v46.annotation.gtf'
updated_annot_gtf="${espresso_dir}/bams_N2_R0_updated.gtf"
abundance_esp="${espresso_dir}/bams_N2_R0_abundance.esp"
group1_txt="samples/${cell_line}_Par.rmats_long.group.txt"
group2_txt="samples/${cell_line}_LM2.rmats_long.group.txt"
group1_name="${cell_line}_Par"
group2_name="${cell_line}_LM2"
output_dir="rmats-long/${cell_line}"
num_threads=32

mkdir -p ${output_dir}

rmats_long.py \
    --abundance ${abundance_esp} \
    --updated-gtf ${updated_annot_gtf} \
    --gencode-gtf ${ref_annot_gtf} \
    --group-1 ${group1_txt} \
    --group-2 ${group2_txt} \
    --group-1-name ${group1_name} \
    --group-2-name ${group2_name} \
    --out-dir ${output_dir} \
    --num-threads ${num_threads} \
    --plot-file-type .pdf \
    --adj-pvalue 0.05 \
    --delta-proportion 0.1
    # --compare-all-within-gene
