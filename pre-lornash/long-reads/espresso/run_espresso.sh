#!/bin/bash
set -eo pipefail

TMPDIR='/scratch/asabe/tmp'
export TMPDIR=${TMPDIR}

ESPRESSO_S=$(which ESPRESSO_S.pl)
ESPRESSO_C=$(which ESPRESSO_C.pl)
ESPRESSO_Q=$(which ESPRESSO_Q.pl)

cell_line=$1
samples_tsv="samples/${cell_line}.espresso.bams.tsv"
ref_genome_fasta='references/GRCh38_no_alt_analysis_set_GCA_000001405.15.id_trimmed.fasta'
ref_annot_gtf='references/gencode.v46.annotation.gtf'
output_dir="espresso/${cell_line}"
num_threads=5
buffer_size=8G

mkdir -p ${output_dir}

perl ${ESPRESSO_S} \
    --list_samples ${samples_tsv} \
    --fa ${ref_genome_fasta} \
    --anno ${ref_annot_gtf} \
    --out ${output_dir} \
    --read_num_cutoff 2 \
    --read_ratio_cutoff 0 \
    --cont_del_max 50 \
    --mapq_cutoff 1 \
    --chrM chrM \
    --num_thread ${num_threads} \
    --sort_buffer_size ${buffer_size}

updated_samples_tsv="${output_dir}/$(basename ${samples_tsv}).updated"

cat ${updated_samples_tsv} \
    | cut -d$'\t' -f3 \
    | while read -r target_id
        do
            perl ${ESPRESSO_C} \
                --in ${output_dir} \
                --fa ${ref_genome_fasta} \
                --target_ID ${target_id} \
                --num_thread ${num_threads} \
                --sort_buffer_size ${buffer_size} \
                &
        done

wait

perl ${ESPRESSO_Q} \
    --list_samples ${updated_samples_tsv} \
    --anno ${ref_annot_gtf} \
    --out_dir ${output_dir}\
    --tsv_compt ${output_dir}/${cell_line}.compatible_isoforms.tsv \
    --num_thread ${num_threads} \
    --read_num_cutoff 2 \
    --read_ratio_cutoff 0
    # --SJ_dist 35 \
    # --internal_boundary_limit 6 \
    # --allow_longer_terminal_exons
