#!/bin/bash
set -eo pipefail

eval "$(conda shell.bash hook)"
CWD='/home/asabe/scratch/projects/cassette-exons'

num_threads=64
ulimit -n 10000

cd ${CWD}

conda activate /scratch/asabe/envs/miso_env
echo $CONDA_PREFIX

#### Calculate insert length distribution
exon_utils \
    --get-const-exons \
    data/references/annotation/gencode.v43.annotation.gff3 \
    --min-exon-size 1000 \
    --output-dir data/references/annotation/exons/

for cell_line in HCC1806
do
    for type in Par LM2
    do
        for replicate in 1 2 3
        do
            sample_prefix="${cell_line}_${type}_rep${replicate}"
            sample_dir="data/samples/${cell_line}/short-reads"
            sample_bam="${sample_dir}/bam/${sample_prefix}.Aligned.sortedByCoord.out.bam"
        
            if [ -f ${sample_bam} ]
            then
                pe_utils \
                    --compute-insert-len \
                    ${sample_bam} \
                    data/references/annotation/exons/gencode.v43.annotation.min_1000.const_exons.gff \
                    --output-dir ${sample_dir}/insert-dist/
            fi
        done
    done
done