#!/bin/bash
set -eo pipefail

eval "$(conda shell.bash hook)"
CWD='/home/asabe/scratch/projects/cassette-exons'

num_threads=64
ulimit -n 10000

cd ${CWD}

# mamba create -p ~/scratch/envs/seqkit_env -y
# mamba install -c bioconda seqkit -y

conda activate /scratch/asabe/envs/seqkit_env
echo $CONDA_PREFIX

for cell_line in MDA231 # HCC1806
do
    for type in Par LM2
    do
        for replicate in 1 2 # 3
        do

            sample_prefix="${cell_line}_${type}_rep${replicate}"
            sample_dir="data/samples/${cell_line}/short-reads"

            if [ -f ${sample_dir}/fastq/${sample_prefix}.fixed.fastq.gz ]
            then          
                echo "SeqKit seq: ${sample_prefix}"
                seqkit seq \
                    --threads ${num_threads} \
                    --min-len 1 \
                    --remove-gaps \
                    ${sample_dir}/fastq/${sample_prefix}.fastq.gz \
                    | gzip \
                        > ${sample_dir}/fastq/${sample_prefix}.fixed.fastq.gz
    
                mv ${sample_dir}/fastq/${sample_prefix}.fastq.gz ${sample_dir}/fastq/${sample_prefix}.old.fastq.gz

                seqkit sana \
                    ${sample_dir}/fastq/${sample_prefix}.fixed.fastq.gz \
                    -o ${sample_dir}/fastq/${sample_prefix}.fixed.rescued.fastq.gz

                zcat ${sample_dir}/fastq/${sample_prefix}.fixed.rescued.fastq.gz \
                | awk '{
                    if (NR % 4 == 1) {
                    id = $0
                    } else if (NR % 4 == 2) {
                    seq = $0
                    } else if (NR % 4 == 3) {
                    plus = $0
                    } else if (NR % 4 == 0) {
                    qual = $0
                    if (length(seq) == length(qual)) {
                        print id
                        print seq
                        print plus
                        print qual
                    }
                    }
                }' \
                | gzip \
                > ${sample_dir}/fastq/${sample_prefix}.fastq.gz


            fi
        done
    done
done

