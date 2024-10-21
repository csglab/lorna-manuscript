#!/bin/bash
set -eo pipefail

num_threads=128
ulimit -n 100000

project_dir='/scratch/asabe/projects/lornash'
read_len=100
overhang=$((read_len - 1))

cd ${project_dir}

if [ ! -d "data/references/star" ]; then
    echo "Generating genome index..."
    STAR \
        --runThreadN ${num_threads} \
        --runMode genomeGenerate \
        --genomeDir data/references/star \
        --genomeFastaFiles data/references/GRCh38.primary_assembly.genome.fa \
        --sjdbGTFfile data/references/c2l2.annotation.gtf \
        --sjdbOverhang ${overhang}
fi

cell_lines=(A2058 A375tr A549 BxPC3 C42B Colo320 H1299 H358 H23 HCC44 HEK293 HEPG2 HS707A HCT116 K562 LNCaP LS174Ttr MCF7 MDA-231u MDA-453 MP2 PANC1tr PC3 SW480tr SW620 ZR-75-1)
replicates=(1 2)

for cell_line in "${cell_lines[@]}"
do
    for replicate in "${replicates[@]}"
    do
        sample_prefix="${cell_line}_${replicate}"
        sample_dir="data/short-reads"

        mkdir -p ${sample_dir}/bam/${sample_prefix}

        if [ -f ${sample_dir}/fastq/${sample_prefix}_1.fastq.gz ]
        then
            echo "STAR: ${sample_prefix}"
            STAR \
                --runThreadN ${num_threads} \
                --genomeDir data/references/star \
                --readFilesIn ${sample_dir}/fastq/${sample_prefix}_1.fastq.gz ${sample_dir}/fastq/${sample_prefix}_2.fastq.gz \
                --outFileNamePrefix ${sample_dir}/bam/${sample_prefix}/${sample_prefix}. \
                --outSAMunmapped Within \
                --outSAMtype BAM SortedByCoordinate \
                --readFilesCommand zcat \
                --twopassMode Basic \
                --outSAMattributes Standard \
                --sjdbGTFfile data/references/c2l2.annotation.gtf

            samtools index -@ 32 ${sample_dir}/bam/${sample_prefix}/${sample_prefix}.Aligned.sortedByCoord.out.bam
        fi
    done
done