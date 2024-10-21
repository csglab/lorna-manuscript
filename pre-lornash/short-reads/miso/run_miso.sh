#!/bin/bash
set -eo pipefail

eval "$(conda shell.bash hook)"
CWD='/home/asabe/scratch/projects/cassette-exons'

num_threads=64
ulimit -n 10000

cd ${CWD}

# mamba create -p /scratch/asabe/envs/miso_env python=2.7 -y
# mamba activate /scratch/asabe/envs/miso_env
# mamba install -c bioconda misopy -y

conda activate /scratch/asabe/envs/miso_env
echo $CONDA_PREFIX

###Calculate insert length distribution
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

## Running MISO with single-end reads
for cell_line in MDA231
do
    for type in Par LM2
    do
        for replicate in 1 2 3
        do
            
            sample_prefix="${cell_line}_${type}_rep${replicate}"
            sample_dir="data/samples/${cell_line}/short-reads"

            if [ -f ${sample_dir}/bam/${sample_prefix}.Aligned.sortedByCoord.out.bam ]
            then

                read_len=76

                echo "MISO: ${sample_prefix}..."
                miso \
                    --run \
                        data/references/miso/gff/commonshortest/indexed_SE_hg38_events \
                        ${sample_dir}/bam/${sample_prefix}.Aligned.sortedByCoord.out.bam \
                    --output-dir ${sample_dir}/miso/${sample_prefix}/ \
                    --read-len ${read_len} \
                    --overhang-len 8 \
                    -p ${num_threads}
                
                summarize_miso \
                    --summarize-samples ${sample_dir}/miso/${sample_prefix}/ \
                    ${sample_dir}/miso/${sample_prefix}/ 
            fi
        done
    done
done


## Running MISO with paired-end reads
for cell_line in HCC1806
do
    for type in Par LM2
    do
        for replicate in 1 2 3
        do
            
            sample_prefix="${cell_line}_${type}_rep${replicate}"
            sample_dir="data/samples/${cell_line}/short-reads"

            if [ -f ${sample_dir}/bam/${sample_prefix}.Aligned.sortedByCoord.out.bam ]
            then
            
                insert_len=${sample_dir}/insert-dist/${sample_prefix}.Aligned.sortedByCoord.out.bam.insert_len
                read -r first_line < ${insert_len}
                mean=$(echo $first_line | sed 's/.*mean=\([^,]*\).*/\1/')
                sdev=$(echo $first_line | sed 's/.*sdev=\([^,]*\).*/\1/')
                read_len=151

                echo "MISO: ${sample_prefix}..."
                miso \
                    --run \
                        data/references/miso/gff/commonshortest/indexed_SE_hg38_events \
                        ${sample_dir}/bam/${sample_prefix}.Aligned.sortedByCoord.out.bam \
                    --output-dir ${sample_dir}/miso/${sample_prefix}/ \
                    --read-len ${read_len} \
                    --paired-end ${mean} ${sdev} \
                    --overhang-len 8 \
                    -p ${num_threads}
                
                summarize_miso \
                    --summarize-samples ${sample_dir}/miso/${sample_prefix}/ \
                    ${sample_dir}/miso/${sample_prefix}/ 
            fi
        done
    done
done


for cell_line in HCC1806 MDA231
do
    for replicate in 1 2 3
    do   
        sample_Par="${cell_line}_Par_rep${replicate}"
        sample_LM2="${cell_line}_LM2_rep${replicate}"
        sample_dir="data/samples/${cell_line}/short-reads"

        if [ -f ${sample_dir}/bam/${sample_LM2}.Aligned.sortedByCoord.out.bam ]
        then
            echo "MISO Compare: ${sample_Par} vs. ${sample_LM2}..."
            compare_miso \
                --compare-samples \
                ${sample_dir}/miso/${sample_Par}/ \
                ${sample_dir}/miso/${sample_LM2}/ \
                ${sample_dir}/miso/comparisons/

        fi
    done
done



for cell_line in HCC1806 MDA231
do
    for type in Par LM2
    do
        for replicate in 1 2 3
        do
            
            sample_dir="data/samples/${cell_line}/long-reads"
            if [[ $cell_line == "HCC1806" ]]
            then
                sample_prefix="m84127_240426_234003_s3.HCC1806_${type}_${replicate}"
                read_len=3500
            elif [[ $cell_line == "MDA231" ]]
            then
                sample_prefix="210222_AH_IsoSeq.MDA_MB_231_${type}_duplicate_${replicate}"
                read_len=2000
            fi

            if [ -f ${sample_dir}/bam/${sample_prefix}.flnc.sorted.bam ]
            then
            
                echo "MISO: ${sample_prefix}..."
                miso \
                    --run \
                        data/references/miso/gff/commonshortest/indexed_SE_hg38_events \
                        ${sample_dir}/bam/${sample_prefix}.flnc.sorted.bam \
                    --output-dir ${sample_dir}/miso/${sample_prefix}/ \
                    --read-len ${read_len} \
                    --overhang-len 8 \
                    -p ${num_threads}
            fi
        done
    done
done


for i in 1 2
do
    miso \
        --run \
            data/references/miso/gff/commonshortest/indexed_SE_hg38_events \
            data/long-read/flnc/m64182_230110_195600.HEK293_${i}.flnc_clean.sorted.bam \
        --output-dir data/long-read/miso/HEK293_${i}/ \
        --read-len 3000 \
        -p 16
    
    summarize_miso \
        --summarize-samples data/long-read/miso/HEK293_${i}/ \
         data/long-read/miso/HEK293_${i}/ 

done

# sashimi_plot_event=''
# sashimi_plot \
#     --plot-event ${sashimi_plot_event} \
#     data/references/miso/gff/commonshortest/indexed_SE_hg38_events \
#     plots/sashimi_plot_settings.txt \
#     --output-dir plots/