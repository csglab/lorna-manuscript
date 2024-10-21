#!/bin/bash
set -eo pipefail

# docker pull guigolab/ggsashimi:v1.1.5
# asabe@emad-2:~/scratch/projects/cassette-exons/data/samples$ find . -type f -name "*.Aligned.sortedByCoord.out.bam*" -exec ln -f {} /home/asabe/scratch/projects/cassette-exons/ggsashimi/bams/ \;
# asabe@emad-2:~/scratch/projects/pacbio/data/cassette-exons$ find . -type f -name "*.flnc.sorted.bam*" -exec ln -f {} /home/asabe/scratch/projects/cassette-exons/ggsashimi/bams/ \;

CWD='/scratch/asabe/projects/cassette-exons/ggsashimi'

cell_lines=('HCC1806' 'MDA231')
ref_annot_gtf='references/gencode.v46.annotation.gtf'
color_palette='ggsashimi/examples/palette.txt'
flank=128

# cell_line='HCC1806'

for cell_line in ${cell_lines[@]}
do
    echo "Processing ${cell_line}"
    bams_tsv="bams/${cell_line}.bams.tsv"
    suppa_dpsi="dpsi/${cell_line}.SE.Par_vs_LM2.dpsi"
    
    output_dir=${suppa_dpsi%.dpsi}
    plots_dir="plots/${output_dir}"
    beds_dir="beds/${output_dir}"
    mkdir -p ${plots_dir} ${beds_dir}
    # docker run -it --entrypoint /bin/bash --volume ${CWD}:${CWD} --workdir ${CWD} guigolab/ggsashimi:v1.1.5

    tail -n +2 "$suppa_dpsi" | while IFS=$'\t' read -r event_id dpsi pval; do
        if (( $(echo "$pval <= 0.05" | bc -l) )); then
            IFS=';' read -r gene_info event_info <<< "$event_id"
            gene_id="${gene_info%%.*}"
            event_type="${event_info%%:*}"
            chr="${event_info#*:}"
            chr="${chr%%:*}"
            rest="${event_info#*:$chr:}"
            e1="${rest%%-*}"
            rest="${rest#*-}"
            s2="${rest%%:*}"
            rest="${rest#*:}"
            e2="${rest%%-*}"
            s3="${rest#*-}"
            s3=${s3%:*}
            strand=${rest#*:}
            
            # if [[ $strand == '+' ]]; then
            #     strand='SENSE'
            # else
            #     strand='ANTISENSE'
            # fi
            
            # echo $chr $e1 $s2 $e2 $s3 $strand
            sashimi_s=$((e1-flank))
            sashimi_e=$((s3+flank))

            event="${chr}:${sashimi_s}-${sashimi_e}"
            event_id="${chr}_${s2}_${e2}"
            echo "${event_id} :: ${event}"

            # ./ggsashimi.py \
            docker run --rm --volume ${CWD}:${CWD} --workdir ${CWD} \
                guigolab/ggsashimi:v1.1.5 \
                    --bam ${bams_tsv} \
                    --gtf ${ref_annot_gtf} \
                    --coordinates ${event} \
                    --out-prefix ${plots_dir}/${event_id} \
                    --min-coverage 10 \
                    --alpha 0.5 \
                    --labels 3 \
                    --color-factor 4 \
                    --base-size 16 \
                    --ann-height 5 \
                    --height 4 \
                    --width 20 \
                    --out-format pdf \
                    --palette ${color_palette} \
                    --junctions-bed ${beds_dir}/${event_id}

                    # --overlay 5 \
                    # --shrink
                    # --fix-y-scale
        fi
    done
done