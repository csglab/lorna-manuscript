#!/bin/bash
set -eo pipefail

# Get event region and output filename from the input arguments
while getopts ":c:o:" opt; do
  case $opt in
    c) event="$OPTARG"
    ;;
    o) output="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

TSV='/scratch/asabe/projects/c2l2/vis/input_bams.tsv'
GTF='/scratch/asabe/projects/c2l2/data/references/c2l2.annotation.gtf'
PAL='/scratch/asabe/projects/c2l2/vis/ggsashimi/examples/palette.txt'
SCR='/scratch/asabe/projects/c2l2/vis/ggsashimi/ggsashimi.py'

echo $event
echo $output

${SCR}  \
        -b ${TSV} \
        -c ${event} \
        -o ${output} \
        -M 1 \
        -A mean \
        -C 3 \
        -O 3 \
        --shrink \
        --alpha 0.5 \
        -j junction \
        -F pdf \
        -P ${PAL} \
        -g ${GTF} \
        --base-size=20 \
        --height=4 \
        --ann-height=9 \
        --width=18

echo "Sashimi plot for ${event} saved to ${output}.png."
