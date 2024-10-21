#!/bin/bash
set -eo pipefail

mkdir -p events psi dpsi

suppa.py generateEvents \
  --input-file references/c2l2.annotation.v2.gtf \
  --output-file events/c2l2 \
  --event-type SE SS MX RI FL \
#   --format ioe
  # # Defaults:
  #   --boundary S \
  #   --threshold 10 \
  #   --exon-length 100 \
  #   --mode INFO \

# SE: Skipping exon (SE) events
# SS: Alternative 5' (A5) and 3' (A3) splice sites (it generates both)
# MX: Mutually Exclusive (MX) exons
# RI: Retained intron (RI)
# FL: Alternative first (AF) and last (AL) exons (it generates both)

events=(SE A5 A3 MX RI AF AL)
# events=(SE)

mkdir -p psi

for event in "${events[@]}"; do
  suppa.py psiPerEvent \
    --ioe-file events/c2l2_${event}_strict.ioe \
    --expression-file quant/mpaqt.tpm.tsv \
    --output-file psi/c2l2.${event} &
done

# mkdir samples
# find /scratch/asabe/csg/mpaqt/MPAQT/projects/cassette-exons/samples \
#   -name "*.MPAQT.LR_SR.tsv" \
#   -type f \
#   -exec cp -v {} ./samples \; 

suppa.py generateEvents \
  --input-file references/gencode.v46.annotation.gtf \
  --output-file events/gencode.v46 \
  --event-type SE SS MX RI FL \
  --format ioe
  ## Defaults:
    # --boundary S \
    # --threshold 10 \
    # --exon-length 100 \
    # --mode INFO \

# # SE: Skipping exon (SE) events
# # SS: Alternative 5' (A5) and 3' (A3) splice sites (it generates both)
# # MX: Mutually Exclusive (MX) exons
# # RI: Retained intron (RI)
# # FL: Alternative first (AF) and last (AL) exons (it generates both)

events=(SE A5 A3 MX RI AF AL)
# events=(SE)
cell_lines=(HCC1806 MDA231)

for cell_line in ${cell_lines[@]}
do
  for type in Par LM2
  do
    bash combine_mpaqt_results.sh ${cell_line}_${type} samples
  done
done

for event in ${events[@]}
do
  for cell_line in ${cell_lines[@]}
  do
    for type in Par LM2
    do
      suppa.py psiPerEvent \
        --ioe-file events/gencode.v46_${event}_strict.ioe \
        --expression-file samples/${cell_line}_${type}.expression.tsv \
        --output-file psi/${cell_line}_${type}.${event}
        ## Defaults:
        # --total-filter 0 \
        # --mode INFO 
    done
  done
done

for event in ${events[@]}
do
  for cell_line in ${cell_lines[@]}
  do
    suppa.py diffSplice \
      --method empirical \
      --psi psi/${cell_line}_Par.${event}.psi psi/${cell_line}_LM2.${event}.psi \
      --tpm samples/${cell_line}_Par.expression.tsv samples/${cell_line}_LM2.expression.tsv \
      --input events/gencode.v46_${event}_strict.ioe \
      --output dpsi/${cell_line}.${event}.Par_vs_LM2 \
      --save_tpm_events

    mv dpsi/${cell_line}.${event}.Par_vs_LM2.dpsi.temp.0 dpsi/${cell_line}.${event}.Par_vs_LM2.dpsi
  done
done



