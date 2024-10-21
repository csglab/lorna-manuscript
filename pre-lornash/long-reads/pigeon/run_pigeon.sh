#!/bin/bash

samtools faidx data/references/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta 
pigeon sort data/long-reads/batches-collapsed/flnc_reads.clustered.hq.aligned.collapsed.gff -o data/long-reads/batches-collapsed/flnc_reads.clustered.hq.aligned.collapsed.sorted.gff
pigeon index data/long-reads/batches-collapsed/flnc_reads.clustered.hq.aligned.collapsed.sorted.gff

pigeon sort data/long-reads/batches-processed/batches-collapsed-extra-exons/flnc_reads.clustered.hq.aligned.collapsed.extra_exons.gff -o data/long-reads/batches-processed/batches-collapsed-extra-exons/flnc_reads.clustered.hq.aligned.collapsed.extra_exons.sorted.gff
pigeon index data/long-reads/batches-processed/batches-collapsed-extra-exons/flnc_reads.clustered.hq.aligned.collapsed.extra_exons.sorted.gff

pigeon sort data/references/annotation/gencode.v43.primary_assembly.annotation.gtf -o data/references/annotation/gencode.v43.primary_assembly.annotation.sorted.gtf
pigeon index data/references/annotation/gencode.v43.primary_assembly.annotation.sorted.gtf

pigeon sort data/references/cage-peaks/refTSS.hg38.4.1.no_alt.bed -o data/references/cage-peaks/refTSS.hg38.4.1.no_alt.sorted.bed
pigeon index data/references/cage-peaks/refTSS.hg38.4.1.no_alt.sorted.bed

pigeon sort data/references/cage-peaks/refTSS_v3.3_human_coordinate.hg38.no_alt.bed -o data/references/cage-peaks/refTSS_v3.3_human_coordinate.hg38.no_alt.sorted.bed
pigeon index data/references/cage-peaks/refTSS_v3.3_human_coordinate.hg38.no_alt.sorted.bed

pigeon sort data/references/intropolis/intropolis.v1.hg19_with_liftover_to_hg38.min_count_10.modified2.tsv -o data/references/intropolis/intropolis.v1.hg19_with_liftover_to_hg38.min_count_10.modified2.sorted.tsv
pigeon index data/references/intropolis/intropolis.v1.hg19_with_liftover_to_hg38.min_count_10.modified2.sorted.tsv

mkdir -p data/long-reads/batches-transcriptome

pigeon classify \
    data/long-reads/batches-collapsed/flnc_reads.clustered.hq.aligned.collapsed.sorted.gff \
    data/references/annotation/gencode.v43.annotation.sorted.gtf \
    data/references/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    --flnc data/long-reads/batches-collapsed/flnc_reads.clustered.hq.aligned.collapsed.abundance.txt \
    --cage-peak data/references/cage-peaks/refTSS_v3.3_human_coordinate.hg38.no_alt.sorted.bed \
    --poly-a data/references/polyA-motifs/polyA_motifs_list.txt \
    --out-dir data/long-reads/batches-transcriptome \
    --num-threads 0 \
    --log-level INFO


pigeon classify \
    data/long-reads/batches-processed/batches-collapsed/flnc_reads.clustered.hq.aligned.collapsed.sorted.gff \
    data/references/annotation/gencode.v43.annotation.sorted.gtf \
    data/references/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    --flnc data/long-reads/batches-processed/batches-collapsed/flnc_reads.clustered.hq.aligned.collapsed.abundance.txt \
    --cage-peak data/references/cage-peaks/refTSS_v3.3_human_coordinate.hg38.no_alt.sorted.bed \
    --poly-a data/references/polyA-motifs/polyA_motifs_list.txt \
    --out-prefix flnc_reads.clustered.hq.aligned.collapsed.sorted \
    --out-dir data/long-reads/batches-transcriptome \
    --num-threads 120 \
    --log-level INFO
