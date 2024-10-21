#!/bin/bash

PROJECT_DIR='/scratch/asabe/projects/miniencode-dev'
REF_DIR="${PROJECT_DIR}/data/references"

#### Reference genome ####
mkdir -p ${REF_DIR}/genome

wget \
    'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz' \
    -P ${REF_DIR}/genome

gunzip --keep \
    ${REF_DIR}/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz \
    > ${REF_DIR}/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta

#### Reference annotation ####
GENCODE_VERSION='43'

mkdir -p ${REF_DIR}/annotation

wget \
    "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}/gencode.v${GENCODE_VERSION}.primary_assembly.annotation.gtf.gz" \
    -P ${REF_DIR}/annotation

gunzip --keep \
    ${REF_DIR}/annotation/gencode.v${GENCODE_VERSION}.primary_assembly.annotation.gtf.gz

#### Reference TSS (CAGE-Peaks) ####
REF_TSS_VERSION='4.1'

mkdir -p ${REF_DIR}/cage-peaks

[[ $REF_TSS_VERSION == '3.3' ]] && \
    REF_TSS_LINK='https://reftss.riken.jp/datafiles/3.3/human/refTSS_v3.3_human_coordinate.hg38.bed.gz'

[[ $REF_TSS_VERSION == '4.1' ]] && \
    REF_TSS_LINK='https://reftss.riken.jp/datafiles/4.1/human/refTSS.hg38.4.1.bed.gz'

wget \
    ${REF_TSS_LINK} \
    -P ${REF_DIR}/cage-peaks/

REF_TSS_FILE=$(basename ${REF_TSS_LINK})

gunzip --keep \
    ${REF_DIR}/cage-peaks/${REF_TSS_FILE}

REF_TSS_FILE=${REF_TSS_FILE%.gz}

## Removing ALT scaffolds to be consistent with GENCODE seqname notation
grep -v "_alt" ${REF_DIR}/cage-peaks/${REF_TSS_FILE} > ${REF_DIR}/cage-peaks/${REF_TSS_FILE/.bed/.no_alt.bed}

#### Intropolis ####
mkdir -p ${REF_DIR}/intropolis

## `intropolis` is a list of exon-exon junctions found across 21,504 human RNA-seq samples on the Sequence Read Archive (SRA) from spliced read alignment to hg19 with Rail-RNA.

## Download instructions:
## - Full dataset: https://github.com/nellore/intropolis/blob/master/README.md
## - min_count_10.modified (into STAR junction format): https://github.com/Magdoll/images_public/tree/master/SQANTI2_support_data
## - min_count_10.modified2: https://downloads.pacbcloud.com/public/dataset/MAS-Seq/REF-pigeon_ref_sets/Human_hg38_Gencode_v39/

wget \
    'https://github.com/Magdoll/images_public/raw/master/SQANTI2_support_data/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified.gz' \
    -O ${REF_DIR}/intropolis/intropolis.v1.hg19_with_liftover_to_hg38.min_count_10.modified.tsv.gz

gunzip --keep \
    ${REF_DIR}/intropolis/intropolis.v1.hg19_with_liftover_to_hg38.min_count_10.modified.tsv.gz

wget \
    'https://downloads.pacbcloud.com/public/dataset/MAS-Seq/REF-pigeon_ref_sets/Human_hg38_Gencode_v39/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified2.sorted.tsv' \
    -O ${REF_DIR}/intropolis/intropolis.v1.hg19_with_liftover_to_hg38.min_count_10.modified2.tsv

#### polyA Motifs ####

mkdir -p ${REF_DIR}/polyA-motifs

## The most common polyA motifs for human and mouse
wget \
    'https://raw.githubusercontent.com/ConesaLab/SQANTI3/master/data/polyA_motifs/mouse_and_human.polyA_motif.txt' \
    -O ${REF_DIR}/polyA-motifs/polyA_motifs_list.txt

## Another link: https://downloads.pacbcloud.com/public/dataset/MAS-Seq/REF-pigeon_ref_sets/Human_hg38_Gencode_v39/polyA.list.txt
