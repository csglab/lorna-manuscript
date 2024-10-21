#!/bin/bash

ref_genome='data/references/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta'
ref_annotation='data/references/annotation/gencode.v43.annotation.gtf'
clustered_hq_fastq='data/long-reads/batches-processed/batches-clustered/flnc_reads.clustered.hq.fasta.gz'
output_dir='data/long-reads/batches-processed/isoquant-results'

isoquant.py \
    --output ${output_dir}/clustered_hq_fastq \
    --data_type pacbio_ccs \
    --fl_data \
    --reference ${ref_genome} \
    --genedb ${ref_annotation} \
    --complete_genedb \
    --fastq ${clustered_hq_fastq} \
    --sqanti_output \
    --threads 248
