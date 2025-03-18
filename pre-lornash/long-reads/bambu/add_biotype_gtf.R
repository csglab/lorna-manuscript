#!/usr/bin/env Rscript

# Rscript lorna-manuscript/pre-lornash/long-reads/bambu/add_biotype_gtf.R \
#   --input_gtf_file data/transcriptome/human.c2l2.revio.gencode_v47.GRCh38_p14/human.c2l2.revio.gencode_v47.GRCh38_p14.extended_annotations.sorted.gtf \
#   --ref_gtf_file data/references/gencode.v47.annotation.gtf \
#   --output_gtf_file data/transcriptome/human.c2l2.revio.gencode_v47.GRCh38_p14/human.c2l2.revio.gencode_v47.GRCh38_p14.extended_annotations.sorted.biotyped.gtf

# Rscript lorna-manuscript/pre-lornash/long-reads/bambu/add_biotype_gtf.R \
#   --input_gtf_file data/transcriptome/neurondiff_isoseq.gencode_v47.GRCh38_p14/neurondiff_isoseq.gencode_v47.GRCh38_p14.extended_annotations.sorted.gtf \
#   --ref_gtf_file data/references/gencode.v47.annotation.gtf \
#   --output_gtf_file data/transcriptome/neurondiff_isoseq.gencode_v47.GRCh38_p14/neurondiff_isoseq.gencode_v47.GRCh38_p14.extended_annotations.sorted.biotyped.gtf

library(optparse)
library(data.table)
library(stringr)
library(magrittr)
library(rtracklayer) %>% suppressPackageStartupMessages()
library(GenomicRanges) %>% suppressPackageStartupMessages()

option_list <- list(
  make_option("--input_gtf_file", type = "character", help = "Path to input GTF file", metavar = "file"),
  make_option("--ref_gtf_file", type = "character", help = "Path to reference GTF file", metavar = "file"),
  make_option("--output_gtf_file", type = "character", help = "Path to output GTF file", metavar = "file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# setwd('/scratch/asabe/projects/lornash')

input_gtf_file <- opt$input_gtf_file
ref_gtf_file <- opt$ref_gtf_file
output_gtf_file <- opt$output_gtf_file

# input_gtf_file <- 'data/transcriptome/human.c2l2.revio.gencode_v47.GRCh38_p14/human.c2l2.revio.gencode_v47.GRCh38_p14.extended_annotations.sorted.gtf'
# ref_gtf_file <- 'data/references/gencode.v47.annotation.gtf'
# output_gtf_file <- 'data/transcriptome/human.c2l2.revio.gencode_v47.GRCh38_p14/human.c2l2.revio.gencode_v47.GRCh38_p14.extended_annotations.sorted.biotyped.gtf'

gtf<- import(input_gtf_file)
ref <- import(ref_gtf_file)

gtf <- as.data.table(gtf)
ref <- as.data.table(ref)

gene2type <- ref[, .(gene_id, gene_type)] %>% unique()
gtf <- merge(gtf, gene2type, by = 'gene_id', all.x = T, sort = F)
gtf[is.na(gene_type), gene_type := 'novel']

tr2type <- ref[, .(transcript_id, transcript_type)] %>% unique()
gtf <- merge(gtf, tr2type, by = 'transcript_id', all.x = T, sort = F)
gtf[is.na(transcript_type), transcript_type := 'novel']

gtf_dt <- gtf

gtf_gr <- GRanges(
  seqnames = gtf_dt$seqnames,
  ranges = IRanges(start = gtf_dt$start, end = gtf_dt$end),
  strand = gtf_dt$strand,
  transcript_id = gtf_dt$transcript_id,
  gene_id = gtf_dt$gene_id,
  source = gtf_dt$source,
  type = gtf_dt$type,
  score = gtf_dt$score,
  phase = gtf_dt$phase,
  exon_number = gtf_dt$exon_number,
  gene_type = gtf_dt$gene_type,
  transcript_type = gtf_dt$transcript_type
)

export(gtf_gr, output_gtf_file, format = "gtf")
