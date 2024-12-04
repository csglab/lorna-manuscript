#!/usr/bin/env Rscript

# Rscript lorna-manuscript/pre-lornash/long-reads/bambu/extract_transcripts_fasta.R \
#   --gtf_file data/transcriptome/human.c2l2.revio.gencode_v47.GRCh38_p14/human.c2l2.revio.gencode_v47.GRCh38_p14.extended_annotations.sorted.gtf \
#   --fasta_file data/transcriptome/human.c2l2.revio.gencode_v47.GRCh38_p14/human.c2l2.revio.gencode_v47.GRCh38_p14.extended_annotations.transcripts.fasta

library(optparse)
library(data.table)
library(stringr)
library(magrittr)
library(rtracklayer) %>% suppressPackageStartupMessages()
library(Biostrings) %>% suppressPackageStartupMessages()
library(BSgenome.Hsapiens.UCSC.hg38) %>% suppressPackageStartupMessages()

option_list <- list(
  make_option("--gtf_file", type = "character", help = "Path to input GTF file", metavar = "file"),
  make_option("--fasta_file", type = "character", help = "Path to output FASTA file", metavar = "file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# setwd('/scratch/asabe/projects/lornash')

gtf_file <- opt$gtf_file
fasta_file <- opt$fasta_file

ref_genome <- BSgenome.Hsapiens.UCSC.hg38

# setwd('/scratch/asabe/projects/amir')
# gtf_file <- 'data/transcriptome/human.c2l2.revio.gencode_v47.GRCh38_p14/human.c2l2.revio.gencode_v47.GRCh38_p14.extended_annotations.sorted.gtf'
# ref_gtf_file <- 'data/references/gencode.v47.annotation.gtf'
# fasta_file <- 'data/transcriptome/human.c2l2.revio.gencode_v47.GRCh38_p14/human.c2l2.revio.gencode_v47.GRCh38_p14.extended_annotations.transcripts.fasta'

gtf <- import(gtf_file)
gtf <- as.data.table(gtf)

columns_to_keep <- c('seqnames', 'type', 'feature', 'start', 'end', 'width', 'strand', 'gene_id', 'transcript_id', 'exon_number')
columns_to_keep <- intersect(columns_to_keep, colnames(gtf))
gtf <- gtf[, columns_to_keep, with = F]
gtf <- gtf[type %in% c('transcript', 'exon')]
gtf[, exon_number := as.integer(exon_number)]
gtf[, num_exons := max(exon_number, na.rm = T), transcript_id]
gtf[strand == '-', exon_number := (num_exons - exon_number + 1)]

gtf[type != 'transcript', seq_ := getSeq(ref_genome,
                                         names = seqnames,
                                         start = start,
                                         end = end,
                                         as.character = TRUE)]


gtf[, type := factor(type, levels = c('exon', 'transcript'), ordered = TRUE)]

tr_seqs <- gtf[type != 'transcript'][
  order(transcript_id, exon_number, type),
  .(seq_ = str_c(seq_, collapse = '')),
  by = .(transcript_id, strand)]

tr_seqs[strand == '-', seq_ := seq_ %>% DNAStringSet() %>% reverseComplement() %>% as.character()]

tr_seqs[, seq_len := nchar(seq_)]

# ref_fasta = 'data/references/gencode.v47.transcripts.fa'
# ref_seqs = readDNAStringSet(ref_fasta)
# names(ref_seqs) <- names(ref_seqs) %>% str_split('[|]', simplify = T) %>% extract(, 1)

# test = tr_seqs[strand == '+' & transcript_id %like% 'ENST'][2]
# as.character(ref_seqs[names(ref_seqs) == test$transcript_id]) == test$seq_
# 
# test = tr_seqs[strand == '-' & transcript_id %like% 'ENST'][2]
# as.character(ref_seqs[names(ref_seqs) == test$transcript_id]) == test$seq_

fasta = tr_seqs$seq_
names(fasta) = tr_seqs$transcript_id
fasta = DNAStringSet(fasta)

writeXStringSet(fasta, fasta_file, format = 'fasta')

