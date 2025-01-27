#!/usr/bin/env Rscript

# Rscript lorna-manuscript/pre-lornash/long-reads/bambu/sort_gtf.R \
#   --gtf_file /home/asabe/scratch/projects/amir/data/transcriptome/human.c2l2.revio.gencode_v47.GRCh38_p14/extended_annotations.gtf \
#   --output_gtf_file /home/asabe/scratch/projects/amir/data/transcriptome/human.c2l2.revio.gencode_v47.GRCh38_p14/human.c2l2.revio.gencode_v47.GRCh38_p14.extended_annotations.sorted.gtf \
#   --assembly_report_file /scratch/asabe/projects/pacbio/data/references/genome/GRCh38_latest_assembly_report.txt \

# Rscript lorna-manuscript/pre-lornash/long-reads/bambu/sort_gtf.R \
#   --gtf_file /home/asabe/scratch/projects/amir/data/transcriptome/mouse.databank.gencode_vM36.GRCm39/extended_annotations.gtf \
#   --output_gtf_file /home/asabe/scratch/projects/amir/data/transcriptome/mouse.databank.gencode_vM36.GRCm39/mouse.databank.gencode_vM36.GRCm39.extended_annotations.sorted.gtf \
#   --assembly_report_file /scratch/asabe/projects/pacbio/data/references/genome/GCF_000001635.27_GRCm39_assembly_report.txt \

library(optparse)
library(data.table)
library(stringr)
library(magrittr)

option_list <- list(
  make_option("--gtf_file", type = "character", help = "Path to input GTF file", metavar = "file"),
  make_option("--output_gtf_file", type = "character", help = "Path to output sorted GTF file", metavar = "file"),
  make_option("--assembly_report_file", type = "character", help = "Path to assembly report file", metavar = "file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# setwd('/scratch/asabe/projects/lornash')

gtf_file <- opt$gtf_file
output_gtf_file <- opt$output_gtf_file
assembly_report_file <- opt$assembly_report_file

gtf_cols <- c('seqnames', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute')
gtf <- fread(gtf_file,
             col.names = gtf_cols)

assembly_info <- fread(assembly_report_file,
                      skip = '# Sequence-Name',
                      header = T,
                      select = c('GenBank-Accn', 'RefSeq-Accn', 'UCSC-style-name', 'Sequence-Length'),
                      col.names = c('genbank_seqnames', 'refseq_seqnames', 'chr', 'chr_len'))

setnames(assembly_info, 'genbank_seqnames', 'seqnames')
gtf = merge(gtf, assembly_info[, .(seqnames, chr)], by = 'seqnames', all.x = TRUE)
gtf[is.na(chr), chr := seqnames]

### Sorting GTF: first primary chromosomes, then random, then unlocalized
chr_order <- gtf[, .(chr)]
chr_order[, num := chr %>% str_extract('^chr([0-9]+|X|Y|M|Un)') %>% str_remove('chr')]

chr_order[, type := fcase(
  str_detect(chr, 'random'), 'random',
  str_detect(chr, 'alt'), 'alt',
  str_detect(chr, 'fix'), 'fix',
  str_detect(chr, 'chrUn'), 'unlocalized',
  default = 'primary'
)]

chr_order[, num := factor(num, levels = c(seq(22), 'X', 'Y', 'M', 'Un'), ordered = TRUE)]
chr_order[, type := factor(type, levels = c('primary', 'fix', 'random', 'unlocalized', 'alt'), ordered = TRUE)]
chr_order <- unique(chr_order)
setorder(chr_order, type, num)

gtf_chrs <- chr_order[, chr %>% unique()]
gtf[, chr := factor(chr, levels = gtf_chrs, ordered = TRUE)]
gtf <- gtf[order(chr)]

### Adding `attributes`
gtf[, gene_id := 
      attribute %>%
      str_extract('gene_id "(BambuGene[0-9]+|ENSG[0-9]+.[0-9]+)";') %>%
      str_remove('gene_id') %>%
      str_remove_all('"') %>%
      str_remove(';') %>%
      str_trim()]

gtf[, transcript_id := 
      attribute %>%
      str_extract('transcript_id "(BambuTx[0-9]+|ENST[0-9]+.[0-9]+)";') %>%
      str_remove('transcript_id') %>%
      str_remove_all('"') %>%
      str_remove(';') %>%
      str_trim()]
      

gene_order <- gtf[order(chr, start), gene_id %>% unique()]
tr_order <- gtf[order(chr, start), transcript_id %>% unique()]
feature_order <- c('transcript', 'exon')

gtf[, gene_id := factor(gene_id, levels = gene_order, ordered = T)]
gtf[, transcript_id := factor(transcript_id, levels = tr_order, ordered = T)]
gtf[, feature := factor(feature, levels = feature_order, ordered = T)]

setorder(gtf, gene_id, transcript_id, feature)

gtf <- gtf[, gtf_cols, with = F]

fwrite(gtf, sep = '\t', quote = F, col.names = F, row.names = F, file = output_gtf_file)
