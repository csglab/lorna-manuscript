#!/usr/bin/env Rscript

# Rscript lorna-manuscript/pre-lornash/long-reads/bambu/postprocess_bambu_counts.R \
#  --counts_file data/transcriptome/human.c2l2.revio.gencode_v47.GRCh38_p14/CPM_transcript.txt

library(optparse)
library(data.table)
library(stringr)
library(magrittr)

option_list <- list(
make_option(c("--counts_file"), type = "character", default = "human", help = "Species: either 'human' or 'mouse' [default= %default]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

counts_file <- opt$counts_file
output_prefix <- str_remove(counts_file, "\\.\\w+$")
output_file <- str_glue("{output_prefix}.processed.csv")

# setwd('/scratch/asabe/projects/lornash')

counts <- fread(counts_file)
sample_full_ids <- colnames(counts)[-c(1, 2)]
sample_ids <- sample_full_ids %>% str_extract('_s{1,4}.*.flnc') %>% str_remove('^_') %>% str_remove('.flnc')
cell_lines <- sample_ids %>% str_remove('^s{1,4}..') %>% str_remove('_1$|_2$') %>% unique() %>% sort()
num_cell_lines <- length(cell_lines)
setnames(counts, sample_full_ids, sample_ids)

cell_line_ids <- sample_ids %>% str_split_fixed('[.]', n = 2) %>% extract(, 2) %>% unique()

for (cl in cell_line_ids) {
  counts[, (cl) := rowSums(.SD), .SDcols = patterns(str_glue("s[1-4][.]{cl}"))]
}

counts <- counts[, c('TXNAME', 'GENEID', cell_line_ids), with = FALSE]
setnames(counts, 'TXNAME', 'transcript_id')
setnames(counts, 'GENEID', 'gene_id')
fwrite(counts, output_file, quote = FALSE)
