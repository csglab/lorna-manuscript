#!/usr/bin/env Rscript

# Rscript lorna-manuscript/pre-lornash/long-reads/bambu/run_bambu.R \
#   --species human \
#   --mode quant \
#   --do_postprocess TRUE \
#   --num_threads 248

# Rscript lorna-manuscript/pre-lornash/long-reads/bambu/run_bambu.R \
#   --species mouse \
#   --mode quant \
#   --do_postprocess FALSE \
#   --num_threads 64

# Rscript lorna-manuscript/pre-lornash/long-reads/bambu/run_bambu.R \
#   --species human \
#   --mode quant \
#   --do_postprocess FALSE \
#   --num_threads 64 \
#   --data_dir /scratch/asabe/projects/pacbio/data/databank/human/bam \
#   --bam_pattern 'neurondiff_isoseq.*.flnc.sorted.bam$' \
#   --prefix neurondiff_isoseq.gencode_v47.GRCh38_p14

library(optparse)
library(data.table)
library(stringr)
library(magrittr)
library(bambu)


option_list <- list(
  make_option(
    c("--species"),
    type = "character",
    default = "human",
    help = "Species: either 'human' or 'mouse' [default= %default]",
    metavar = "character"
  ),
  make_option(
    c("--mode"),
    type = "character",
    help = "Mode: either 'discovery' or 'quant'",
    metavar = "character"
  ),
  make_option(
    c("--do_postprocess"),
    type = "logical",
    default = FALSE,
    help = "Perform post-processing [default= %default]",
    metavar = "logical"
  ),
  make_option(
    c("--num_threads"),
    type = "integer",
    default = 32,
    help = "Number of threads to use [default= %default]",
    metavar = "integer"
  ),
  make_option(
    c("--data_dir"),
    type = "character",
    help = "Directory containing BAM files",
    metavar = "character"
  ),
  make_option(
    c("--bam_pattern"),
    type = "character",
    default = "*.flnc.sorted.bam$",
    help = "Pattern to match BAM files [default= %default]",
    metavar = "character"
  ),
  make_option(
    c("--prefix"),
    type = "character",
    help = "Prefix for output files",
    metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

species <- opt$species
mode <- opt$mode
do_postprocess <- opt$do_postprocess
num_threads <- opt$num_threads
data_dir <- opt$data_dir
bam_pattern <- opt$bam_pattern
prefix <- opt$prefix

# setwd('/scratch/asabe/projects/lornash')
setDTthreads(num_threads)
bambu_dir <- 'data/transcriptome'
reads_bam <- list.files(data_dir, bam_pattern, recursive = TRUE, full.names = TRUE)

if (species == 'human') {
  
  annotation_gtf <- 'data/references/gencode.v47.annotation.gtf'
  genome_fasta <- 'data/references/GRCh38.p14.genome.fa'

} else if (species == 'mouse') {

  annotation_gtf <- 'data/references/gencode.vM36.annotation.gtf'
  genome_fasta <- 'data/references/GRCm39.genome.fa'

}

annotation_bambu <- prepareAnnotations(annotation_gtf)

cache_dir = str_glue("{bambu_dir}/cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

if (mode == 'quant') {

  se <- bambu(
    reads = reads_bam,
    annotations = annotation_bambu,
    genome = genome_fasta,
    rcOutDir = cache_dir,
    opt.discovery = list(
      min.readCount = 2,
      min.readFractionByGene = 0.05,
      min.sampleNumber = 2
    ),
    ncore = num_threads
  )
  
  output_prefix <- str_glue("{bambu_dir}/{prefix}")
  dir.create(output_prefix, recursive = TRUE, showWarnings = FALSE)
  saveRDS(se, str_glue("{output_prefix}.rds"))
  writeBambuOutput(se, output_prefix)

  if (do_postprocess) {
    bambu_counts_file <- str_glue("{output_prefix}/counts_transcript.txt")
    counts <- fread(bambu_counts_file)
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
    fwrite(counts, str_glue("{output_prefix}.counts.csv"), quote = FALSE)
  }

} else if (mode == 'discovery') {
  se_discovery <- bambu(
    reads = reads_bam,
    annotations = annotation_bambu,
    genome = genome_fasta,
    rcOutDir = 'data/bambu/cache',
    quant = FALSE,
    NDR = 1,
    ncore = num_threads
  )
  
  output_prefix <- str_glue("{bambu_dir}/databank_{species}_bambu_se_discovery")
  saveRDS(se_discovery, str_glue("{output_prefix}.rds"))
  writeToGTF(se_discovery, str_glue("{output_prefix}.gtf"))
}