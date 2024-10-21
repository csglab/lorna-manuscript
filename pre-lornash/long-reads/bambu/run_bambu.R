library(optparse)
library(data.table)
library(stringr)
library(magrittr)
library(bambu)

option_list <- list(
  make_option(c("--species"), type = "character", default = "human", help = "Species: either 'human' or 'mouse' [default= %default]", metavar = "character")
  make_option(c("--mode"), type = "character", help = "Mode: either 'discovery' or 'quant'", metavar = "character"),
  make_option(c("--do_postprocess"), type = "logical", default = FALSE, help = "Perform post-processing [default= %default]", metavar = "logical"),
  make_option(c("--num_threads"), type = "integer", default = 32, help = "Number of threads to use [default= %default]", metavar = "integer"),
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

species <- opt$species
mode <- opt$mode
do_postprocess <- opt$do_postprocess
num_threads <- opt$num_threads

CWD <- '/scratch/asabe/projects/lornash'
bambu_dir <- 'data/bambu'

setwd(CWD)
setDTthreads(num_threads)

if (species == 'human') {
  annotation_gtf <- 'data/references/annotation/gencode.v46.annotation.gtf'
  genome_fasta <- 'data/references/genome/GRCh38.p14.genome.fa'
  reads_bam <- list.files('data/long-reads/bam', '*.flnc.sorted.bam$', recursive = TRUE, full.names = TRUE)
} else if (species == 'mouse') {
  annotation_gtf <- 'data/references/annotation/gencode.vM34.annotation.gtf'
  genome_fasta <- 'data/references/genome/GRCm39.genome.fa'
  reads_bam <- list.files(str_glue("data/databank/{species}"), '*.sorted.bam$', recursive = TRUE, full.names = TRUE)
}

annotation_bambu <- prepareAnnotations(annotation_gtf)

if (mode == 'quant') {

  se <- bambu(
    reads = reads_bam,
    annotations = annotation_bambu,
    genome = genome_fasta,
    rcOutDir = 'data/bambu/cache',
    opt.discovery = list(
      min.readCount = 2,
      min.readFractionByGene = 0.05,
      min.sampleNumber = 2
    ),
    ncore = num_threads
  )
  
  output_prefix <- str_glue("{bambu_dir}/databank_{species}_bambu_se")
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
    fwrite(counts, str_glue("{output_prefix}.counts.tsv"), sep = '\t', quote = FALSE)
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