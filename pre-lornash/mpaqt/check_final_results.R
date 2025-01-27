# setwd('/scratch/asabe/csg/mpaqt/lornash')

library(data.table)
library(stringr)
library(magrittr)

quant_files = list.files('data/transcriptome', '.csv', full.names = T)
quants = lapply(quant_files, fread)
# names(quants) = 
quant_ids = quant_files %>% basename() %>% str_remove('.*.extended_annotations.') %>% str_remove('.csv')

names(quants) = quant_ids


long_only = setdiff(quants$long_reads.cpm$transcript_id, quants$mpaqt.tpm$transcript_id)
mpaqt_only = setdiff(quants$mpaqt.tpm$transcript_id, quants$long_reads.cpm$transcript_id)


long_only[long_only %like% 'Bambu']
length(long_only)

gtf = rtracklayer::import('data/transcriptome/human.c2l2.revio.gencode_v47.GRCh38_p14/human.c2l2.revio.gencode_v47.GRCh38_p14.extended_annotations.sorted.biotyped.gtf')
gtf = as.data.table(gtf)

long_only_trs = gtf[type == 'transcript'][.(long_only), on = .(transcript_id)]

long_only_trs[, .N, gene_type]
long_only_trs[, .N, transcript_type]

all_trs = gtf[type == 'transcript'][, .(seqnames, start, end, width, strand, transcript_id, gene_id, gene_type, transcript_type)]

long_only_trs = long_only_trs[, .(seqnames, start, end, width, strand, transcript_id, gene_id, gene_type, transcript_type)]

all_trs[, .N]
quants$mpaqt.tpm[, .N]

fwrite(all_trs, 'data/transcriptome/human.c2l2.revio.gencode_v47.GRCh38_p14.extended_annotations.all_transcripts.csv.gz')

fwrite(long_only_trs, 'data/transcriptome/human.c2l2.revio.gencode_v47.GRCh38_p14.extended_annotations.excluded_in_mpaqt.csv.gz')
