setwd('/scratch/asabe/projects/foundation-model/downstream/motif-discovery')

library(data.table)
library(stringr)
library(tidyr)
library(magrittr)
extract <- magrittr::extract

seqs <- fread('sequences/skipped_exons.motif_discovery.csv', select = c('seq_id', 'replaced_part'))

# Read the data
md_lls <- fread('lls/se_md_dataset.likelihoods.csv')
md_lls <- md_lls %>% 
  separate(seq_id, c('transcript_id', 'event_id', 'region', 'reg_start', 'reg_end'), sep = '@', remove = F) %>% 
  as.data.table()

ref_trs <- fread('lls/downstream_trs_dataset.likelihoods.csv')
ref_trs <- ref_trs[, .(transcript_id, likelihood, tts_logprobs)]
setnames(ref_trs, c('likelihood', 'tts_logprobs'), c('ref_likelihood', 'ref_tts_logprobs'))

# Merge datasets and compute differences
md_lls <- merge(md_lls, ref_trs, by = 'transcript_id')
md_lls[, dll := ref_likelihood - likelihood]
md_lls[, summary(dll)]
md_lls[, dtts := ref_tts_logprobs - tts_logprobs]

md_lls <- md_lls[, .(seq_id, dll)]
md_lls <- merge(md_lls, seqs, by = 'seq_id')

md_lls = md_lls[, .(dll, replaced_part)]

fwrite(md_lls, 'lls/md_lls.csv')