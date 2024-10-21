setwd('/scratch/asabe/projects/foundation-model/downstream/splice-site-mutations')

library(data.table)
library(stringr)
library(stringi)
library(tidyr)
library(magrittr)

se_events_ioe_file = '/scratch/asabe/projects/foundation-model/downstream/narrowpeak/events/c2l2_SE_strict.ioe'
se_events_psi_file = '/scratch/asabe/projects/foundation-model/downstream/narrowpeak/psi/c2l2.psi'
counts_file = '/scratch/asabe/projects/foundation-model/downstream/narrowpeak/quant/mpaqt.csv'
tr_seqeunces_csv = '/scratch/asabe/projects/foundation-model/preprocess/pre-mrna/data/databank_human_104_bambu_se_quant.all_chrs.preprocessed.updated.csv.gz'

events = fread(se_events_ioe_file)

events[, num_alt_trs := str_count(alternative_transcripts, ',') + 1]
events[, num_total_trs := str_count(total_transcripts, ',') + 1]
events[, frac_alt_trs := round(num_alt_trs / num_total_trs, 2)]

id_splitted = events[, event_id %>%
                       str_remove('.*;SE:') %>%
                       str_split_fixed(':', 4)]

chromosome = id_splitted[, 1]

left_coords = id_splitted %>%
  extract(, 2) %>%
  str_split_fixed('-', 2)

right_coords = id_splitted %>%
  extract(, 3) %>%
  str_split_fixed('-', 2)

left_ss =
  left_coords %>%
  extract(, 2) %>%
  as.integer() %>%
  subtract(1)

right_ss =
  right_coords %>%
  extract(, 1) %>%
  as.integer() %>%
  add(1)

strand = id_splitted[, 4]

events[, chr := chromosome]
events[, left_ss := left_ss]
events[, right_ss := right_ss]
events[, strand := strand]

events[, acceptor := fifelse(strand == '+', left_ss, right_ss)]
events[, donor := fifelse(strand == '+', right_ss, left_ss)]

### Filtering events with no coverage
psi = fread(se_events_psi_file)
setnames(psi, 'V1', 'event_id')

psi[, total := rowSums(.SD, na.rm = T), .SDcols = patterns('_1|_2')]
psi = psi[total > 0]
psi = psi[, .(event_id, total)]

events = events[psi, on = .(event_id)]

### Getting major transcripts
counts = fread(counts_file)

counts[, total := rowMeans(.SD), .SDcols = patterns('_1|_2')]
counts[, tpm := total / sum(total) * 1e6]

counts = counts[, .(transcript_id, tpm)]

event_inclusion_trs = events[, .(event_id, gene_id, alternative_transcripts)] %>%
  separate_longer_delim(alternative_transcripts, delim = ',') %>%
  as.data.table() %>%
  unique()
setnames(event_inclusion_trs, 'alternative_transcripts', 'transcript_id')

event_all_trs = events[, .(event_id, gene_id, total_transcripts)] %>%
  separate_longer_delim(total_transcripts, delim = ',') %>%
  as.data.table() %>%
  unique()
setnames(event_all_trs, 'total_transcripts', 'transcript_id')

event_exclusion_trs = fsetdiff(event_all_trs, event_inclusion_trs)

event_inclusion_trs = merge(event_inclusion_trs, counts, by = 'transcript_id')
event_exclusion_trs = merge(event_exclusion_trs, counts, by = 'transcript_id')

major_events_inclusion = event_inclusion_trs[, .(total_tpm = sum(tpm)), event_id][total_tpm >= 1, event_id]
major_events_exclusion = event_exclusion_trs[, .(total_tpm = sum(tpm)), event_id][total_tpm >= 1, event_id]

major_events = intersect(major_events_inclusion, major_events_exclusion)
event_inclusion_trs = event_inclusion_trs[event_id %in% major_events]
event_exclusion_trs = event_exclusion_trs[event_id %in% major_events]

major_inclusion_trs = event_inclusion_trs[, .SD[which.max(tpm)], .(event_id, gene_id)]
major_exclusion_trs = event_exclusion_trs[, .SD[which.max(tpm)], .(event_id, gene_id)]

major_inclusion_trs[, type := 'inclusion']
major_exclusion_trs[, type := 'exclusion']

major_trs = rbind(major_inclusion_trs, major_exclusion_trs)

selected_events = merge(major_trs, events[, .(event_id, gene_id, chr, donor, acceptor, strand)],
      by = c('event_id', 'gene_id'))

selected_events = melt(selected_events, 
     measure.vars = c("acceptor", "donor"), 
     variable.name = "ss_type", 
     value.name = "pos")

#### Masking splice sites
seqs = fread(tr_seqeunces_csv)
seqs_coords = seqs[, .(gene_id, transcript_id, chr, strand, gene_start, gene_end)]

selected_events = merge(selected_events, seqs_coords, by = c('gene_id', 'transcript_id', 'chr', 'strand'))

selected_events[, rel_pos := fifelse(strand == '+',
                                     pos - gene_start + 1,
                                     gene_end - pos + 1 )]

selected_events[, rel_pos := rel_pos + 2]  ## +2 for H and S

selected_events[, rel_start := fcase(
  ss_type == 'donor' & strand == '+', rel_pos,
  ss_type == 'acceptor' & strand == '+', rel_pos - 1,
  ss_type == 'donor' & strand == '-', rel_pos,
  ss_type == 'acceptor' & strand == '-', rel_pos - 1,
  default = NA
)]

selected_events = selected_events[!is.na(rel_start)]
selected_events[, rel_end := fifelse(rel_start == rel_pos,
                                     rel_pos + 1,
                                     rel_pos)]

selected_events = merge(selected_events, seqs[, .(transcript_id, seq)], by = 'transcript_id')

selected_events[, `:=`(
  replaced_part = mapply(function(seq, start, end) {
    str_sub(seq, start, end)
  }, seq, rel_start, rel_end),
  
  seq = mapply(function(seq, start, end) {
    str_sub(seq, start, end) <- strrep("N", end - start + 1)
    seq
  }, seq, rel_start, rel_end)
)]

selected_events[, seq_id := str_c(transcript_id, event_id, type, ss_type, rel_pos, sep = '@')]
selected_events[, seq_len := nchar(seq)]

fwrite(selected_events, 'sequences/skipped_exons.splice_sites.csv.gz')

system('scp sequences/skipped_exons.splice_sites.csv.gz chimera:/home/saberi/projects/lornash/data/')
