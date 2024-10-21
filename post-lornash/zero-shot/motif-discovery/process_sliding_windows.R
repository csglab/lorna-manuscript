setwd('/scratch/asabe/projects/foundation-model/downstream/motif-discovery')

library(data.table)
library(stringr)
library(tidyr)
library(magrittr)

events_ioe_file = '/scratch/asabe/projects/foundation-model/downstream/narrowpeak/events/c2l2_SE_strict.ioe'
suppa_psi_file = '/scratch/asabe/projects/foundation-model/downstream/narrowpeak/psi/c2l2.psi'
counts_file = '/scratch/asabe/projects/foundation-model/downstream/narrowpeak/quant/mpaqt.csv'
tr_seqeunces_csv = '/scratch/asabe/projects/foundation-model/preprocess/pre-mrna/data/databank_human_104_bambu_se_quant.all_chrs.preprocessed.updated.csv.gz'

# flank = 4096
ss_flank = 5
intron_flank = 256

#### Reading SUPPA2's generated SE events
events = fread(events_ioe_file)

events[, num_alt_trs := str_count(alternative_transcripts, ',') + 1]
events[, num_total_trs := str_count(total_transcripts, ',') + 1]
events[, frac_alt_trs := round(num_alt_trs / num_total_trs, 2)]

# events[, summary(frac_alt_trs)]

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

left_intron_start =
  left_coords %>%
  extract(, 1) %>%
  as.integer() %>%
  add(ss_flank + 1)

left_intron_end =
  left_coords %>%
  extract(, 2) %>%
  as.integer() %>%
  subtract(ss_flank + 1)

left_intron_start = pmax(left_intron_start, left_intron_end - intron_flank)

middle_exon_start =
  left_coords %>%
  extract(, 2) %>%
  as.integer() %>%
  add(ss_flank)

middle_exon_end =
  right_coords %>%
  extract(, 1) %>%
  as.integer() %>%
  subtract(ss_flank)

right_intron_start =
  right_coords %>%
  extract(, 1) %>%
  as.integer() %>%
  add(ss_flank + 1)

right_intron_end =
  right_coords %>%
  extract(, 2) %>%
  as.integer() %>%
  subtract(ss_flank + 1)

right_intron_end = pmin(right_intron_end, right_intron_start + intron_flank)

strand = id_splitted[, 4]

events[, chr := chromosome]
events[, left_intron_start := left_intron_start]
events[, left_intron_end := left_intron_end]
events[, middle_exon_start := middle_exon_start]
events[, middle_exon_end := middle_exon_end]
events[, right_intron_start := right_intron_start]
events[, right_intron_end := right_intron_end]
events[, strand := strand]

### Filtering events with no coverage
psi = fread(suppa_psi_file)
setnames(psi, 'V1', 'event_id')

# psi[, avg := rowMeans(.SD, na.rm = T), .SDcols = patterns('_1|_2')]
psi[, total := rowSums(.SD, na.rm = T), .SDcols = patterns('_1|_2')]

psi = psi[total > 0]
psi = psi[, .(event_id, total)]

events = events[psi, on = .(event_id)]

### Getting major transcripts
counts = fread(counts_file)
counts[, total := rowMeans(.SD), .SDcols = patterns('_1|_2')]
counts[, tpm := total / sum(total) * 1e6]

counts = counts[, .(transcript_id, tpm)]

# event_inclusion_trs = events[, .(event_id, gene_id, alternative_transcripts)] %>%
#   separate_longer_delim(alternative_transcripts, delim = ',') %>%
#   as.data.table() %>%
#   unique()
# setnames(event_inclusion_trs, 'alternative_transcripts', 'transcript_id')
# 
# event_all_trs = events[, .(event_id, gene_id, total_transcripts)] %>%
#   separate_longer_delim(total_transcripts, delim = ',') %>%
#   as.data.table() %>%
#   unique()
# setnames(event_all_trs, 'total_transcripts', 'transcript_id')
# 
# event_exclusion_trs = fsetdiff(event_all_trs, event_inclusion_trs)
# 
# event_inclusion_trs = merge(event_inclusion_trs, counts, by = 'transcript_id')
# event_exclusion_trs = merge(event_exclusion_trs, counts, by = 'transcript_id')
# 
# major_events_inclusion = event_inclusion_trs[, .(total_tpm = sum(tpm)), event_id][total_tpm >= 1, event_id]
# major_events_exclusion = event_exclusion_trs[, .(total_tpm = sum(tpm)), event_id][total_tpm >= 1, event_id]
# 
# major_events = intersect(major_events_inclusion, major_events_exclusion)
# event_inclusion_trs = event_inclusion_trs[event_id %in% major_events]
# event_exclusion_trs = event_exclusion_trs[event_id %in% major_events]
# 
# major_inclusion_trs = event_inclusion_trs[, .SD[which.max(tpm)], .(event_id, gene_id)]
# major_exclusion_trs = event_exclusion_trs[, .SD[which.max(tpm)], .(event_id, gene_id)]
# 
# major_inclusion_trs[, type := 'inclusion']
# major_exclusion_trs[, type := 'exclusion']
# 
# major_trs = rbind(major_inclusion_trs, major_exclusion_trs)
# 
# major_trs = major_trs[type == 'inclusion']

event_trs = events[, .(event_id, gene_id, alternative_transcripts)] %>%
  separate_longer_delim(alternative_transcripts, delim = ',') %>%
  as.data.table() %>%
  unique()
setnames(event_trs, 'alternative_transcripts', 'transcript_id')

event_trs = merge(event_trs, counts, by = 'transcript_id')
major_trs = event_trs[, .SD[which.max(tpm)], .(event_id, gene_id)][tpm >= 1]

event_coords = events[, .SD, .SDcols = patterns('event_id|_start|_end|chr|strand')]
regions = c('left_intron', 'middle_exon', 'right_intron')
region_coords = lapply(regions, function(region) {
  region_coords = event_coords[, c('event_id', str_c(region, c('start', 'end'), sep = '_'), 'chr', 'strand'), with = F]
  region_coords[, region := region]
  setnames(region_coords, str_c(region, c('start', 'end'), sep = '_'), c('start', 'end'))
  setcolorder(region_coords, c('event_id', 'region', 'chr', 'start', 'end', 'strand'))
  region_coords
})

event_coords = rbindlist(region_coords)
event_coords = event_coords[start < end]


selected_events = merge(major_trs, event_coords,
                        by = 'event_id')

seqs = fread(tr_seqeunces_csv)
seqs_coords = seqs[, .(gene_id, transcript_id, chr, strand, gene_start, gene_end)]

selected_events = merge(selected_events, seqs_coords, by = c('gene_id', 'transcript_id', 'chr', 'strand'))

selected_events[, rel_start := fifelse(strand == '+',
                                       start - gene_start + 1,
                                       gene_end - end + 1)]

selected_events[, rel_end := fifelse(strand == '+',
                                     end - gene_start + 1,
                                     gene_end - start + 1)]

selected_events[, rel_start := rel_start + 2]  ## +2 for H and S
selected_events[, rel_end := rel_end + 2]  ## +2 for H and S

region_len = 20
region_flank = 10
region_coords = selected_events[, mapply(function(event_id, region, start, end, flank, len) {
  reg_starts = seq(start, end, flank)
  reg_starts = reg_starts[reg_starts <= end - len]
  reg_ends = reg_starts + len
  list(event_id, region, start, end, reg_starts, reg_ends)
  }, event_id, region, rel_start, rel_end, region_flank, region_len, SIMPLIFY = F)]

region_coords = lapply(region_coords, function(x)
  data.table(event_id = x[[1]],
             region = x[[2]],
             rel_start = x[[3]],
             rel_end = x[[4]],
             reg_start = x[[5]],
             reg_end = x[[6]]))

region_coords = rbindlist(region_coords)

selected_events = merge(selected_events, region_coords, by = c('event_id', 'region', 'rel_start', 'rel_end'))

selected_events = merge(selected_events, seqs[, .(transcript_id, seq)], by = 'transcript_id')

selected_events[, `:=`(
  replaced_part = mapply(function(seq, start, end) {
    str_sub(seq, start, end)
  }, seq, reg_start, reg_end),
  
  seq = mapply(function(seq, start, end) {
    str_sub(seq, start, end) <- strrep("N", end - start + 1)
    seq
  }, seq, reg_start, reg_end)
)]

selected_events[, seq_id := str_c(transcript_id, event_id, region, reg_start, reg_end, sep = '@')]
selected_events[, seq_len := nchar(seq)]

fwrite(selected_events, 'sequences/skipped_exons.motif_discovery.csv.gz')
system('scp sequences/skipped_exons.motif_discovery.csv.gz chimera:/home/saberi/projects/lornash/data/')

