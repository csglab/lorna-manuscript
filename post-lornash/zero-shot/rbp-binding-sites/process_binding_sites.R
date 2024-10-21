setwd('/scratch/asabe/projects/foundation-model/downstream/narrowpeak')

library(data.table)
library(stringr)
library(stringi)
library(tidyr)
library(magrittr)

metadata_file = 'data/metadata.tsv'
events_ioe_file = 'events/c2l2_SE_strict.ioe'
suppa_psi_file = 'psi/c2l2.psi'
counts_file = 'quant/mpaqt.csv'
tr_seqeunces_csv = '/scratch/asabe/projects/foundation-model/preprocess/pre-mrna/data/databank_human_104_bambu_se_quant.all_chrs.preprocessed.updated.csv.gz'

# flank = 4096
ss_flank = 5

meta = fread(metadata_file)
meta[, .N, `File assembly`]
meta = meta[`File assembly` == 'GRCh38']
meta = meta[`Biological replicate(s)` == '1, 2']
meta = meta[`Library strand specific` == 'strand-specific']

files = meta[, `File accession`]

files = str_c('data/', files, '.bed.gz')

narrowPeak_columns = c('chromosome', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak')

beds = lapply(files, fread, col.names = narrowPeak_columns)
beds = rbindlist(beds)

name_splitted = beds[, name %>% str_split('_', simplify = T, n = 3)]

beds[, rbp := name_splitted[, 1]]
beds[, cell_line := name_splitted[, 2]]
beds[, .N, cell_line]

beds[cell_line == 'SM-9MVZL', cell_line := 'AdrenalGland']
beds = beds[cell_line != '']

#### Reading SUPPA2's generated SE events
events = fread(events_ioe_file)

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

psi[, avg := rowMeans(.SD, na.rm = T), .SDcols = patterns('_1|_2')]
psi[, total := rowSums(.SD, na.rm = T), .SDcols = patterns('_1|_2')]

psi = psi[total > 0]
psi = psi[, .(event_id, total)]

events = events[psi, on = .(event_id)]

### Getting major transcripts
counts = fread(counts_file)
counts = counts[, .SD, .SDcols = patterns('transcript|K562|HEPG2')]

counts[, HepG2 := rowMeans(.SD), .SDcols = patterns('HEPG2')]
counts[, K562 := rowMeans(.SD), .SDcols = patterns('K562')]

counts = counts[, .(transcript_id, HepG2, K562)]

counts = melt(counts, id.vars = 'transcript_id', measure.vars = c('HepG2', 'K562'), variable.name = 'cell_line', value.name = 'tpm')

event_trs = events[, .(event_id, gene_id, alternative_transcripts)] %>%
  separate_longer_delim(alternative_transcripts, delim = ',') %>%
  as.data.table() %>%
  unique()
setnames(event_trs, 'alternative_transcripts', 'transcript_id')

event_trs = merge(event_trs, counts, by = 'transcript_id')
major_trs = event_trs[, .SD[which.max(tpm)], .(event_id, gene_id, cell_line)]

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

setkey(event_coords, chr, strand, start, end)
setnames(beds, 'chromosome', 'chr')
setkey(beds, chr, strand, start, end)

beds_within = foverlaps(beds, event_coords, type = 'within', nomatch = NULL)
beds_within = merge(beds_within, major_trs, by = c('event_id', 'cell_line'))
beds_within = beds_within[tpm >= 1]

beds_within[, i.len := i.end - i.start + 1]
beds_within = beds_within[i.len >= 10]

find_available_ranges <- function(occupied, start, end) {
  occupied <- occupied[order(start)]
  available_ranges <- list()
  current_start <- start
  
  for (i in 1:nrow(occupied)) {

    if (current_start < occupied$start[i]) {
      available_ranges <- append(available_ranges, list(c(current_start, occupied$start[i] - 1)))
    }
    current_start <- max(current_start, occupied$end[i] + 1)
  }
  
  if (current_start <= end) {
    available_ranges <- append(available_ranges, list(c(current_start, end)))
  }
  
  return(available_ranges)
}

generate_random_regions <- function(group, start, end) {
  set.seed(42)
  occupied <- data.table(start = group$i.start, end = group$i.end)
  
  lengths <- occupied$end - occupied$start + 1
  num_regions <- length(lengths)
  
  available_ranges <- find_available_ranges(occupied, start, end)
  
  new_regions <- list()
  
  for (length in lengths) {
    valid_range_found <- FALSE
    
    available_ranges <- Filter(function(x) !is.null(x) && length(x) == 2, available_ranges)

    for (range in available_ranges) {
      range_start <- range[1]
      range_end <- range[2]
      
      if (!is.na(range_start) && !is.na(range_end) && (range_end - range_start + 1 >= length)) {
    
        new_start <- sample(range_start:(range_end - length + 1), 1)
        new_end <- new_start + length - 1
    
        new_regions <- append(new_regions, list(c(new_start, new_end)))

        updated_ranges <- list()
        for (r in available_ranges) {
          if (new_end < r[1] || new_start > r[2]) {
            updated_ranges <- append(updated_ranges, list(r))
          } else {
            if (new_start > r[1] && new_end < r[2]) {
              updated_ranges <- append(updated_ranges, list(c(r[1], new_start - 1), c(new_end + 1, r[2])))
            } else if (new_start > r[1]) {
              updated_ranges <- append(updated_ranges, list(c(r[1], new_start - 1)))
            } else if (new_end < r[2]) {
              updated_ranges <- append(updated_ranges, list(c(new_end + 1, r[2])))
            }
          }
        }
        
        available_ranges <- updated_ranges
        valid_range_found <- TRUE
        break
      }
    }
    if (!valid_range_found) {
      warning(paste("No valid range found for length", length, "in event", .BY$event_id, ", region", .BY$region, ", cell line", .BY$cell_line))
      break
    }
  }
  
  return(new_regions)
}

result <- beds_within[, .(new_regions = generate_random_regions(.SD, .BY$start, .BY$end)), 
                      by = .(event_id, region, cell_line, start, end)]

expanded_result <- result[, .(i.start_new = sapply(new_regions, `[`, 1), 
                              i.end_new = sapply(new_regions, `[`, 2)), 
                          by = .(event_id, region, cell_line, start, end)]

beds_unique = beds_within[, c('event_id', 'cell_line', 'chr', 'strand', 'region', 'start', 'end', 'gene_id', 'transcript_id'), with = F]
beds_unique = unique(beds_unique)

beds_within_expanded = merge(
  beds_unique,
  expanded_result,
  by = c('event_id', 'region', 'cell_line', 'start', 'end'),
  all.x = T)

setnames(beds_within_expanded, c('i.start_new', 'i.end_new'), c('i.start', 'i.end'))

cols_to_add = setdiff(colnames(beds_within), colnames(beds_within_expanded))
beds_within_expanded[, (cols_to_add) := NA]
beds_within_expanded[, type := 'non_rbp']
beds_within[, type := 'rbp']
setcolorder(beds_within_expanded, colnames(beds_within))

beds_within_all = rbind(beds_within, beds_within_expanded)
selected_events = beds_within_all[, .(event_id, chr, i.start, i.end, strand, region, cell_line, rbp, type, gene_id, transcript_id, tpm, signalValue, pValue)]
setnames(selected_events, c('i.start', 'i.end'), c('start', 'end'))
selected_events = selected_events[!(is.na(start) | is.na(end))]

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

selected_events[, seq_id := str_c(transcript_id, event_id, region, cell_line, rbp, type, rel_start, rel_end, sep = '@')]
selected_events[, seq_len := nchar(seq)]

fwrite(selected_events, 'sequences/narrowpeak.selected_events.csv.gz')
system('scp sequences/narrowpeak.selected_events.csv.gz chimera:/home/saberi/projects/lornash/data/')
