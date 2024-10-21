setwd('/scratch/asabe/projects/foundation-model/downstream/differential-analysis')

suppressMessages({
  
  library(data.table)
  library(stringr)
  library(magrittr)
  library(tidyr)
  library(ggplot2)  
  library(rtracklayer)
  library(BSgenome.Hsapiens.UCSC.hg38)
  
})

extract <- magrittr::extract

depmap_mutations = fread('/scratch/asabe/projects/foundation-model/downstream/depmap/data/DepMap_Mutation_Annotations.tsv')
depmap_mutations[, str_c(`locus annotation`, `sample id`, sep = '__') %>% unique() %>% length()]

depmap_cell_lines = fread('data/20240913T150850_RD01_DepMap_23Q4_sampleids_mappedto_samplenames.tsv')

depmap_mutations = merge(depmap_mutations, depmap_cell_lines, by = 'sample id')
depmap_sample2cell = depmap_mutations[, c('sample id', 'sample name', 'locus mutation', 'transcript id'), with = F] %>% unique()

c2l2_to_depmap = fread('data/c2l2_to_depmap.csv')
c2l2_to_depmap[depmap_name == '', depmap_name := NA]

depmap_mutations = merge(depmap_mutations, c2l2_to_depmap, by.x = 'sample name', by.y = 'depmap_name')

depmap_mutations = depmap_mutations[, .(tissue, cell_line_id, `gene id`, chrom, start, end, strand, ref, alt)]
depmap_mutations = unique(depmap_mutations)
setnames(depmap_mutations, 'gene id', 'gene_id')

depmap_mutations[, is_snv :=
                   (start == end) &
                   (nchar(ref) == 1) &
                   (nchar(alt) == 1)]

depmap_mutations = depmap_mutations[is_snv == T]

### 420 single nucleotide variants on c2l2 cell lines
fwrite(depmap_mutations, 'analysis/depmap_mutations.processed.csv')

depmap_mutations <- fread('analysis/depmap_mutations.processed.csv')
depmap_mutations[, pos := start]

splicing_event_files <- list.files(path = '/scratch/asabe/projects/foundation-model/downstream/narrowpeak/events', pattern = "*.ioe", full.names = T)
splicing_events <- lapply(splicing_event_files, fread)
splicing_events <- rbindlist(splicing_events, use.names = TRUE)

splicing_events[, event_type := sub(".*;([[:alnum:]_]+):.*", "\\1", event_id)]

splicing_events[, .N, event_type]

splicing_events <- splicing_events[grepl("Bambu", alternative_transcripts)]
splicing_events <- splicing_events[!grepl("ENST", alternative_transcripts)]


splicing_events[, `:=`(chrom = NA_character_, start = NA_real_, end = NA_real_)]

splicing_events[, c("chrom", "start", "end") := {

  if (event_type == "SE") {
  
    tstrsplit(event_id, "[:|-]", keep = c(2, 3, 6))
  } else if (event_type == "MX") {
  
  
    tstrsplit(event_id, "[:|-]", keep = c(2, 3, 10))
  } else if (event_type == "A5" | event_type == "A3") {
  
  
    tstrsplit(event_id, "[:|-]", keep = c(2, 3, 6))
  } else if (event_type == "RI") {
  
    tstrsplit(event_id, "[:|-]", keep = c(2, 3, 6))
  } else if (event_type == "AF" | event_type == "AL") {
  
    tstrsplit(event_id, "[:|-]", keep = c(2, 3, 8))
  } else {
  
    tstrsplit(event_id, "[:|-]", keep = c(2, 3, 6))
  }
}, event_id]

splicing_events[, `:=`(start = as.numeric(start),
                       end = as.numeric(end))]

depmap_mutations[, `:=`(start = pos, end = pos)]

setkey(depmap_mutations, chrom, start, end)
setkey(splicing_events, chrom, start, end)

overlapping_mutations <- foverlaps(depmap_mutations, splicing_events, nomatch = 0L)
overlapping_mutations = overlapping_mutations[, .(tissue, cell_line_id, chrom, pos, ref, alt, event_type, event_id, alternative_transcripts, total_transcripts)] %>% unique()

trs1 = overlapping_mutations[, .(alternative_transcripts, chrom, pos, ref, alt)] %>% 
  separate_rows(alternative_transcripts) %>% 
  as.data.table()

trs2 = overlapping_mutations[, .(total_transcripts, chrom, pos, ref, alt)] %>% 
  separate_rows(total_transcripts) %>% 
  as.data.table()

colnames(trs1) = colnames(trs2) = c('transcript_id', 'chrom', 'pos', 'ref', 'alt')
trs = rbind(trs1, trs2)



flank_len = 16
max_seq_len = 2^16 ## 64K

ref_genome = BSgenome.Hsapiens.UCSC.hg38
gtf_file = '../narrowpeak/references/c2l2.annotation.v2.gtf'
output_file = 'data/c2l2.depmap.overlapping_mutations.updated.csv.gz'
  
gtf = import(gtf_file)
gtf = as.data.table(gtf)

columns_to_keep = c('seqnames', 'type', 'feature', 'start', 'end', 'width', 'strand', 'gene_id', 'transcript_id', 'exon_number')
columns_to_keep = intersect(columns_to_keep, colnames(gtf))
gtf = gtf[, columns_to_keep, with = F]
gtf = gtf[type %in% c('transcript', 'exon')]
gtf[, exon_number := as.integer(exon_number)]

gtf = gtf[transcript_id %in% unique(trs$transcript_id)]
### Removing transcripts with 1 exon only
gtf[, num_exons := max(exon_number, na.rm = T), transcript_id]

gtf = gtf[num_exons > 1]

gtf[, gene_start := min(start), gene_id]
gtf[, gene_end := max(end), gene_id]

### Adding flanking regions to the gene coords
gtf[, gene_start := gene_start - flank_len]
gtf[, gene_end := gene_end + flank_len]

print(overlapping_mutations)

gtf[, tokenized_gene_len := (gene_end - gene_start + 1) + 1 + 1 + 1 + 1 + 1]
gtf = gtf[tokenized_gene_len <= max_seq_len]
gtf[, tokenized_gene_len := NULL]

### Reversing exon_number on the negative strand
gtf[strand == '-', exon_number := (num_exons - exon_number + 1)]

introns = gtf[type == 'exon']
introns = introns[, start_ := start]
introns[order(end), start := c(end[-.N] + 1, NA), transcript_id]
introns[order(end), end := c(start_[-1] - 1, NA), transcript_id]
introns[, start_ := NULL]
introns[, width := end - start + 1]
introns[, type := 'intron']
introns = introns[!is.na(start)]
gtf = rbind(gtf, introns)

non_rna = gtf[type == 'exon']
non_rna[, first_exon_num := min(exon_number), transcript_id]
non_rna[, last_exon_num := max(exon_number), transcript_id]
non_rna[, is_first_exon := fifelse(exon_number == first_exon_num, T, F)]
non_rna[, is_last_exon := fifelse(exon_number == last_exon_num, T, F)]
non_rna = non_rna[is_first_exon | is_last_exon]
non_rna[is_first_exon == T, end := start-1]
non_rna[is_first_exon == T, start := gene_start]
non_rna[is_last_exon == T, start := end+1]
non_rna[is_last_exon == T, end := gene_end]
non_rna[, width := end - start + 1]
non_rna[is_first_exon == T, type := 'non_rna_upstream']
non_rna[is_last_exon == T, type := 'non_rna_downstream']
non_rna[, c('first_exon_num', 'last_exon_num', 'is_first_exon', 'is_last_exon') := NULL]
gtf = rbind(gtf, non_rna)

gtf = gtf[strand != '*']

invalid_trs = gtf[(start <= 0) | (end <= 0) | (start > end), transcript_id]
print(str_glue('{length(invalid_trs)} invalid transcripts (due to negative coordinates)'))
gtf = gtf[!(transcript_id %in% invalid_trs)]

## if the strand is '-', it'll revComp. the sequence, which we don't want at this step.
gtf[type != 'transcript', seq := getSeq(ref_genome,
                                        names = seqnames,
                                        start = start,
                                        end = end,
                                        as.character = TRUE)]
# strand = strand


invalid_trs = gtf[str_count(seq, 'N') > 0, transcript_id %>% unique()]
print(str_glue('{length(invalid_trs)} invalid transcripts (due to N)'))
gtf = gtf[!(transcript_id %in% invalid_trs)]

non_rna_transformer = c('A' = 'W', 'C' = 'X', 'G' = 'Y', 'T' = 'Z')
transform_to_non_rna = function(seq) {
  non_rna_transformer[
    seq %>% str_split('') %>% unlist()
  ] %>% str_c(collapse = '')
}

transform_to_non_rna_v = Vectorize(transform_to_non_rna, 'seq')
gtf[type %like% 'non_rna', seq := transform_to_non_rna_v(seq)]

gtf[type == 'non_rna_upstream',  seq := str_c(seq, 'S')] ## Adding seq, S (TSS), will add the H/M token later
gtf[type == 'non_rna_downstream', seq := str_c('E', seq)] ## Adding E (TTS), flank_end
gtf[type == 'exon', seq := str_to_lower(seq)]

gtf[, type := factor(type, levels = c('non_rna_upstream', 'exon', 'intron', 'non_rna_downstream', 'transcript'), ordered = TRUE)]

tr_seqs = gtf[type != 'transcript'][
  order(transcript_id, exon_number, type),
  .(seq = str_c(seq, collapse = '')),
  by = .(transcript_id, strand)]

complement_transformer = c('A' = 'T', 'T' = 'A', 'C' = 'G', 'G' = 'C',
                           'a' = 't', 't' = 'a', 'c' = 'g', 'g' = 'c',
                           'W' = 'Z', 'Z' = 'W', 'X' = 'Y', 'Y' = 'X',
                           'D' = 'R', 'R' = 'D',
                           'S' = 'E', 'E' = 'S')

transform_to_reverse_complement = function(seq) {
  complement_transformer[
    seq %>% str_split('') %>% unlist() %>% rev()
  ] %>% str_c(collapse = '')
}

transform_to_reverse_complement_v = Vectorize(transform_to_reverse_complement, 'seq')

tr_seqs[strand == '-', seq := transform_to_reverse_complement_v(seq)]
species_token = 'H'
tr_seqs[, seq := str_c(species_token, seq)]

tr_seqs[, seq_len := nchar(seq)]
tr_seqs[, strand := NULL]

tr_info = gtf[type == 'transcript',
              c('seqnames', 'start', 'end', 'width', 'strand', 'gene_id', 'transcript_id', 'num_exons', 'gene_end', 'gene_start')] 

colnames(tr_info) = c('chr', 'tr_start', 'tr_end', 'tr_len', 'strand', 'gene_id', 'transcript_id', 'num_exons', 'gene_end', 'gene_start')

tr_seqs = merge(tr_seqs, tr_info, by = 'transcript_id', all.x = TRUE)

tr_seqs[, species := 'human']

setcolorder(tr_seqs, c('species', 'chr', 'strand', 'gene_start', 'gene_end', 'gene_id', 'tr_start', 'tr_end', 'transcript_id', 'num_exons', 'tr_len', 'seq_len', 'seq'))

setorderv(tr_seqs, c('species', 'chr', 'gene_start', 'gene_id', 'tr_start', 'transcript_id'))

seqs = tr_seqs 
seqs = seqs[species == 'human']
# 
seqs_ = merge(seqs, trs, by = 'transcript_id')
seqs_[, seq := str_remove(seq, '^H')]

complement_transformer = c('A' = 'T', 'T' = 'A', 'C' = 'G', 'G' = 'C',
                           'a' = 't', 't' = 'a', 'c' = 'g', 'g' = 'c',
                           'W' = 'Z', 'Z' = 'W', 'X' = 'Y', 'Y' = 'X',
                           'D' = 'R', 'R' = 'D',
                           'S' = 'E', 'E' = 'S')

transform_to_reverse_complement = function(seq) {
  complement_transformer[
    seq %>% str_split('') %>% unlist() %>% rev()
  ] %>% str_c(collapse = '')
}

transform_to_reverse_complement_v = Vectorize(transform_to_reverse_complement, 'seq')

seqs_[strand == '-', seq := transform_to_reverse_complement_v(seq)]

seqs_[, seq := str_c('H', seq)]

seqs_[, tss_offset := str_locate(seq, 'S') %>% extract(, 1)] ### TSS position
seqs_[, tts_offset := str_locate(seq, 'E') %>% extract(, 1)] ### TTS position
seqs_[, start := pos]
seqs_[, end := pos]
seqs_[, pos := NULL]
 
seqs_[, pos := fcase(
  start < gene_start + tss_offset, (start - gene_start + 1) + 1,
  start %between% list(gene_start + tss_offset, gene_start + tts_offset), (start - gene_start + 1) + 1 + 1,
  start > gene_start + tts_offset, (start - gene_start + 1) + 1 + 1 + 1,
  default = NA
)]

non_rna_transformer = c('A' = 'W', 'C' = 'X', 'G' = 'Y', 'T' = 'Z')

seqs_[, ref_tok := str_sub(seq, pos, pos)]
seqs_[, alt_tok := fcase(
  ref_tok %in% c('A', 'C', 'G', 'T'), alt,
  ref_tok %in% c('a', 'c', 'g', 't'), str_to_lower(alt),
  ref_tok %in% c('W', 'X', 'Y', 'Z'), non_rna_transformer[alt],
  default = NA
)]

# seqs_[, .(ref, alt, ref_tok, alt_tok)]

seqs_[, seq_alt := seq]
str_sub(seqs_$seq_alt, seqs_$pos, seqs_$pos) <- seqs_$alt_tok

seqs_[, seq_alt_id := str_c(gene_id, transcript_id, pos, ref_tok, alt_tok, chr, start, end, strand, ref, alt, sep = '@')]
seqs_[, seq_id := str_c(gene_id, transcript_id, pos, ref_tok, ref_tok, chr, start, end, strand, ref, ref, sep = '@')]

seqs_orig = seqs_[, .(seq_id, chr, strand, seq_len, seq)] %>% unique()
seqs_mut = seqs_[, .(seq_alt_id, chr, strand, seq_len, seq_alt)] %>% unique()
colnames(seqs_orig) = colnames(seqs_mut) = c('seq_id', 'chr', 'strand', 'seq_len', 'seq')

seqs__ = rbind(seqs_orig, seqs_mut)
seqs__[strand == '-', seq := seq %>% str_remove('^H') %>% transform_to_reverse_complement_v() %>% str_c('H', .)]

fwrite(seqs__, 'analysis/c2l2.depmap.overlapping_mutations.updated.csv.gz')
system('scp analysis/c2l2.depmap.overlapping_mutations.updated.csv.gz chimera:/home/saberi/projects/lornash/data/')
