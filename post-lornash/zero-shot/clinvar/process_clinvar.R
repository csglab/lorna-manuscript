setwd('/scratch/asabe/projects/foundation-model/downstream/clinvar')

library(data.table)
library(stringr)
library(magrittr)

variant_summary_file = 'variant_summary.txt'
assembly_report_file = 'GRCh38_latest_assembly_report.txt'
annotation_gff_file = 'GRCh38_latest_genomic.gtf'

tr_seqs_file = '/scratch/asabe/projects/foundation-model/preprocess/pre-mrna/data/refseq.GRCh38_latest_genomic.preprocessed.updated.csv'
tr_seqs_output_file = '/scratch/asabe/projects/foundation-model/preprocess/pre-mrna/data/clinvar_2stars.preprocessed.updated.csv.gz'
vars_metadata_output_file = '/scratch/asabe/projects/foundation-model/preprocess/pre-mrna/data/clinvar_2stars.preprocessed.updated.metadata.csv.gz'

chr_info = fread(assembly_report_file, skip = 63, select = c('RefSeq-Accn', 'UCSC-style-name'), col.names = c('refseq_chr', 'ucsc_chr'))
chr_info = chr_info[refseq_chr != 'na']

tr_info = fread(annotation_gff_file, skip = '#', header = F,
                col.names = c('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'))

tr_info[, transcript_id := attribute %>% str_extract('transcript_id "[^"]+"') %>% str_remove('transcript_id ') %>% str_remove_all('"')]
num_exons = tr_info[feature == 'exon', .(num_exons = .N), transcript_id]
tr_info = tr_info[feature == 'transcript']
tr_info = tr_info[, .(seqname, start, end, strand, transcript_id)]
tr_info = merge(tr_info, num_exons, by = 'transcript_id')
 
# # saveRDS(tr_info, 'process_clinvar.tr_info.rds')
# tr_info = readRDS('process_clinvar.tr_info.rds')

vars = fread(variant_summary_file)
setnames(vars, '#AlleleID', 'AlleleID')
vars = vars[Assembly == 'GRCh38']
vars = vars[ClinSigSimple != -1]
vars = merge(vars, chr_info, by.x = 'ChromosomeAccession', by.y = 'refseq_chr')
vars[, transcript_id := Name %>% str_split_fixed(':', n = 2) %>% extract(, 1) %>% str_remove(GeneSymbol) %>% str_remove('[(][)]')]
vars[, is_coding := transcript_id %>% like('^NM_|^XM_')]

vars[, .N, is_coding]
vars[, .N, Type]

vars = vars[is_coding == TRUE]

snps = vars[Type == 'single nucleotide variant']

snps = snps[, .(ChromosomeAccession, AlleleID, Name, ClinicalSignificance, PhenotypeList, OriginSimple, ReviewStatus, NumberSubmitters,
                VariationID, PositionVCF, ReferenceAlleleVCF, AlternateAlleleVCF, ucsc_chr, transcript_id)]

snps = merge(snps, tr_info, by = 'transcript_id')

min_num_exons = 2
max_seq_len = 2^16
snps = snps[num_exons >= min_num_exons]

### from https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
# snps[, .N, ReviewStatus]
snps[, star := fcase(
  ReviewStatus == 'practice guideline', 4,
  ReviewStatus == 'reviewed by expert panel', 3,
  ReviewStatus == 'criteria provided, multiple submitters, no conflicts', 2,
  ReviewStatus %in% c('criteria provided, conflicting classifications', 'criteria provided, single submitter'), 1,
  ReviewStatus %in% c('no assertion criteria provided', 'no classification provided', 'no classifications from unflagged records'), 0)]

# snps[, .N, star][order(star)]
min_clinvar_stars = 2
snps = snps[star >= min_clinvar_stars]

# snps[, .N, ClinicalSignificance]
snps[, label := fcase(
  ClinicalSignificance %like% 'athogenic', 'pathogenic',
  ClinicalSignificance %like% 'ncertain', 'vus',
  ClinicalSignificance %like% 'enign', 'benign',
  default = NA_character_
)]
# snps = snps[label == 'benign' | label == 'pathogenic']
snps = snps[!is.na(label)]

snps[, .N, label]
snps[, transcript_id %>% unique() %>% length()]

tr_ids = fread(tr_seqs_file, header = T, select = 'transcript_id')

snps = snps[tr_ids, on = .(transcript_id), nomatch = NULL]
snps = snps[ChromosomeAccession == seqname]

tr_seqs = fread(tr_seqs_file, header = T, select = c('transcript_id', 'gene_id', 'seq_len', 'chr', 'gene_start', 'gene_end', 'seq'))
tr_seqs = merge(tr_seqs, snps, by = 'transcript_id')

tr_seqs[, start_offset := str_locate(seq, 'S') %>% extract(, 1)] ### TSS position
tr_seqs[, pos := fifelse(strand == '+',
                         start_offset + (PositionVCF - start + 1),
                         start_offset + (end - PositionVCF + 1)
                         )]
tr_seqs = tr_seqs[start_offset < pos] ## variants on transcript sequence only, not on flanking regions

# ## Or: +1 for H and +1 for S
# tr_seqs[, pos := fifelse(strand == '+',
#                          (PositionVCF - gene_start + 1) + 1 + 1,
#                          (gene_end - PositionVCF + 1)  + 1 + 1)
#         ]

## transcripts' up/downstream variants
tr_seqs = tr_seqs[pos > 0]

# tr_seqs[nchar(ReferenceAlleleVCF) != 1]
# tr_seqs[nchar(AlternateAlleleVCF) != 1]

tr_seqs[, ref := str_sub(seq, pos, pos)]

# tr_seqs[str_to_upper(ref) != ReferenceAlleleVCF, .N, strand]
# complement = c(A = 'T', T = 'A', C = 'G', G = 'C')
# tr_seqs[strand == '+'][str_to_upper(ref) != ReferenceAlleleVCF]
# tr_seqs[strand == '-'][complement[str_to_upper(ref)] != ReferenceAlleleVCF]

tr_seqs[, alt := AlternateAlleleVCF]
tr_seqs[, alt := fifelse(ref == str_to_upper(ref), alt, str_to_lower(alt))]

complement = c(A = 'T', T = 'A', C = 'G', G = 'C',
               a = 't', t = 'a', c = 'g', g = 'c')

tr_seqs[strand == '-', alt := complement[alt]]

tr_seqs[, seq_alt := seq]

str_sub(tr_seqs$seq_alt, tr_seqs$pos, tr_seqs$pos) <- tr_seqs$alt

tr_seqs[, seq_id := transcript_id]
tr_seqs[, seq_alt_id := str_c(transcript_id, VariationID, sep = '_')]

variation_info = tr_seqs[, -c('seq', 'seq_alt')]

columns_to_keep = c('chr', 'gene_id', 'seq_len')
ref_seqs = tr_seqs[, c('seq_id', 'seq', columns_to_keep), with = F] %>% unique()
alt_seqs = tr_seqs[, c('seq_alt_id', 'seq_alt', columns_to_keep), with = F] %>% unique()

setnames(ref_seqs, 'seq_id', 'transcript_id')
colnames(alt_seqs) = colnames(ref_seqs) 
all_seqs = rbind(ref_seqs, alt_seqs)

fwrite(all_seqs, tr_seqs_output_file)
fwrite(variation_info, vars_metadata_output_file)
