setwd('/scratch/asabe/projects/foundation-model/downstream/clinvar')

library(data.table)
library(stringr)
library(magrittr)
library(ggplot2)
library(PRROC)

### Variant Metadata
metadata_file = '/scratch/asabe/projects/foundation-model/preprocess/pre-mrna/data/clinvar_2stars.preprocessed.updated.metadata.csv.gz'

meta = fread(metadata_file)
meta[, chr := str_remove(ucsc_chr, 'chr')]
meta[, var_id := str_c(chr, PositionVCF, ReferenceAlleleVCF, AlternateAlleleVCF, sep = '-')]
setnames(meta, 'seq_alt_id', 'mutant_id')
setkey(meta, transcript_id, mutant_id)

### Variant Positions
variant_summary = fread('/scratch/asabe/projects/foundation-model/downstream/clinvar/data/variant_summary.txt', select = c('#AlleleID', 'Name', 'Assembly'))
variant_summary = variant_summary[Assembly == 'GRCh38']

extract_numeric_position <- function(hgvs) {
  # Regular expression to capture the intronic position and direction
  matches <- str_match(hgvs, "c\\.(\\d+)([+-])(\\d+)([ACGT])")
  
  if (!is.na(matches[1])) {
    direction <- matches[3]
    
    if (direction == "+") {
      return(as.integer(matches[4]))  # Return the intronic position as integer
    } else if (direction == "-") {
      return(-as.integer(matches[4])) # Return the negative intronic position
    }
  }
  return(NA_integer_)  # Return NA if no match is found
}

colnames(variant_summary) = c('allele_id', 'hgvs', 'assembly')
variant_summary = variant_summary[, .(allele_id, hgvs)]

### LoRNA^{SH}

ll_files = list.files('data/likelihoods', '*likelihoods.csv', full.names = T)
ll_file_ids = ll_files %>% basename() %>% str_remove('[.].*')

.read_lls = function(ll_file, ll_file_id) {
  
  lls = fread(ll_file)
  
  lls[, var_id := str_split(transcript_id, '_', simplify = T) %>% extract(, 3)]
  trs = lls[var_id == '']
  trs = trs[, .(transcript_id, likelihood)]
  setnames(trs, 'likelihood', 'wild_ll')
  setkey(trs, transcript_id)
  
  vars = lls[var_id != '']
  vars = vars[, .(transcript_id, likelihood)]
  setnames(vars, c('transcript_id', 'likelihood'), c('mutant_id', 'mutant_ll'))
  vars[, transcript_id := str_extract(mutant_id, 'N._.*_') %>% str_remove('_$')]
  setkey(vars, transcript_id)
  
  vars = merge(vars, trs)
  vars[, diff_ll := mutant_ll - wild_ll]
  ll_cols = vars[, colnames(.SD), .SDcols = patterns('_ll')]
  setnames(vars, ll_cols, str_c(ll_file_id, ll_cols, sep = '_'))
  setkey(vars, transcript_id, mutant_id)
  
  vars
}

lls_ = mapply(.read_lls, ll_files, ll_file_ids, SIMPLIFY = F)
lls = merge(lls, meta, by = c('transcript_id', 'mutant_id'))

lornash = lls[, .(transcript_id, var_id, label, lornash_wild_ll, lornash_mutant_ll, lornash_diff_ll, AlleleID)]
hyenadna = lls[, .(transcript_id, var_id, label, hyenadna_wild_ll, hyenadna_mutant_ll, hyenadna_diff_ll, AlleleID)]
nuct = lls[, .(transcript_id, var_id, label, nuct_wild_ll, nuct_mutant_ll, nuct_diff_ll, AlleleID)]

### SpliceAI
vcf_file = 'spliceai/clinvar_20240331.non_coding.annotated.vcf'

vcf_data <- fread(vcf_file, skip = "#CHROM")

# Extract the INFO column and split it into individual fields
info_fields <- vcf_data[, tstrsplit(INFO, ";", fixed = TRUE)]

# Extract unique field names dynamically from the INFO fields
info_names <- unique(unlist(lapply(info_fields, function(x) sub("=.*", "", x))))

# Assign temporary column names to the split INFO fields
setnames(info_fields, paste0("INFO_", seq_along(info_fields)))
vcf_data <- cbind(vcf_data, info_fields)

# Convert the INFO fields into separate columns by extracting their values
vcf_data[, (info_names) := lapply(info_names, function(name) {
  sapply(1:.N, function(i) {
    # Extract the value for the specific INFO field
    info_list <- unlist(strsplit(INFO[i], ";"))
    value <- grep(paste0("^", name, "="), info_list, value = TRUE)
    ifelse(length(value) > 0, sub(".*=", "", value), NA_character_)
  })
})]

# Remove the original INFO fields
vcf_data[, grep("^INFO_", names(vcf_data)) := NULL]
vcf_data[, var_id := str_c(`#CHROM`, POS, REF, ALT, sep = '-')]
vcf = merge(vcf_data, meta, by = 'var_id')
spliceai_info = vcf$SpliceAI %>% str_split('[|]', n = 10, simplify = T)

spliceai_ds = spliceai_info[, c(3,4,5,6)]
spliceai_ds = spliceai_ds %>% apply(2, as.numeric)
spliceai_scores = rowSums(spliceai_ds)
vcf$score = spliceai_scores

vcf = vcf[!is.na(score)]
vcf[, ALLELEID := as.integer(ALLELEID)]

spliceai =  vcf[, .(transcript_id, var_id, label, score, ALLELEID)]

setnames(spliceai, c('score', 'ALLELEID'), c('spliceai_score', 'allele_id'))
scores = merge(lornash, spliceai, by = c('transcript_id', 'var_id', 'label', 'allele_id'))

scores = merge(scores, variant_summary, by = 'allele_id')
scores[, rel_pos := sapply(hgvs, extract_numeric_position)]
scores[, hgvs := NULL]

scores[, label_binary := fifelse(label == "pathogenic", 1, 0)]

spliceai_pr <- pr.curve(
  scores.class0 = scores[label_binary == 1, spliceai_score],
  scores.class1 = scores[label_binary == 0, spliceai_score],
  curve = TRUE
)

spliceai_roc <- roc.curve(
  scores.class0 = scores[label_binary == 1, spliceai_score],
  scores.class1 = scores[label_binary == 0, spliceai_score],
  curve = TRUE
)

plot(spliceai_pr, main = 'SpliceAI - PR')
plot(spliceai_roc, main = 'SpliceAI - ROC')

lornash_pr <- pr.curve(
  scores.class0 = scores[label_binary == 1, lornash_dll],
  scores.class1 = scores[label_binary == 0, lornash_dll],
  curve = TRUE
)

lornash_roc <- roc.curve(
  scores.class0 = scores[label_binary == 1, lornash_dll],
  scores.class1 = scores[label_binary == 0, lornash_dll],
  curve = TRUE
)

plot(lornash_pr, main = 'LoRNA^{SH} - PR')
plot(lornash_roc, main = 'LoRNA^{SH} - ROC')

spliceai_pr <- pr.curve(
  scores.class0 = scores[abs(rel_pos) > 5][label_binary == 1, spliceai_score],
  scores.class1 = scores[abs(rel_pos) > 5][label_binary == 0, spliceai_score],
  curve = TRUE
)

spliceai_roc <- roc.curve(
  scores.class0 = scores[abs(rel_pos) > 5][label_binary == 1, spliceai_score],
  scores.class1 = scores[abs(rel_pos) > 5][label_binary == 0, spliceai_score],
  curve = TRUE
)

plot(spliceai_pr, main = 'SpliceAI - PR - Excluding SS±5')
plot(spliceai_roc, main = 'SpliceAI - ROC - Excluding SS±5')

lornash_pr <- pr.curve(
  scores.class0 = scores[abs(rel_pos) > 5][label_binary == 1, lornash_dll],
  scores.class1 = scores[abs(rel_pos) > 5][label_binary == 0, lornash_dll],
  curve = TRUE
)

lornash_roc <- roc.curve(
  scores.class0 = scores[abs(rel_pos) > 5][label_binary == 1, lornash_dll],
  scores.class1 = scores[abs(rel_pos) > 5][label_binary == 0, lornash_dll],
  curve = TRUE
)

plot(lornash_pr, main = 'LoRNA^{SH} - PR - Excluding SS±5')
plot(lornash_roc, main = 'LoRNA^{SH} - ROC - Excluding SS±5')

spliceai_pr <- pr.curve(
  scores.class0 = scores[rel_pos < -5][label_binary == 1, spliceai_score],
  scores.class1 = scores[rel_pos < -5][label_binary == 0, spliceai_score],
  curve = TRUE
)

spliceai_roc <- roc.curve(
  scores.class0 = scores[rel_pos < -5][label_binary == 1, spliceai_score],
  scores.class1 = scores[rel_pos < -5][label_binary == 0, spliceai_score],
  curve = TRUE
)

plot(spliceai_pr, main = 'SpliceAI - PR - Donor Side Only - Excluding SS±5')
plot(spliceai_roc, main = 'SpliceAI - ROC - Donor Side Only - Excluding SS±5')

lornash_pr <- pr.curve(
  scores.class0 = scores[rel_pos < -5][label_binary == 1, lornash_dll],
  scores.class1 = scores[rel_pos < -5][label_binary == 0, lornash_dll],
  curve = TRUE
)

lornash_roc <- roc.curve(
  scores.class0 = scores[rel_pos < -5][label_binary == 1, lornash_dll],
  scores.class1 = scores[rel_pos < -5][label_binary == 0, lornash_dll],
  curve = TRUE
)

plot(lornash_pr, main = 'LoRNA^{SH} - PR - Donor Side Only - Excluding SS±5')
plot(lornash_roc, main = 'LoRNA^{SH} - ROC - Donor Side Only - Excluding SS±5')

spliceai_pr <- pr.curve(
  scores.class0 = scores[rel_pos > 5][label_binary == 1, spliceai_score],
  scores.class1 = scores[rel_pos > 5][label_binary == 0, spliceai_score],
  curve = TRUE
)

spliceai_roc <- roc.curve(
  scores.class0 = scores[rel_pos > 5][label_binary == 1, spliceai_score],
  scores.class1 = scores[rel_pos > 5][label_binary == 0, spliceai_score],
  curve = TRUE
)

plot(spliceai_pr, main = 'SpliceAI - PR - Acceptor Side Only - Excluding SS±5')
plot(spliceai_roc, main = 'SpliceAI - ROC - Acceptor Side Only - Excluding SS±5')

lornash_pr <- pr.curve(
  scores.class0 = scores[rel_pos > 5][label_binary == 1, lornash_dll],
  scores.class1 = scores[rel_pos > 5][label_binary == 0, lornash_dll],
  curve = TRUE
)

lornash_roc <- roc.curve(
  scores.class0 = scores[rel_pos > 5][label_binary == 1, lornash_dll],
  scores.class1 = scores[rel_pos > 5][label_binary == 0, lornash_dll],
  curve = TRUE
)

plot(lornash_pr, main = 'LoRNA^{SH} - PR - Acceptor Side Only - Excluding SS±5')
plot(lornash_roc, main = 'LoRNA^{SH} - ROC - Acceptor Side Only - Excluding SS±5')
