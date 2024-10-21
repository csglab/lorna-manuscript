setwd('/scratch/asabe/projects/foundation-model/downstream/clinvar')

library(data.table)
library(stringr)
library(magrittr)
library(pROC)
library(ggplot2)
library(PRROC)

metadata_file = '/scratch/asabe/projects/foundation-model/preprocess/pre-mrna/data/clinvar_2stars.preprocessed.updated.metadata.csv.gz'
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
lls = Reduce(merge, lls_)

metadata = fread(metadata_file)
metadata = metadata[, .(seq_id, seq_alt_id, is_coding, label)]
setnames(metadata, c('seq_id', 'seq_alt_id'), c('transcript_id', 'mutant_id'))
setkey(metadata, transcript_id, mutant_id)

lls = merge(lls, metadata)

data = lls
data[, label_binary := fifelse(label.x == "pathogenic", 1, 0)]

ggplot(data) +
  aes(x = label.x, y = lornash_diff_ll, fill = label.x) +
  geom_boxplot()

lornash_roc <- roc(data$label_binary, data$lornash_diff_ll, plot=T)
lornash_auc <- auc(lornash_roc)

hyenadna_roc <- roc(data$label_binary, data$hyenadna_diff_ll)
hyenadna_auc <- auc(hyenadna_roc)

nuct_roc <- roc(data$label_binary, data$nuct_diff_ll)
nuct_roc <- auc(nuct_roc)

lornash_pr = pr.curve(scores.class0 = data[label_binary == 1, lornash_diff_ll],
                     scores.class1 = data[label_binary == 0, lornash_diff_ll],
                     curve = TRUE)

lornash_roc <- roc.curve(
  scores.class0 = data[label_binary == 1, lornash_diff_ll],
  scores.class1 = data[label_binary == 0, lornash_diff_ll],
  curve = TRUE
)

lornash_aucpr = lornash_pr$auc.integral

hyenadna_pr = pr.curve(scores.class0 = data[label_binary == 1, hyenadna_diff_ll],
                     scores.class1 = data[label_binary == 0, hyenadna_diff_ll],
                     curve = TRUE)

hyenadna_roc <- roc.curve(
  scores.class0 = data[label_binary == 1, hyenadna_diff_ll],
  scores.class1 = data[label_binary == 0, hyenadna_diff_ll],
  curve = TRUE
)

hyenadna_aucpr = hyenadna_pr$auc.integral

nuct_pr = pr.curve(scores.class0 = data[label_binary == 1, nuct_diff_ll],
                     scores.class1 = data[label_binary == 0, nuct_diff_ll],
                     curve = TRUE)

nuct_roc <- roc.curve(
  scores.class0 = data[label_binary == 1, nuct_diff_ll],
  scores.class1 = data[label_binary == 0, nuct_diff_ll],
  curve = TRUE
)

nuct_aucpr = nuct_pr$auc.integral

lornash_curve <- data.frame(Recall = lornash_pr$curve[, 1], Precision = lornash_pr$curve[, 2], Model = "lornash")
hyenadna_curve <- data.frame(Recall = hyenadna_pr$curve[, 1], Precision = hyenadna_pr$curve[, 2], Model = "hyenadna")
nuct_curve <- data.frame(Recall = nuct_pr$curve[, 1], Precision = nuct_pr$curve[, 2], Model = "nuct")
combined_curve <- rbind(lornash_curve, hyenadna_curve, nuct_curve)

pr_plot <- ggplot(combined_curve, aes(x = Recall, y = Precision, color = Model)) +
  geom_line(size = 1) +
  labs(title = "Precision-Recall Curve", x = "Recall", y = "Precision") +
  theme_minimal() +
  scale_color_manual(values = c("lornash" = "#DA6959", "hyenadna" = "#B8B8B8", "nuct" = "#70AD47"))

