setwd('/scratch/asabe/projects/foundation-model/downstream/differential-analysis')

library(data.table)
library(stringr)
library(magrittr)
library(tidyr)
extract <- magrittr::extract

lls = fread('analysis/c2l2_depmap_dataset.likelihoods.csv')

lls = lls %>%
  separate(seq_id, c('gene_id', 'transcript_id', 'pos', 'ref_tok', 'alt_tok', 'chr', 'start', 'end', 'strand', 'ref', 'alt'), sep = '@') %>%
  as.data.table()

ref_lls = lls[ref == alt, .(gene_id, transcript_id, likelihood, tss_logprobs, tts_logprobs)] %>% unique()
setnames(ref_lls, c('likelihood', 'tss_logprobs', 'tts_logprobs'), c('ref_ll', 'ref_tss_ll', 'ref_tts_ll'))

mut_lls = lls[ref != alt]
mut_lls[, var_id := str_c(chr, start, end, ref, alt, sep = '_')]
mut_lls = mut_lls[, .(gene_id, transcript_id, var_id, likelihood, tss_logprobs, tts_logprobs)] %>% unique()
setnames(mut_lls, c('likelihood', 'tss_logprobs', 'tts_logprobs'), c('mut_ll', 'mut_tss_ll', 'mut_tts_ll'))

lls_ = merge(mut_lls, ref_lls, by = c('gene_id', 'transcript_id'))
lls_[, dll := mut_ll - ref_ll]
lls_[, tss_dll := mut_tss_ll - ref_tss_ll]
lls_[, tts_dll := mut_tts_ll - ref_tts_ll]

depmap_muts = fread('analysis/depmap_mutations.processed.csv')
depmap_muts[, var_id := str_c(chrom, start, end, ref, alt, sep = '_')]
depmap_muts = depmap_muts[, .(tissue, cell_line_id, var_id)] %>% unique()
dat = merge(lls_, depmap_muts, by = 'var_id')

mpaqt_tpm = fread('analysis/mpaqt.tpm.csv')
mpaqt_tpm = melt(mpaqt_tpm, id.vars = 'transcript_id', variable.name = 'cell_line_id', value.name = 'tpm')

dat_ = merge(dat, mpaqt_tpm, by = c('transcript_id', 'cell_line_id'), all.x = T)
overlapping_mutations[, var_id := str_c(chrom, pos, pos, ref, alt, sep = '_')]
dat_ = merge(dat_, overlapping_mutations, by = c('tissue', 'cell_line_id', 'var_id'), allow.cartesian = T)
dat_ = dat_[str_detect(alternative_transcripts, transcript_id)]

psi = lapply(list.files('../narrowpeak/psi', pattern = '*.psi', full.names = T), fread)
psi = lapply(psi, setnames, 'V1', 'event_id')
psi = rbindlist(psi)
psi = psi[event_id %in% unique(dat_$event_id)]
psi = melt(psi, id.vars = 'event_id', variable.name = 'cell_line_id', value.name = 'psi')

dat_ = merge(dat_, psi, by = c('cell_line_id', 'event_id'))

col_data = fread('../abundance/data/miniencode.col_data.csv')
col_data = col_data[, .(tissue, cell_line_id)]
col_data = unique(col_data)
cl_pairs = col_data[, expand.grid(cell_line_id, cell_line_id), tissue]
cl_pairs = cl_pairs[Var1 != Var2]
setnames(cl_pairs, c('Var1', 'Var2'), c('cell_line_id', 'cell_line_2'))

dat_ = merge(dat_, cl_pairs, by = c('tissue', 'cell_line_id'), allow.cartesian = T)
dat_ = merge(dat_, mpaqt_tpm, by.x = c('transcript_id', 'cell_line_2'), by.y =  c('transcript_id', 'cell_line_id'), all.x = T)
dat_ = merge(dat_, psi, by.x = c('event_id', 'cell_line_2'), by.y =  c('event_id', 'cell_line_id'))

setnames(dat_, c('cell_line_id', 'tpm.x', 'psi.x'), c('cell_line_1', 'tpm_1', 'psi_1'))
setnames(dat_, c('tpm.y', 'psi.y'), c('tpm_2', 'psi_2'))

dat_[, dtpm := tpm_1 - tpm_2]
dat_[, dpsi := psi_1 - psi_2]

setcolorder(dat_, c('tissue', 'cell_line_1', 'cell_line_2', 'event_type', 'dtpm', 'dpsi', 'dll', 'tpm_1', 'tpm_2', 'psi_1', 'psi_2', 'mut_ll', 'ref_ll', 'event_id', 'var_id'))

normalize_vector <- function(x) {
  return (2 * (x - min(x)) / (max(x) - min(x)) - 1)
}

dat_[, dll := normalize_vector(dll)]

pearson_cor <- cor(dat_$dpsi, dat_$dll, method = "pearson")
plot(dat_$dpsi, dat_$dll, xlim = c(-1, 1), ylim = c(-1, 1), 
     xlab = "dpsi", ylab = "dll", main = "Scatter Plot with Trend Line")
abline(lm(dll ~ dpsi, data = dat_), col = "blue")
text(x = 0.5, y = -0.8, labels = paste("Pearson:", round(pearson_cor, 2)))

pearson_cor <- cor.test(dat_$dpsi, dat_$dll, method = "pearson")
p_value <-  pearson_cor$p.value
pearson_cor <- pearson_cor$estimate

ggplot(dat_, aes(x = dpsi, y = dll, color = transcript_id, shape = event_type)) +
  geom_point(size = 5, stroke = 2) +
  geom_smooth(method = "lm", color = "red", se = FALSE, size = 1, aes(group = 1)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  xlim(-1, 1) +ylim(-1, 1) +
  labs(title = "dPSI (C2L2+MPAQT+SUPPA) vs. dLL (LoRNA^{SH})",
       x = expression(Delta * PSI), y = expression(Delta * LL), color = "Transcript ID", shape = "Event Type") +
  annotate("text", x = -0.9, y = 0.9, label = paste("Pearson: ", round(pearson_cor, 2)), 
           size = 5, color = "black", hjust = 0) +
  annotate("text", x = -0.9, y = 0.8, label = paste("P-value: ", round(p_value, 4)), 
           size = 5, color = "black", hjust = 0) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 16),               
        axis.title = element_text(size = 16)) +
  scale_shape_manual(values = c(2, 6, 0, 5, 1))
