setwd('/scratch/asabe/projects/foundation-model/downstream/exon-trap')

library(data.table)
library(stringr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(ggExtra)
library(fields)

exon_trap_counts_file = 'paper/data/data_zip/2D.txt'
exon_trap_likelihoods_file = 'data/lornash.exon_trap_dataset.likelihoods.csv'
exon_trap_likelihoods_metadata_file = 'data/exon_trap.preprocessed.updated.metadata.csv.gz'

counts = fread(exon_trap_counts_file)
lls = fread(exon_trap_likelihoods_file)

lls = separate(lls, exon_id, c('exon_id', 'vec_id', 'type'), sep = '_') %>% as.data.table()
lls = merge(lls, counts[, .(exon_id, genomic_category, gc_category, count)])

ggplot(lls) +
  aes(x = type, y = likelihood, fill = type) +
  geom_violin() +
  theme_bw() +
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text=element_text(size=12),
        axis.text.x = element_text(angle = 90)))

lls[, cor(likelihood, log(count)), .(type)]
lls[, cor(likelihood, log(count)), .(type, genomic_category)]

lls_exon = lls[type == 'exon', .(exon_id, likelihood)]
lls_intron = lls[type == 'intron']
lls_merged = merge(lls_exon, lls_intron, by = 'exon_id', suffixes = c('_exon', '_intron'))
lls_merged[, delta_ll := likelihood_exon - likelihood_intron]
lls_merged[, delta_tts := tts_logprobs_exon - tts_logprobs_intron]

lls_merged[, .(corr = cor(delta_ll, log(count)))]
lls_merged[, cor(delta_ll, log(count)), .(genomic_category)]

smoothScatter(lls_merged$delta_ll, log(lls_merged$count), 
              xlab = "Delta Log-Likelihood", 
              ylab = "Log(Counts)", 
              cex.axis = 1.5, 
              cex.lab = 1.7)
par(new = TRUE)
image.plot(legend.only = TRUE, 
           zlim = c(0, max(lls_merged$count, na.rm = TRUE)),
           col = colorRampPalette(c("white", "lightblue", "blue", "darkblue", "black"))(256),
           legend.width = 1.5, legend.mar = 5,
           axis.args = list(cex.axis = 1.5),
           legend.args = list(text = "", side = 4, font = 2, line = 3, cex = 1.7))

pearson_corr <- cor.test(lls_merged$delta_ll, log(lls_merged$count), method = "pearson")
spearman_corr <- cor.test(lls_merged$delta_ll, log(lls_merged$count), method = "spearman")
pearson_text <- paste0("Pearson : rho = ", round(pearson_corr$estimate, 2), 
                       ", p = ", format.pval(pearson_corr$p.value, digits = 2))
spearman_text <- paste0("Spearman: rho = ", round(spearman_corr$estimate, 2), 
                        ", p = ", format.pval(spearman_corr$p.value, digits = 2))
text(x = -33, y = 11.6, labels = pearson_text, pos = 4, col = "darkred", cex = 1.6)
text(x = -33, y = 11.3, labels = spearman_text, pos = 4, col = "darkred", cex = 1.6)

hyd = fread('data/hyenadna.exon_trap_dataset.likelihoods.csv')
hyd = merge(hyd, counts[, .(exon_id, count)], by = 'exon_id')
hyd[, cor(likelihood, log(count))]

nuct = fread('data/nuct.exon_trap_dataset.likelihoods.csv')
nuct = merge(nuct, counts[, .(exon_id, count)], by = 'exon_id')
nuct[, cor(likelihood, log(count))]

spliceai = fread('data/spliceai.exon_trap_dataset.likelihoods.csv')
spliceai[, likelihood := log(acceptor_prob) + log(donor_prob)] 
spliceai = merge(spliceai, counts[, .(exon_id, count)], by = 'exon_id')
spliceai[, cor(likelihood, log(count))]

ggplot(lls_merged[, .(corr = cor(delta_ll, log(count)))]) +
  aes(x = genomic_category, y = corr) +
  geom_col(alpha = 0.5, position='dodge') +
  theme_bw() +
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text=element_text(size=12)) +
  scale_y_continuous(name = 'Correlation: dLL vs. logCounts', limits = c(0, 0.2), breaks = seq(0, 0.2, 0.01))

(ggplot(lls_merged) +
  aes(x = delta_ll, y = log(count)) +
  geom_point() +
  geom_density_2d() +
  theme_bw() +
  geom_vline(xintercept = -50, linetype = 'dashed', color = 'red')) %>%
  ggMarginal(type = 'histogram', bins = 100)