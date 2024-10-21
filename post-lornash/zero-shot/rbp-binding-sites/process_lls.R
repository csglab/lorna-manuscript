setwd('/scratch/asabe/projects/foundation-model/downstream/narrowpeak')

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

np_lls <- fread('data/narrowpeak_dataset.likelihoods.csv')
np_lls <- np_lls %>% 
  separate(seq_id, c('transcript_id', 'event_id', 'region', 'cell_line', 'rbp', 'type', 'rel_start', 'rel_end'), sep = '@') %>% 
  as.data.table()

ref_trs <- fread('data/downstream_trs_dataset.likelihoods.csv')
ref_trs <- ref_trs[, .(transcript_id, likelihood, tts_logprobs)]
setnames(ref_trs, c('likelihood', 'tts_logprobs'), c('ref_likelihood', 'ref_tts_logprobs'))

np_lls <- merge(np_lls, ref_trs, by = 'transcript_id')
np_lls[, dll := ref_likelihood - likelihood]
np_lls[, dtts := ref_tts_logprobs - tts_logprobs]
np_lls[, Type := fifelse(type == 'rbp', 'RBP', 'NON_RBP')]

rbp_order <- np_lls[, .N, by = rbp][order(-N)]$rbp
np_lls[, rbp := factor(rbp, levels = rbp_order)]

np_lls[, Type2 := fifelse(Type == 'NON_RBP', 'NON_RBP', cell_line)]
np_lls[, Type2 := factor(Type2, levels = c('NON_RBP', 'HepG2', 'K562'), ordered = T)]
np_lls[, num_trs := .N, rbp]

theme_publication <- function() {
  theme_minimal() +
    theme(
      text = element_text(size = 14),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
}

p1 <- ggplot(np_lls, aes(x = Type, y = dll)) +
  geom_boxplot(outlier.shape = NA) +
  coord_flip() +
  labs(x = 'Type', y = 'DLL', title = 'Comparison of DLL by Type') +
  theme_publication() 

p2 <- ggplot(np_lls, aes(x = Type2, y = dll)) +
  geom_boxplot(outlier.shape = NA) +
  coord_flip() +
  labs(x = 'Type', y = 'DLL', title = 'Comparison of DLL by Type and Cell Line') +
  theme_publication() +
  theme(legend.position = 'bottom')

p3 <- ggplot(np_lls[rbp!='NON_RBP'], aes(x = reorder(rbp, num_trs), y = dll, fill = rbp)) +
  geom_boxplot(outlier.shape = NA) +
  coord_flip() +
  labs(x = 'Type', y = 'DLL', title = 'Comparison of DLL by Type and RBP') +
  theme_publication() +
  guides(fill = 'none')

p4 <- np_lls[rbp!='NON_RBP'] %>%
  group_by(rbp) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = reorder(rbp, count), y = count, fill = rbp)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  labs(x = 'RBP', y = 'Number of Transcript IDs', title = 'Number of Transcript IDs per RBP') +
  theme_publication() +
  theme(legend.position = 'none')

print(p1)
print(p2)
print(p3)
print(p4)

np_lls[, summary(dll), rbp]

result <- np_lls[, .(
  min_dll = min(dll, na.rm = TRUE),
  Q1_dll = quantile(dll, 0.25, na.rm = TRUE),
  med_dll = median(dll, na.rm = TRUE),
  mean_dll = mean(dll, na.rm = TRUE),
  Q3_dll = quantile(dll, 0.75, na.rm = TRUE),
  max_dll = max(dll, na.rm = TRUE)
), by = rbp]

fwrite(result, 'rbps_deletions_dlls_summary.csv')

p3 <- ggplot(np_lls, aes(x = reorder(rbp, num_trs), y = dll, fill = rbp)) +
  geom_boxplot(outlier.shape = NA) +
  coord_flip() +
  labs(x = 'RBP', y = 'DLL', title = 'Comparison of DLL by Type and RBP') +
  theme_publication() +
  guides(fill = 'none')

p4 <- np_lls %>%
  group_by(rbp) %>%
  summarise(count = n()) %>%
  as.data.table() %>% 
  ggplot(aes(x = reorder(rbp, count), y = count, fill = rbp)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  labs(x = 'RBP', y = 'Number of Transcript IDs', title = 'Number of Transcript IDs per RBP') +
  theme_publication() +
  theme(legend.position = 'none')

combined_plot <- p3 * p4 

print(combined_plot)
