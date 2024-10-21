setwd('/scratch/asabe/projects/foundation-model/downstream/abundance')

library(data.table)
library(stringr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(cocor)
library(patchwork)


### short-read
# counts_files = list.files('/home/asabe/scratch/projects/miniencode-dev/data/short-reads/isoform-abundance', '.abundance.tsv', full.names = T)

# counts_ = lapply(counts_files, fread, select = c('target_id', 'tpm'))
# .merge = function(x, y) merge(x, y, by = 'target_id', suffix = c(str_c(sample(LETTERS, 4),collapse=''), str_c(sample(LETTERS, 4),collapse='')))
# counts = Reduce(.merge, counts_)

# counts[, read_count := rowMeans(.SD), .SDcols = patterns('tpm')]
# counts = counts[, .(target_id, read_count)]
# setnames(counts, 'target_id', 'tr_id')

### long-reads
# counts_file = 'data/miniencode.count_data.csv'
# row_data_file = 'data/miniencode.row_data.csv'

# counts = fread(counts_file)
# counts = counts[, .(TXNAME, GENEID)]
# colnames(counts) = c('tr_id', 'gene_id')
# counts = counts[is.na(readCount), readCount := 0]
# counts = counts[readCount == 0, relReadCount := 0]
# counts = counts[, .(TXNAME, GENEID, NDR, readCount, relReadCount)]
# colnames(counts) = c('tr_id', 'gene_id', 'ndr', 'read_count', 'rel_read_count')
# setkey(counts, tr_id)

### mpaqt
counts_file = 'data/mpaqt.updated.csv'
row_data_file = 'data/miniencode.row_data.csv'

counts = fread(counts_file)
row_data = fread(row_data_file)

counts[, read_count := rowMeans(.SD), .SDcols = patterns('_1$|_2$')]

counts = merge(counts, row_data, by = 'transcript_id')
setnames(counts, 'transcript_id', 'tr_id')
setkey(counts, tr_id)

ll_files = list.files('data', '*likelihoods.csv', full.names = T)
ll_file_ids = ll_files %>% basename() %>% str_remove('[.].*')

.read_lls = function(ll_file, ll_file_id) {
  lls = fread(ll_file, select = c('transcript_id', 'seq_len', 'tts_logprobs'))
  colnames(lls) = c('tr_id', str_c(ll_file_id, c('seq_len', 'tts_logprobs'), sep = '_'))
  setkey(lls, tr_id)
  lls
}

lls_ = mapply(.read_lls, ll_files, ll_file_ids, SIMPLIFY = F)
lls = Reduce(merge, lls_)

lapply(counts_files, function(x) {
  reads = fread(x, select = c('target_id', 'tpm'))
  temp = merge(reads, lls, by.x = 'target_id', 'tpm')
  })

count_lls = merge(counts, lls)

count_lls[, .SD %>% as.matrix() %>% cor(), .SDcols = patterns('_tts_logprobs')]
count_lls[, num_trs := .N, gene_id]

lornash_cor = count_lls[, cor(lornash_tts_logprobs, log(read_count))]
hyenadna_cor = count_lls[, cor(hyenadna_tts_logprobs, log(read_count))]
nuct_cor = count_lls[, cor(nuct_tts_logprobs, log(read_count))]

lornash_gene_cors = count_lls[num_trs > 2][, .(corr = cor(lornash_tts_logprobs, log(read_count))), gene_id][, corr]
hyenadna_gene_cors = count_lls[num_trs > 2][, .(corr = cor(hyenadna_tts_logprobs, log(read_count))), gene_id][, corr]
nuct_gene_cors = count_lls[num_trs > 2][, .(corr = cor(nuct_tts_logprobs, log(read_count))), gene_id][, corr]

bar_data <- data.frame(
  Method = factor(c("LoRNA^{SH}", "HyenaDNA", "NucT"), levels = c("LoRNA^{SH}", "HyenaDNA", "NucT")),
  Correlation = c(lornash_cor, hyenadna_cor, nuct_cor)
)

plot_bar <- ggplot(bar_data, aes(x = Method, y = Correlation, fill = Method)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("LoRNA^{SH}" = "#DA6959", "HyenaDNA" = "#B8B8B8", "NucT" = "#B8B8B8")) +
  scale_y_continuous(n.breaks = 5) +
  theme_minimal() +
  theme(
    legend.position = "none",                        
    axis.title.x = element_blank(),                  
    axis.title.y = element_blank(),                  
    axis.text.x = element_text(size = 20, angle = 90),           
    axis.text.y = element_text(size = 20),           
    axis.ticks.x = element_blank(),                  
    axis.ticks.length.y = unit(0, "pt"),             
    panel.grid.minor = element_blank(),              
    plot.title = element_text(hjust = 0.5, size = 20)
  )

violin_data <- data.frame(
  Correlation = c(lornash_gene_cors, hyenadna_gene_cors, nuct_gene_cors),
  Method = factor(rep(c("LoRNA^{SH}", "HyenaDNA", "NucT"),
                      c(length(lornash_gene_cors), length(hyenadna_gene_cors), length(nuct_gene_cors))))
)

plot_violin <- ggplot(violin_data, aes(x = Method, y = Correlation, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = c("LoRNA^{SH}" = "#DA6959", "HyenaDNA" = "#B8B8B8", "NucT" = "#B8B8B8")) +
  theme_minimal() +
  theme(
    legend.position = "none",                        
    axis.title.x = element_blank(),                  
    axis.title.y = element_blank(),                  
    axis.text.x = element_text(size = 20),           
    axis.text.y = element_text(size = 20),           
    axis.ticks.x = element_blank(),                  
    axis.ticks.length.y = unit(0, "pt"),             
    panel.grid.minor = element_blank(),              
    plot.title = element_text(hjust = 0.5, size = 20)
  ) + ylim(-1, 1)


boxplot_data <- data.frame(
  Correlation = c(lornash_gene_cors, hyenadna_gene_cors, nuct_gene_cors),
  Method = factor(rep(c("LoRNA^{SH}", "HyenaDNA", "NucT"),
                      c(length(lornash_gene_cors), length(hyenadna_gene_cors), length(nuct_gene_cors))),
                  levels = c("LoRNA^{SH}", "HyenaDNA", "NucT")) 
)

p_value_hyena <- wilcox.test(lornash_gene_cors, hyenadna_gene_cors)$p.value
p_value_nuct <- wilcox.test(lornash_gene_cors, nuct_gene_cors)$p.value

plot_boxplot <- ggplot(boxplot_data, aes(x = Method, y = Correlation, fill = Method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = c("LoRNA^{SH}" = "#DA6959", "HyenaDNA" = "#B8B8B8", "NucT" = "#B8B8B8")) +
  theme_minimal() +
  theme(
    legend.position = "none",                        
    axis.title.x = element_blank(),                  
    axis.title.y = element_blank(),                  
    axis.text.x = element_text(size = 20, angle = 90),           
    axis.text.y = element_text(size = 20),           
    axis.ticks.x = element_blank(),                  
    axis.ticks.length.y = unit(0, "pt"),             
    panel.grid.minor = element_blank(),              
    plot.title = element_text(hjust = 0.5, size = 20)
  ) + ylim(-1, 1.35) + 
  stat_compare_means(comparisons = list(c("LoRNA^{SH}", "HyenaDNA"), c("LoRNA^{SH}", "NucT")),
                     method = "wilcox.test", label = "p.format", vjust = -0.5)


bar_data <- data.frame(
  Method = factor(c("LoRNA^{SH}", "HyenaDNA", "NucT"), levels = c("LoRNA^{SH}", "HyenaDNA", "NucT")),
  Correlation = c(lornash_cor, hyenadna_cor, nuct_cor)
)

plot_bar <- ggplot(bar_data, aes(x = Method, y = Correlation, fill = Method)) +
  geom_bar(stat = "identity", width = 0.75, alpha = 0.7, colour="black") +
  scale_fill_manual(values = c("LoRNA^{SH}" = "#DA6959", "HyenaDNA" = "#B8B8B8", "NucT" = "#B8B8B8")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.ticks.x = element_blank(),
    axis.ticks.length.y = unit(0, "pt"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) + ylim(0.0, 0.16)# 

boxplot_data <- data.frame(
  Correlation = c(lornash_gene_cors, hyenadna_gene_cors, nuct_gene_cors),
  Method = factor(rep(c("LoRNA^{SH}", "HyenaDNA", "NucT"),
                      c(length(lornash_gene_cors), length(hyenadna_gene_cors), length(nuct_gene_cors))),
                  levels = c("LoRNA^{SH}", "HyenaDNA", "NucT"))
)

p_value_hyena <- wilcox.test(lornash_gene_cors, hyenadna_gene_cors)$p.value
p_value_nuct <- wilcox.test(lornash_gene_cors, nuct_gene_cors)$p.value

plot_boxplot <- ggplot(boxplot_data, aes(x = Method, y = Correlation, fill = Method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = c("LoRNA^{SH}" = "#DA6959", "HyenaDNA" = "#B8B8B8", "NucT" = "#B8B8B8")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.ticks.x = element_blank(),
    axis.ticks.length.y = unit(0, "pt"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) + ylim(-1, 1) +
  stat_compare_means(comparisons = list(c("LoRNA^{SH}", "HyenaDNA"), c("LoRNA^{SH}", "NucT")),
                     method = "wilcox.test", label = "p.format", vjust = -0.5)

combined_plot <- plot_bar + plot_boxplot + 
  plot_layout(ncol = 2, widths = c(1, 1))