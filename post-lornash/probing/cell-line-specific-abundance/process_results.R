setwd('/scratch/asabe/projects/foundation-model/downstream/fine-tune/abundance-prediction')

library(data.table)
library(stringr)
library(magrittr)
library(ggplot2)
library(ggpubr)

result_files = list.files('results', '.csv', full.names = T)
res_ = lapply(result_files, fread)
res = rbindlist(res_)
res = res[, .(cell_line, pearson_test, pvalue_test, model)]
colnames(res) = c('cell_line', 'test_pearson_cor', 'test_pvalue', 'model')

col_data = fread('/scratch/asabe/projects/pacbio/data/bambu/databank_human_104_bambu_se/miniencode.data.updated/miniencode.col_data.csv')

col_data = col_data[, .(tissue, cell_line, cell_line_id)]
setnames(res, 'cell_line', 'cell_line_id')
res = merge(res, col_data, by = 'cell_line_id')
res = unique(res)
res[tissue %like% 'Kidney', tissue := 'Kidney']

tissue_order <- c("Skin", "Lung", "Breast", "Liver", "Kidney", "Pancreas", "Colon", "Prostate", "Bone")
tissue_order <- res[order(-test_pearson_cor)][, unique(tissue)]

color_palette <- c(
  "Skin" = "#5D6EAB",
  "Lung" = "#5E8731",
  "Breast" = "#A54747",
  "Liver" = "#835C94",
  "Kidney" = "#C29D4C",
  "Pancreas" = "#964F4F",
  "Colon" = "#35644C",
  "Prostate" = "#A65F8A",
  "Bone" = "#D2802E"      
)
res[, tissue := factor(tissue, levels = tissue_order)]
setorder(res, tissue)

ggbarplot(res,
          x = "cell_line", y = "test_pearson_cor",
          fill = "tissue",             
          color = "white",          
          palette = color_palette,          
          sort.val = "desc",         
          sort.by.groups = T,    
          x.text.angle = 90,
          xlab = 'Cell Line',
          ylab = 'Pearson Correlation Coefficient',
          legend.title = 'Tissue',
          width = 0.95
)
