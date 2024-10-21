setwd('/scratch/asabe/projects/foundation-model/downstream/logfc-classification')

library(data.table)
library(stringr)
library(magrittr)
library(ggplot2)

res_files = list.files('results', '.csv', full.names = T)
res = lapply(res_files, fread)
res = rbindlist(res)

res = res[Set == 'Test']
res = res[Metric == 'AUROC']

res = dcast(res, Tissue + Cell_Line_1 + Cell_Line_2 + Task + Set ~ Metric, value.var = "Value")

res[, cell_line_1 := fifelse(Task == 'Task 1', Cell_Line_1, Cell_Line_2)]
res[, cell_line_2 := fifelse(Task == 'Task 1', Cell_Line_2, Cell_Line_1)]
res[, tissue := Tissue]

res = res[, .(tissue, cell_line_1, cell_line_2, AUROC)]
setnames(res, 'AUROC', 'test_auroc')

ordered_cell_lines <- res[, .(cell_lines = unique(c(cell_line_1, cell_line_2))), by = tissue]
ordered_cell_lines <- ordered_cell_lines[order(tissue), cell_lines]

res[, cell_line_1 := factor(cell_line_1, levels = rev(ordered_cell_lines))]
res[, cell_line_2 := factor(cell_line_2, levels = ordered_cell_lines)]

res <- na.omit(res)

p <- ggplot(res, aes(x = cell_line_1, y = cell_line_2, fill = test_auroc)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(test_auroc, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "red", high = "green", mid = "yellow", midpoint = 0.75, limit = c(0.5, 1), space = "Lab", name = "AUROC") +
  theme_minimal() +
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Global Heatmap for All Tissues", x = "Cell Line 1", y = "Cell Line 2") +
  coord_fixed()