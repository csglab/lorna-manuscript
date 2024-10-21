setwd('/scratch/asabe/projects/foundation-model/downstream/logfc-classification')

library(data.table)
library(stringr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(caret)
library(pROC)
library(PRROC)
library(doParallel)

num_cores <- 64

extract <- magrittr::extract

embeddings_file = '/scratch/asabe/projects/foundation-model/downstream/embeddings/miniencode.all_chrs.tts_original_embeddings.layer_15.csv'

embeds = fread(embeddings_file, header = T)
colnames(embeds) = c('transcript_id', str_c('dim', seq(128)))

### col_data
col_data = fread('/scratch/asabe/projects/foundation-model/downstream/differential-analysis/data/miniencode.col_data.csv')
col_data[, deseq_id := cell_line %>%
           str_replace_all('-', '_')  %>%
           str_replace_all(' ', '_') %>%
           str_replace_all('[.]', '_') %>%
           str_remove_all('[(]|[)]')]


unique_tissues <- unique(col_data$tissue)

for (tis in unique_tissues) {
  
  print(tis)
  cell_lines <- unique(col_data[tissue == tis, deseq_id])
  
  if (length(cell_lines) < 2)
    next()
  
  cell_line_pairs <- combn(cell_lines, 2, simplify = FALSE)
  
  for (pair in cell_line_pairs) {

    cell_line_1 <- col_data[deseq_id == pair[1], cell_line][1]
    cell_line_2 <- col_data[deseq_id == pair[2], cell_line][1]
    
    de_file <- paste0("/scratch/asabe/projects/foundation-model/downstream/differential-analysis/data/by_cell_line/", pair[1], "_vs_", pair[2], "/", pair[1], "_vs_", pair[2], "_DET_by_DESeq2_TRX_support_level_0_or_more.csv")
    
    if (file.exists(de_file)) {
      pdf_output = paste0("AUROC_AUPRC_", tis, "_", cell_line_1, "_vs_", cell_line_2, ".pdf")
      if (file.exists(pdf_output)) {
        print(str_glue("SKIPPING {pdf_output}"))
        next()
      }
      
      de_data = fread(de_file)
      setnames(de_data, 'TRX.ID', 'transcript_id', skip_absent = T)
      
      de_data = merge(de_data, embeds, by = 'transcript_id')
      
      cl <- makeCluster(num_cores)
      registerDoParallel(cl)

      de_data[, label1 := ifelse(log2FoldChange > 1 & padj < 0.05, "Class1", "Class0")]
      de_data[, label2 := ifelse(log2FoldChange < -1 & padj < 0.05, "Class1", "Class0")]

      de_data$label1 <- factor(de_data$label1, levels = c("Class0", "Class1"))
      de_data$label2 <- factor(de_data$label2, levels = c("Class0", "Class1"))

      X <- de_data[, .SD, .SDcols = patterns("^dim")]
      y1 <- de_data$label1
      y2 <- de_data$label2

      control <- trainControl(method = "cv", number = 10, summaryFunction = twoClassSummary, classProbs = TRUE, allowParallel = TRUE)

      model1 <- train(X, y1, method = "rf", trControl = control, metric = "ROC")
      pred1 <- predict(model1, X, type = "prob")

      conf_matrix1 <- confusionMatrix(predict(model1, X), y1)
      print(conf_matrix1)

      roc1 <- roc(response = y1, predictor = pred1[, "Class1"])
      auc1 <- auc(roc1)
      pr1 <- pr.curve(scores.class0 = pred1[y1 == "Class1", "Class1"], scores.class1 = pred1[y1 == "Class0", "Class1"], curve = TRUE)
      
      cat("AUROC for Task 1:", auc1, "\n")
      cat("AUPRC for Task 1:", pr1$auc.integral, "\n")

      model2 <- train(X, y2, method = "rf", trControl = control, metric = "ROC")
      pred2 <- predict(model2, X, type = "prob")

      conf_matrix2 <- confusionMatrix(predict(model2, X), y2)
      print(conf_matrix2)

      roc2 <- roc(response = y2, predictor = pred2[, "Class1"])
      auc2 <- auc(roc2)
      pr2 <- pr.curve(scores.class0 = pred2[y2 == "Class1", "Class1"], scores.class1 = pred2[y2 == "Class0", "Class1"], curve = TRUE)
      
      cat("AUROC for Task 2:", auc2, "\n")
      cat("AUPRC for Task 2:", pr2$auc.integral, "\n")

      pdf(pdf_output)

      plot(roc1, main = paste("AUROC for Task 1: ", tis, cell_line_1, " vs ", cell_line_2))

      plot(pr1, main = paste("AUPRC for Task 1: ", tis, cell_line_1, " vs ", cell_line_2), auc.main = TRUE)

      draw_confusion_matrix <- function(cm, title) {
        cm_table <- as.table(cm$table)
        fourfoldplot(cm_table, color = c("#CC6666", "#99CC99"),
                     conf.level = 0, margin = 1, main = title)
      }
      draw_confusion_matrix(conf_matrix1, paste("Confusion Matrix for Task 1: ", tis, cell_line_1, " vs ", cell_line_2))

      plot(roc2, main = paste("AUROC for Task 2: ", tis, cell_line_1, " vs ", cell_line_2))

      plot(pr2, main = paste("AUPRC for Task 2: ", tis, cell_line_1, " vs ", cell_line_2), auc.main = TRUE)

      draw_confusion_matrix(conf_matrix2, paste("Confusion Matrix for Task 2: ", tis, cell_line_1, " vs ", cell_line_2))
      
      dev.off()

      stopCluster(cl)
      registerDoSEQ()
      
    }
  }
}
