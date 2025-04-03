if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "edgeR"))

library(data.table)
library(rtracklayer)
library(DESeq2)
library(edgeR)


RNAseq_counts <- read.csv("/data3/users/Qiuqoo/CM_ZJU/DEG_DMG/DEG/RNAseq_counts_1027.csv")


#### ###########

counts <- RNAseq_counts[,c("NSE1", "NSE2", "NSE3", "NSE4",
                           "GE1", "GE2", "GE3", "GE4",
                           "HE1", "HE2", "HE3",
                           "PYE1", "PYE3", "PYE4",
                           "LTE1", "LTE3", "LTE4",   
                           "LCE1", "LCE2", "LCE3", "LCE4"   ,
                           "PCE1", "PCE2", "PCE3", "PCE4")]
#################################

samples <- list("NSE" = c("NSE1", "NSE2", "NSE3", "NSE4"),
                "GE" = c("GE1", "GE2", "GE3", "GE4"),
                "HE" = c("HE1", "HE2", "HE3"),
                "PYE" = c("PYE1", "PYE3", "PYE4"),
                "LTE" = c("LTE1", "LTE3", "LTE4"),
                "PCE" = c("PCE1", "PCE2", "PCE3", "PCE4"),
                "LCE" = c("LCE1", "LCE2", "LCE3", "LCE4"))


for (sample in names(samples)) {

  sampleCondition <- factor(ifelse(colnames(counts) %in% samples[[sample]], "T", "C"),
                            levels = c("C", "T"))
  

  exp_count <- DGEList(counts = counts, group = sampleCondition, genes = rownames(counts))
  

  exp_count <- estimateDisp(exp_count)
  

  DEGs <- exactTest(exp_count, pair = c("C", "T"))
  

  idx <- which((abs(DEGs$table$logFC) >= 1) & DEGs$table$PValue < 0.01)
  DEGs_filtered <- DEGs[idx,]
  

  DEGs_filtered <- as.data.frame(DEGs_filtered$table)
  DEGs_filtered$ID <- rownames(DEGs_filtered)
  

  write.csv(DEGs_filtered, file = paste0("/data3/users/Qiuqoo/CM_ZJU/DEG_DMG/DEG/DEGs_", sample, "_logFC_1_P_0.01.csv"))
}
