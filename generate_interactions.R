
# 加载必要包
library(Hmisc)
library(dplyr)


generate_interactions <- function(fullExpr, targetExpr, threshold_r = 0.8, threshold_p = 0.01, topN = NULL) {

  common_samples <- intersect(colnames(fullExpr), colnames(targetExpr))
  if(length(common_samples) == 0) {
    stop("fullExpr and targetExpr do not have common samples, please check the data!")
  }
  

  fullExpr <- fullExpr[, common_samples, drop = FALSE]
  targetExpr <- targetExpr[, common_samples, drop = FALSE]
  

  fullExpr <- fullExpr[apply(fullExpr, 1, sd) > 0, ]
  

  cor_res <- rcorr(as.matrix(t(fullExpr)), as.matrix(t(targetExpr)), type = "pearson")
  cor_matrix <- cor_res$r
  p_matrix   <- cor_res$P
  

  target_genes <- rownames(targetExpr)
  

  coexpressed_genes <- list()
  

  for (gene in target_genes) {

    gene_cor <- cor_matrix[, gene]
    gene_p   <- p_matrix[, gene]
    

    selected <- names(which(abs(gene_cor) > threshold_r & gene_p < threshold_p))
    selected <- setdiff(selected, target_genes)
    
    if(length(selected) > 0) {
      gene_interactions <- data.frame(
        TargetGene      = gene,
        CoexpressedGene = selected,
        Correlation     = gene_cor[selected],
        PValue          = gene_p[selected],
        stringsAsFactors = FALSE
      )
      coexpressed_genes[[gene]] <- gene_interactions
    }
  }
  
  all_interactions <- do.call(rbind, coexpressed_genes)
  all_interactions <- unique(all_interactions)
  
  if(!is.null(topN)) {
    sorted_interactions <- all_interactions[order(-abs(all_interactions$Correlation)), ]
    all_interactions <- head(sorted_interactions, topN)
  }
  
  return(all_interactions)
}





fullExpr <- read.csv("/data3/users/Qiuqoo/CM_ZJU/article/GRN/cluster_out_FPKM.csv", row.names = 1)
WOXExpr <- read.csv("/data3/users/Qiuqoo/CM_ZJU/article/GRN/WOX.csv", row.names = 1)
AGLExpr <- read.csv("/data3/users/Qiuqoo/CM_ZJU/article/GRN/AGL.csv", row.names = 1)
LBDExpr <- read.csv("/data3/users/Qiuqoo/CM_ZJU/article/GRN/LBD.csv", row.names = 1)
ltExpr <- read.csv("/data3/users/Qiuqoo/CM_ZJU/article/GRN/lt.csv", row.names = 1)
ERF1Expr <- read.csv("/data3/users/Qiuqoo/CM_ZJU/article/GRN/ERF1.csv", row.names = 1)
AUXINExpr <- read.csv("/data3/users/Qiuqoo/CM_ZJU/article/GRN/AUXIN.csv", row.names = 1)
CKExpr <- read.csv("/data3/users/Qiuqoo/CM_ZJU/article/GRN/CK.csv", row.names = 1)
EthyleneExpr <- read.csv("/data3/users/Qiuqoo/CM_ZJU/article/GRN/Ethylene.csv", row.names = 1)


result_interactions <- generate_interactions(fullExpr[,4:10], WOXExpr, threshold_r = 0.8, threshold_p = 0.01, topN = 300)
AGL_result_interactions <- generate_interactions(fullExpr[,4:10], AGLExpr, threshold_r = 0.8, threshold_p = 0.01, topN = 300)
LBD_result_interactions <- generate_interactions(fullExpr[,4:10], LBDExpr, threshold_r = 0.8, threshold_p = 0.01, topN = 300)
lt_result_interactions <- generate_interactions(fullExpr[,4:10], ltExpr, threshold_r = 0.8, threshold_p = 0.01, topN = 300)
ERF1_result_interactions <- generate_interactions(fullExpr[,4:10], ERF1Expr, threshold_r = 0.8, threshold_p = 0.01, topN = 300)
AUXIN_result_interactions <- generate_interactions(fullExpr[,4:10], AUXINExpr, threshold_r = 0.8, threshold_p = 0.01, topN = 300)
CK_result_interactions <- generate_interactions(fullExpr[,4:10], CKExpr, threshold_r = 0.8, threshold_p = 0.01, topN = 300)
Ethylene_result_interactions <- generate_interactions(fullExpr[,4:10], EthyleneExpr, threshold_r = 0.8, threshold_p = 0.01, topN = 300)

head(result_interactions)
dim(result_interactions)

write.csv(result_interactions, file = "/data3/users/Qiuqoo/CM_ZJU/article/GRN/WOX_top300_interactions.csv", row.names = FALSE)
write.csv(AGL_result_interactions, file = "/data3/users/Qiuqoo/CM_ZJU/article/GRN/AGL_top300_interactions.csv", row.names = FALSE)
write.csv(LBD_result_interactions, file = "/data3/users/Qiuqoo/CM_ZJU/article/GRN/LBD_top300_interactions.csv", row.names = FALSE)
write.csv(lt_result_interactions, file = "/data3/users/Qiuqoo/CM_ZJU/article/GRN/PLT2_top300_interactions.csv", row.names = FALSE)
write.csv(ERF1_result_interactions, file = "/data3/users/Qiuqoo/CM_ZJU/article/GRN/ERF1_top300_interactions.csv", row.names = FALSE)
write.csv(AUXIN_result_interactions, file = "/data3/users/Qiuqoo/CM_ZJU/article/GRN/AUXIN_top300_interactions.csv", row.names = FALSE)
write.csv(CK_result_interactions, file = "/data3/users/Qiuqoo/CM_ZJU/article/GRN/CK_top300_interactions.csv", row.names = FALSE)
write.csv(Ethylene_result_interactions, file = "/data3/users/Qiuqoo/CM_ZJU/article/GRN/Ethylene_result_interactions.csv", row.names = FALSE)

