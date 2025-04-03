library(ComplexUpset)
library(ggplot2)
library(dplyr)
library(tidyr)

file_path <- "/data3/users/Qiuqoo/CM_ZJU/DEG_DMG/DEG/nearly/"

samples <- c("NSEvsGE", "GEvsHE", "HEvsPTE", "PTEvsLTE", "LTEvsPCE", "PCEvsLCE")

DEG_lists <- list()
for (sample in samples) {
  file_name <- paste0(file_path, sample, "_deg_all-0824.csv")  #please confirm
  df <- read.csv(file_name)  
  DEG_lists[[sample]] <- df$gene_id  #please confirm
}

all_genes <- unique(unlist(DEG_lists))

DEG_matrix <- data.frame(Gene = all_genes)
for (sample in samples) {
  DEG_matrix[[sample]] <- ifelse(DEG_matrix$Gene %in% DEG_lists[[sample]], 1, 0)
}


DEG_matrix <- DEG_matrix[, comparison_order] #sort

upset_plot <- upset(
  DEG_matrix, 
  comparison_order,  
  name = "DEGs Intersections",
  width_ratio = 0.2, 
  height_ratio = 0.4,
  min_size = 10,  
  sort_sets = FALSE, 
  sort_intersections = FALSE 
)



