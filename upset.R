library(ComplexUpset)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)  # 用于遍历组合

# 设定差异基因文件路径
file_path <- "/data3/users/Qiuqoo/CM_ZJU/DEG_DMG/DEG/nearly/"

# 样品组名称
samples <- c("NSEvsGE", "GEvsHE", "HEvsPTE", "PTEvsLTE", "LTEvsPCE", "PCEvsLCE")

# 读取所有 DEGs 并存入列表
DEG_lists <- list()
for (sample in samples) {
  file_name <- paste0(file_path, sample, "_deg_all-0824.csv")
  df <- read.csv(file_name)  # 读取 CSV 文件
  DEG_lists[[sample]] <- df$gene_id  # 提取基因 ID 并存入列表
}

# 创建一个包含所有基因的全集
all_genes <- unique(unlist(DEG_lists))

# 构建二进制矩阵，表示基因在不同样品组中的分布
DEG_matrix <- data.frame(Gene = all_genes)
for (sample in samples) {
  DEG_matrix[[sample]] <- ifelse(DEG_matrix$Gene %in% DEG_lists[[sample]], 1, 0)
}

# 为保证后续操作不丢失 Gene 列，将列顺序设置为：Gene 列在第一位，后面是各个样品
comparison_order <- c("Gene", samples)
DEG_matrix <- DEG_matrix[, comparison_order]

### 生成 UpSet 图（原始代码段）
BiocManager::install("forcats")  # 如需安装该包，否则可省略

upset_plot <- upset(
  DEG_matrix, 
  comparison_order,  # 按指定顺序排列
  name = "DEGs Intersections",
  width_ratio = 0.2, 
  height_ratio = 1,
  min_size = 20,  
  sort_sets = FALSE,  # 不自动排序行
  sort_intersections = FALSE  # 不自动排序列
)

### 提取 UpSet 图中的各个交集（子集）对应的基因集合
subsets_list <- list()

library(purrr)

# 遍历所有可能的组合（1 到 length(samples) 的所有子集）
for (i in 1:length(samples)) {
  
  combn(samples, i, simplify = FALSE) %>%
    walk(function(combination) {
      # 定义组合名称，多个样本用 & 连接
      combo_name <- paste(combination, collapse = "&")
      
      # 在这些样本中，该基因均为 1，同时在其他样本中均为 0
      # 使用 Base R 进行逻辑筛选，注意 select 的范围需要加上 drop = FALSE
      included <- apply(DEG_matrix[, combination, drop = FALSE], 1, function(x) all(x == 1))
      excluded <- rowSums(DEG_matrix[, setdiff(samples, combination), drop = FALSE]) == 0
      
      selected_genes <- DEG_matrix$Gene[included & excluded]
      
      if (length(selected_genes) > 0) {
        subsets_list[[combo_name]] <<- selected_genes
      }
    })
}



write.csv(subsets_list[["NSEvsGE"]],"/data3/users/Qiuqoo/CM_ZJU/DEG_DMG/DEG/nearly/upset/NSEvsGE_upset_unique.csv")
write.csv(subsets_list[["GEvsHE"]],"/data3/users/Qiuqoo/CM_ZJU/DEG_DMG/DEG/nearly/upset/GEvsHE_upset_unique.csv")
write.csv(subsets_list[["HEvsPTE"]],"/data3/users/Qiuqoo/CM_ZJU/DEG_DMG/DEG/nearly/upset/HEvsPTE_upset_unique.csv")
write.csv(subsets_list[["PTEvsLTE"]],"/data3/users/Qiuqoo/CM_ZJU/DEG_DMG/DEG/nearly/upset/PTEvsLTE_upset_unique.csv")
write.csv(subsets_list[["LTEvsPCE"]],"/data3/users/Qiuqoo/CM_ZJU/DEG_DMG/DEG/nearly/upset/LTEvsPCE_upset_unique.csv")
write.csv(subsets_list[["PCEvsLCE"]],"/data3/users/Qiuqoo/CM_ZJU/DEG_DMG/DEG/nearly/upset/PCEvsLCE_upset_unique.csv")


# 若想查看各个组合中基因的数量
subset_counts <- sapply(subsets_list, length)
print(subset_counts)



