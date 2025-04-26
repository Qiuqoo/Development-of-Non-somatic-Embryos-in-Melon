library(readr)
library(dplyr)
library(ggplot2)

setwd("/data3/users/Qiuqoo/CM_ZJU/GO")


# 设置文件名
file <- "ERF1_target_gene.xls"

# 读取文件
df <- read_tsv(file)

# 更改列名
colnames(df)[3:4] <- c("out", "all")

# 计算 Foldenrichment，筛选并排序
df <- df %>%
  mutate(Foldenrichment = out / all) %>%
  filter(out > 5) %>%
  arrange(desc(Foldenrichment))

# 选前5个结果
df_top <- head(df, 20)

# 绘图，固定 x 轴范围
p <- ggplot(df_top, aes(x = Foldenrichment, y = reorder(Description, Foldenrichment), fill = Foldenrichment)) +
  geom_col() +
  scale_fill_gradient(low = "#008744", high = "#673ab7") +  # 黄色到绿色渐变
  labs(title = file, x = "Fold enrichment", y = "GO term") +
  theme_minimal(base_size = 12)


# 显示图像
print(p)

# 保存 PDF
ggsave(paste0(tools::file_path_sans_ext(file), "_fold_plot.pdf"), plot = p, width = 8, height = 6)
