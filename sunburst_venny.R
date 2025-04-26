## R script to draw a sunburst chart showing presence/absence intersections
# of four omics layers: epiDEG, epiDMG, siRNA, accessibility

# 1. Install and load necessary packages (run once)
# install.packages(c("dplyr", "tidyr", "sunburstR"))
library(dplyr)
library(tidyr)
library(sunburstR)
# 2. Read in your summary data (replace with your actual file or data frame)
# sb <- read.table("sb_summary.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# 安装并加载必要包
if (!require(plotly)) install.packages("plotly")
library(plotly)
if (!require(pagedown)) install.packages("pagedown")
library(pagedown)
# 1. 读入四个基因列表（无表头，单列ID）
path <- "/data3/users/Qiuqoo/CM_ZJU/epiregeneration/"
epiDEG           <- scan(file.path(path, "epiDEG.txt"),           what = "", sep="\n")
epiDMG           <- scan(file.path(path, "epiDMG.txt"),           what = "", sep="\n")
#gene_with_siRNA  <- scan(file.path(path, "gene_with_siRNA.txt"),  what = "", sep="\n")
epiAccessibility <- scan(file.path(path, "epiAccessibility.txt"), what = "", sep="\n")

# 2. 构建全集并打标
all_genes <- sort(unique(c(epiDEG, epiDMG, epiAccessibility)))
df <- data.frame(
  gene   = all_genes,
  DEG    = all_genes %in% epiDEG,
  DMG    = all_genes %in% epiDMG,
  Acc    = all_genes %in% epiAccessibility,
  stringsAsFactors = FALSE
)

# 3. 为每个基因生成“层级路径”：每层是该Omics的 yes/no
df$path <- with(df, 
                paste0(
                  ifelse(DEG,    "epiDEG",    "noDEG"),    "/", 
                  ifelse(DMG,    "epiDMG",    "noDMG"),    "/", 
                  ifelse(Acc,    "accessible","noAcc")
                )
)

# 4. 统计每个路径的基因数
library(dplyr)
sb <- df %>% 
  count(path) %>% 
  # plotly需要把分隔符换成“-”或空格，且不能有R的逻辑值
  mutate(
    labels = gsub("/", "-", path),
    values = n
  )

# Example: manually defining sb for demonstration
# sb <- data.frame(
#   path = c(
#     "epiDEG/epiDMG/noRNA/accessible",    "epiDEG/epiDMG/noRNA/noAcc",
#     "epiDEG/epiDMG/with_siRNA/accessible","epiDEG/epiDMG/with_siRNA/noAcc",
#     "epiDEG/noDMG/noRNA/accessible",     "epiDEG/noDMG/noRNA/noAcc",
#     "epiDEG/noDMG/with_siRNA/accessible","epiDEG/noDMG/with_siRNA/noAcc",
#     "noDEG/epiDMG/noRNA/accessible",     "noDEG/epiDMG/noRNA/noAcc",
#     "noDEG/epiDMG/with_siRNA/accessible","noDEG/epiDMG/with_siRNA/noAcc",
#     "noDEG/noDMG/noRNA/accessible",      "noDEG/noDMG/with_siRNA/accessible",
#     "noDEG/noDMG/with_siRNA/noAcc"
#   ),
#   n = c(341, 1113, 219, 401, 565, 1781, 382, 717, 885, 4149, 587, 992, 1493, 693, 1713),
#   stringsAsFactors = FALSE
# )

# 3. Split the path into four layers and label presence/absence
sb_long <- sb %>%
  separate(path, into = c("DEG", "DMG", "Accessibility"), sep = "/") %>%
  mutate(
    DEG = ifelse(DEG == "epiDEG", "Present", "Absent"),
    DMG = ifelse(DMG == "epiDMG", "Present", "Absent"),
    Accessibility = ifelse(Accessibility == "accessible", "Present", "Absent")
  ) %>%
  # 4. Create sequence column for sunburstR
  mutate(sequence = paste(DEG, DMG,  Accessibility, sep = "-"))

# 5. Prepare data frame for sunburstR (sequence, count)
sb_sunburst <- sb_long %>% select(sequence, count = n)

# 6. Draw the sunburst chart
# This will generate an interactive HTML widget
sunburst_widget <- sunburst(
  sb_sunburst,
  count = TRUE,             # 显示计数
  percent = TRUE,           # 显示整体占比
  # JS 函数：在图中心显示名称、计数及百分比
  # explanation = JS(
  #   "function(d){return d.data.name + ': ' + d.value + ' (' + (d.percent*100).toFixed(1) + '%)';}"
  # ),
  # 自定义颜色映射，按 Present/Absent 着色
  # colors = JS(
  #   "d3.scaleOrdinal().domain(['Present','Absent']).range(['#1f78b4','#33a02c'])"
  # ),
  legend = TRUE           # 显示状态图例（Present/Absent）
)

# 7. 添加各层图注（Ring Legend）
library(htmltools)
legend_html <- tags$div(
  style = "font-size:14px; margin-bottom:10px;", 
  tags$strong("图注："), 
  tags$ul(
    tags$li("第一环(中心): epiDEG"),
    tags$li("第二环: epiDMG"),
    tags$li("第三环: siRNA"),
    tags$li("第四环: Accessibility")
  )
)

# 8. 在浏览器中展示

