
if (!require(plotly)) install.packages("plotly")
library(plotly)


path <- "/data3/users/Qiuqoo/CM_ZJU/epiregeneration/"
epiDEG           <- scan(file.path(path, "epiDEG.txt"),           what = "", sep="\n")
epiDMG           <- scan(file.path(path, "epiDMG.txt"),           what = "", sep="\n")
gene_with_siRNA  <- scan(file.path(path, "gene_with_siRNA.txt"),  what = "", sep="\n")
epiAccessibility <- scan(file.path(path, "epiAccessibility.txt"), what = "", sep="\n")


all_genes <- sort(unique(c(epiDEG, epiDMG, gene_with_siRNA, epiAccessibility)))
df <- data.frame(
  gene   = all_genes,
  DEG    = all_genes %in% epiDEG,
  DMG    = all_genes %in% epiDMG,
  RNA    = all_genes %in% gene_with_siRNA,
  Acc    = all_genes %in% epiAccessibility,
  stringsAsFactors = FALSE
)


df$path <- with(df, 
                paste0(
                  ifelse(DEG,    "epiDEG",    "noDEG"),    "/", 
                  ifelse(DMG,    "epiDMG",    "noDMG"),    "/", 
                  ifelse(RNA,    "with_siRNA","noRNA"),   "/", 
                  ifelse(Acc,    "accessible","noAcc")
                )
)


library(dplyr)
sb <- df %>% 
  count(path) %>% 

  mutate(
    labels = gsub("/", "-", path),
    values = n
  )


fig <- plot_ly(
  labels = sb$labels,
  parents = sub("^[^/]+/(.*)$", "\\1", sub("^[^/]+/(.*)$", "\\1", sub("^[^/]+/(.*)$", "\\1", sb$labels ))),
  values = sb$values,
  type = 'sunburst',
  branchvalues = 'total',
  maxdepth = 4,
  textinfo = 'label+value'
)


if (!require(plotly)) install.packages("plotly")
library(plotly)


labels <- sb$labels      
values <- sb$values      


parents <- sapply(labels, function(lbl) {
  parts <- strsplit(lbl, "-", fixed = TRUE)[[1]]
  if (length(parts) == 1) {
    return(NA)           
  } else {
    return(paste(parts[-length(parts)], collapse = "-"))
  }
})



if (!require(sunburstR)) install.packages("sunburstR")
library(sunburstR)


df_sb <- data.frame(
  path  = sb$path,
  count = sb$n,
  stringsAsFactors = FALSE
)


sunburst(
  df_sb,
  count = TRUE,           
  percent = FALSE,         
  legend = TRUE,           
  withD3 = TRUE          
)
