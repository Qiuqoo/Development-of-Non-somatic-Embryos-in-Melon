library(tidyverse)
library(igraph)
library(ggraph)


df <- read_excel("/data3/users/Qiuqoo/CM_ZJU/cluster_regulation/cluster_regulation.xls", sheet = 1)%>%
  mutate(
    TFcluster = as.character(TFcluster),
    targetcluster = as.character(targetcluster)
  )


cluster_size <- c(
  C1 = 1360,
  C2 = 2037,
  C3 = 2467,
  C4 = 1709,
  C5 = 2722,
  C6 = 5153,
  C7 = 3417,
  C8 = 2944,
  C9 = 2534
)


df2 <- df %>% rowwise() %>% mutate({

  a <- target_numer

  b <- total_targert_numer - a

  c <- cluster_size[paste0("C", targetcluster)] - a

  total_DEGs <- sum(cluster_size)
  d <- (total_DEGs - cluster_size[paste0("C", targetcluster)]) - b
  

  ft <- fisher.test(matrix(c(a, b, c, d), nrow=2),
                    alternative = "greater")
  tibble(pvalue = ft$p.value)
}) %>% ungroup()


sig_edges <- df2 %>% filter(pvalue < 0.05)

sig_edges2 <- sig_edges %>%
  mutate(
    from = paste0("C", TFcluster),
    to   = paste0("C", targetcluster)
  )


g <- graph_from_data_frame(
  d = sig_edges2 %>% dplyr::select(from, to, weight = target_numer),
  vertices = tibble(name = paste0("C", 1:9)),
  directed = TRUE
)



cluster_colors <- c(
  "C1" = "#D62728","C2" = "#FFEB59","C3" = "#6AC4A6",
  "C4" = "#1F77B4","C5" = "#FF5733","C6" = "#FDB462",
  "C7" = "#FFA07A","C8" = "#B0E57C","C9" = "#9370DB"
)

p <- ggraph(g, layout="circle") +
  geom_edge_arc(aes(width = weight),
                arrow      = arrow(length = unit(3, 'mm')),
                end_cap    = circle(2, 'mm'),
                edge_colour= "grey50",
                alpha      = 0.7) +
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_point(aes(fill = name), shape = 21, size = 8) +
  geom_node_text(aes(label = name), color = "white", size = 4) +
  scale_fill_manual(values = cluster_colors) +
  theme_void() +
  labs(title    = "Significant TF Cluster Regulatory Network",
       subtitle = "Filtered by one-sided Fisher’s exact test (P < 0.05)") +
  theme(
    plot.title    = element_text(hjust=0.5, size=16),
    plot.subtitle = element_text(hjust=0.5, size=12)
  )

print(p)

