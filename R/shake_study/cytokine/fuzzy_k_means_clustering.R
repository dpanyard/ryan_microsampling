##
no_function()

sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())

###load data
load("data/shake_study/cytokine_data_analysis/data_preparation/expression_data")
load("data/shake_study/cytokine_data_analysis/data_preparation/sample_info")
load("data/shake_study/cytokine_data_analysis/data_preparation/variable_info")

sxtTools::setwd_project()
setwd("data/shake_study/cytokine_data_analysis/k_means_clustering")

dim(expression_data)
dim(sample_info)
dim(variable_info)

library(Mfuzz)
library(e1071)

sample_info$sample_id == colnames(expression_data)
sample_info$TP


temp_data <- 
  log(expression_data + 1, 2)

sample_info$TP

tp <- c(0, 30, 60, 120, 240)

temp_data_mean <- 
purrr::map(tp, function(x) {
  temp_idx <- which(sample_info$TP == x)
  temp_data[,temp_idx] %>% 
    apply(1, mean)
}) %>% 
  do.call(cbind, .) %>% 
  as.data.frame()

colnames(temp_data_mean) <- paste("Time", tp, sep = "_")

library(Mfuzz)
#first get the time point data together:
# bind that to the data frame
##scale
temp_data <- 
  temp_data_mean %>% 
  apply(1, function(x) (x - mean(x))/sd(x)) %>% 
  t() %>% 
  as.data.frame()

time <- c(0, 30, 60, 120, 240)

temp_data <- rbind(time, temp_data)

row.names(temp_data)[1] <- "time"

write.table(
  temp_data,
  file = "temp_data.txt",
  sep = '\t',
  quote = FALSE,
  col.names = NA
)

#read it back in as an expression set
data <- table2eset(filename = "temp_data.txt")

data.s <- standardise(data)
m1 <- mestimate(data.s)
m1

plot <-
Dmin(data.s, m=m1, crange=seq(2,22,1), repeats=3, visu=TRUE)

plot <-
plot %>%
  data.frame(distance = plot,
             k = seq(2,22,1)) %>%
  ggplot(aes(k, distance)) +
  geom_point(shape = 21, size = 4, fill = "black") +
  geom_smooth() +
  geom_segment(aes(x = k, y = 0, xend = k, yend = distance)) +
  theme_bw() +
  theme(
    # legend.position = c(0, 1),
    # legend.justification = c(0, 1),
    panel.grid = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  labs(
    x = "Cluster number",
    y = "Min. centroid distance"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))

plot

ggsave(plot, filename = "distance_k_number.pdf", width = 7, height = 7)

cluster = 4

c <- mfuzz(data.s, c = cluster, m = m1)

mfuzz.plot(eset = data.s,
           # min.mem = 0.6,
           cl = c,
           mfrow=c(4,4),
           time.labels = time,
           new.window = FALSE)

save(c, file = "c")
load("c")


####any two clusters with correlation > 0.8 should be considered as one
library(corrplot)
layout(1)
center <- c$centers
rownames(center) <- paste("Cluster", rownames(center), sep = ' ')

corrplot::corrplot(
  corr = cor(t(center)),
  type = "full",
  diag = TRUE,
  order = "hclust", 
  hclust.method = "ward.D", 
  # addrect = 5, 
  col = colorRampPalette(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(n=100),
  number.cex = .7, 
  addCoef.col = "black"
)

###cor_plot_cluster

center %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "cluster") %>% 
  tidyr::pivot_longer(cols = -cluster, names_to = "point", values_to = "value") %>% 
  # dplyr::filter(cluster %in% c("Cluster 2", "Cluster 4")) %>% 
  dplyr::mutate(point = factor(point, levels = unique(point))) %>% 
  ggplot(aes(point, value, group = cluster)) +
  geom_point(aes(color = cluster)) +
  geom_line(aes(group = cluster, color = cluster)) 
# geom_smooth(aes(group = cluster), se = FALSE)

library(ComplexHeatmap)

plot <- 
  Heatmap(
    matrix = center,
    cluster_columns = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D",
    name = "Z-score",
    border = TRUE,
    rect_gp = gpar(col= "white")
  )


plot <- ggplotify::as.ggplot(plot)
plot
ggsave(plot, filename = "cluster_heatmap.pdf", width = 9, height = 7)

###
cluster_color <- 
  colorRampPalette(colors = RColorBrewer::brewer.pal(n = 11, name = "Spectral"))(8)

names(cluster_color) <- as.character(1:8)

plot <- 
  center %>%
  as.data.frame() %>%
  tibble::rowid_to_column(var = "cluster") %>%
  tidyr::pivot_longer(cols = -cluster,
                      names_to = "time",
                      values_to = "value") %>%
  dplyr::mutate(time = factor(time, levels = unique(time))) %>%
  dplyr::mutate(cluster = as.character(cluster)) %>%
  ggplot(aes(time, value)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(group = cluster, color = cluster), show.legend = FALSE) +
  # geom_point(aes(group = cluster, fill = cluster), show.legend = FALSE, shape = 21, size = 3) +
  # geom_smooth(aes(color = cluster, group = cluster), se = FALSE) +
  scale_color_manual(values = cluster_color) +
  scale_fill_manual(values = cluster_color) +
  # facet_grid(rows = vars(class)) +
  theme_bw() +
  labs(x = "", y = "Z-score") +
  theme(
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 12
    ),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

ggsave(plot, filename = "cluster_center_plot.pdf", width = 7, height = 4)

centers <- c$centers

cluster_info <-
  data.frame(
    variable_id = names(c$cluster),
    c$membership,
    cluster = c$cluster,
    stringsAsFactors = FALSE
  ) %>%
  arrange(cluster)

xlsx::write.xlsx(cluster_info,
                 "cluster_info.xlsx",
                 row.names = FALSE)


####plot for each cluster
idx <- 3
for(idx in 1:4) {
  cat(idx, " ")
  cluster_data <-
    cluster_info %>%
    dplyr::filter(cluster == idx) %>%
    dplyr::select(1, 1 + idx, cluster)
  
  colnames(cluster_data)[2] <- c("membership")
  
  cluster_data <-
    cluster_data %>%
    dplyr::filter(membership > 0.2)
  path <- paste("cluster", idx, sep = "_")
  dir.create(path)
  
  xlsx::write.xlsx(cluster_data,
                   file = file.path(path, paste("cluster", idx, ".xlsx", sep = "")),
                   row.names = FALSE)
  
  temp_center <-
    centers[idx, , drop = TRUE] %>%
    data.frame(time = names(.),
               value = .,
               stringsAsFactors = FALSE) %>%
    dplyr::mutate(time = factor(time, levels = time))
  
  plot <-
    temp_data[cluster_data$variable_id, ] %>%
    data.frame(
      membership = cluster_data$membership,
      .,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ) %>%
    tibble::rownames_to_column(var = "variable_id") %>%
    tidyr::pivot_longer(
      cols = -c(variable_id, membership),
      names_to = "time",
      values_to = "value"
    ) %>%
    dplyr::mutate(time = factor(time, levels = unique(time))) %>%
    ggplot(aes(time, value, group = variable_id)) +
    geom_line(aes(color = membership)) +
    scale_color_gradientn(colours = c(RColorBrewer::brewer.pal(n = 9, name = "OrRd")[c(1:7)])) +
    theme_bw() +
    theme(
      legend.position = c(0, 1),
      legend.justification = c(0, 1),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 12
      ),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA)
    ) +
    labs(
      x = "",
      y = "Z-score",
      title = paste(
        "Cluster ",
        idx,
        " (",
        nrow(cluster_data),
        " metabolic peaks)",
        sep = ""
      )
    ) +
    geom_line(
      mapping = aes(time, value, group = 1),
      data = temp_center,
      size = 2
    )
  
  plot

  ggsave(plot, filename = file.path(path, paste("cluster",idx, ".pdf", sep = "")),
         width = 8, height = 7)

}


# dim(cluster_data)

###annotation for each cluster
cluster1 <- readxl::read_xlsx("cluster_1/cluster1.xlsx") 
cluster2 <- readxl::read_xlsx("cluster_2/cluster2.xlsx") 
cluster3 <- readxl::read_xlsx("cluster_3/cluster3.xlsx") 
cluster4 <- readxl::read_xlsx("cluster_4/cluster4.xlsx") 
cluster5 <- readxl::read_xlsx("cluster_5/cluster5.xlsx") 
cluster6 <- readxl::read_xlsx("cluster_6/cluster6.xlsx") 
cluster7 <- readxl::read_xlsx("cluster_7/cluster7.xlsx") 
cluster8 <- readxl::read_xlsx("cluster_8/cluster8.xlsx") 

nrow(cluster1) +
  nrow(cluster2) +
  nrow(cluster3) +
  nrow(cluster4) +
  nrow(cluster5) +
  nrow(cluster6) +
  nrow(cluster7) +
  nrow(cluster8) 

##metabolites that not in clusters
non_metabolite <- 
  setdiff(cluster_info$variable_id,
          c(cluster4$variable_id,
            cluster2$variable_id,
            cluster3$variable_id,
            cluster4$variable_id,
            cluster5$variable_id,
            cluster6$variable_id,
            cluster7$variable_id))

plot <-
  temp_data[non_metabolite,] %>%
  tibble::rownames_to_column(var = "variable_id") %>%
  tidyr::pivot_longer(cols = -c(variable_id),
                      names_to = "time", 
                      values_to = "value") %>%
  dplyr::mutate(time = factor(time, levels = unique(time))) %>%
  ggplot(aes(time, value, group = variable_id)) +
  geom_line(alpha = 1) +
  # scale_color_gradientn(colours = c(
  #   RColorBrewer::brewer.pal(n = 9, name = "OrRd")[c(1:7)]
  # )) +
  theme_bw() +
  theme(
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 12
    ),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  labs(
    x = "",
    y = "Z-score",
    title = paste("Unclustered metabolic peaks (", length(non_metabolite), " metabolic peaks)", sep = "")
  ) 
  # geom_rect(
  #   aes(
  #     xmin = 14.5,
  #     xmax = Inf,
  #     ymin = -Inf,
  #     ymax = Inf
  #   ),
  #   fill = ggsci::pal_aaas()(n = 10)[5],
  #   alpha = 0.5,
  #   data = data.frame(),
  #   inherit.aes = FALSE
  # )

plot

ggsave(plot, filename = "non_cluster_metabolites.pdf",
       width = 8, height = 7)

dim(cluster1)
dim(cluster2)
dim(cluster3)
dim(cluster4)
dim(cluster5)
dim(cluster6)
dim(cluster7)
dim(cluster8)


###heatmap for all the samples
library(ComplexHeatmap)
temp_data <-
  temp_data_mean[c(
    cluster4$variable_id,
    cluster8$variable_id,
    cluster3$variable_id,
    cluster6$variable_id,
    cluster7$variable_id,
    cluster1$variable_id,
    cluster2$variable_id,
    cluster5$variable_id
  ), ]

temp_data <- 
  apply(temp_data, 1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

range(temp_data)
library(circlize)
col_fun = colorRamp2(c(-3, 0, 3), c("#4292C6", "white", "red"))

cluster <-
  c(
    rep(4, nrow(cluster4)),
    rep(8, nrow(cluster8)),
    rep(3, nrow(cluster3)),
    rep(6, nrow(cluster6)),
    rep(7, nrow(cluster7)),
    rep(1, nrow(cluster1)),
    rep(2, nrow(cluster2)),
    rep(5, nrow(cluster5))
  )

cluster_color <- 
  colorRampPalette(colors = RColorBrewer::brewer.pal(n = 11, name = "Spectral"))(8)

names(cluster_color) <- as.character(1:8)

ha <-
  rowAnnotation(
    Cluster = cluster,
    col = list(
      Cluster = cluster_color
    )
    # annotation_name_side = c("left")
  )

plot <- 
  Heatmap(
    temp_data,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_names = FALSE,
    border = TRUE,
    col = col_fun,
    right_annotation = ha,
    name = "Z score",
    column_names_rot = 0
    # row_split = factor(as.character(cluster), levels = c(6,3,8,1,4,5,2,9,7))
  )

plot <- ggplotify::as.ggplot(plot)

plot

ggsave(plot = plot, filename = "heatmap_for_all_cluster.pdf", width = 7, height = 7)



###functional annotation for different cluster
cluster1 %>% 
  dplyr::filter(!is.na(Level)) %>% 
  dplyr::filter(Level == 1 | Level == 2) %>% 
  dplyr::pull(Compound.name)

cluster2 %>% 
  dplyr::filter(!is.na(Level)) %>% 
  dplyr::filter(Level == 1 | Level == 2) %>% 
  dplyr::pull(Compound.name)

cluster3 %>% 
  dplyr::filter(!is.na(Level)) %>% 
  dplyr::filter(Level == 1 | Level == 2) %>% 
  dplyr::pull(Compound.name)

cluster4 %>% 
  dplyr::filter(!is.na(Level)) %>% 
  dplyr::filter(Level == 1 | Level == 2) %>% 
  dplyr::pull(Compound.name)

cluster5 %>% 
  dplyr::filter(!is.na(Level)) %>% 
  dplyr::filter(Level == 1 | Level == 2) %>% 
  dplyr::pull(Compound.name)

cluster6 %>% 
  dplyr::filter(!is.na(Level)) %>% 
  dplyr::filter(Level == 1 | Level == 2) %>% 
  dplyr::pull(Compound.name)

cluster7 %>% 
  dplyr::filter(!is.na(Level)) %>% 
  dplyr::filter(Level == 1 | Level == 2) %>% 
  dplyr::pull(Compound.name)


####output file for PIUMet
dir.create("PIUMet")

cluster1_piumet <-
  cluster1 %>%
  dplyr::select(variable_id, mz, fdr) %>%
  dplyr::mutate(fdr = -log(fdr, 10)) %>%
  dplyr::mutate(polarity = case_when(
    stringr::str_detect(variable_id, "POS") ~ "positive",
    stringr::str_detect(variable_id, "NEG") ~ "negative"
  )) %>%
  dplyr::select(mz, polarity, fdr) %>%
  dplyr::arrange(desc(fdr))

colnames(cluster1_piumet) <- NULL

# write.table(
#   x = cluster1_piumet,
#   file = file.path("PIUMet","cluster1_piumet.txt"),
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t"
# )

cluster2_piumet <-
  cluster2 %>%
  dplyr::select(variable_id, mz, fdr) %>%
  dplyr::mutate(fdr = -log(fdr, 10)) %>%
  dplyr::mutate(polarity = case_when(
    stringr::str_detect(variable_id, "POS") ~ "positive",
    stringr::str_detect(variable_id, "NEG") ~ "negative"
  )) %>%
  dplyr::select(mz, polarity, fdr) %>%
  dplyr::arrange(desc(fdr))

colnames(cluster2_piumet) <- NULL

# write.table(
#   x = cluster2_piumet,
#   file = file.path("PIUMet","cluster2_piumet.txt"),
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t"
# )





cluster3_piumet <-
  cluster3 %>%
  dplyr::select(variable_id, mz, fdr) %>%
  dplyr::mutate(fdr = -log(fdr, 10)) %>%
  dplyr::mutate(polarity = case_when(
    stringr::str_detect(variable_id, "POS") ~ "positive",
    stringr::str_detect(variable_id, "NEG") ~ "negative"
  )) %>%
  dplyr::select(mz, polarity, fdr) %>%
  dplyr::arrange(desc(fdr))

colnames(cluster3_piumet) <- NULL

# write.table(
#   x = cluster3_piumet,
#   file = file.path("PIUMet","cluster3_piumet.txt"),
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t"
# )



cluster4_piumet <-
  cluster4 %>%
  dplyr::select(variable_id, mz, fdr) %>%
  dplyr::mutate(fdr = -log(fdr, 10)) %>%
  dplyr::mutate(polarity = case_when(
    stringr::str_detect(variable_id, "POS") ~ "positive",
    stringr::str_detect(variable_id, "NEG") ~ "negative"
  )) %>%
  dplyr::select(mz, polarity, fdr) %>%
  dplyr::arrange(desc(fdr))

colnames(cluster4_piumet) <- NULL

# write.table(
#   x = cluster4_piumet,
#   file = file.path("PIUMet","cluster4_piumet.txt"),
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t"
# )


cluster5_piumet <-
  cluster5 %>%
  dplyr::select(variable_id, mz, fdr) %>%
  dplyr::mutate(fdr = -log(fdr, 10)) %>%
  dplyr::mutate(polarity = case_when(
    stringr::str_detect(variable_id, "POS") ~ "positive",
    stringr::str_detect(variable_id, "NEG") ~ "negative"
  )) %>%
  dplyr::select(mz, polarity, fdr) %>%
  dplyr::arrange(desc(fdr))

colnames(cluster5_piumet) <- NULL

# write.table(
#   x = cluster5_piumet,
#   file = file.path("PIUMet","cluster5_piumet.txt"),
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t"
# )




cluster6_piumet <-
  cluster6 %>%
  dplyr::select(variable_id, mz, fdr) %>%
  dplyr::mutate(fdr = -log(fdr, 10)) %>%
  dplyr::mutate(polarity = case_when(
    stringr::str_detect(variable_id, "POS") ~ "positive",
    stringr::str_detect(variable_id, "NEG") ~ "negative"
  )) %>%
  dplyr::select(mz, polarity, fdr) %>%
  dplyr::arrange(desc(fdr))

colnames(cluster6_piumet) <- NULL

# write.table(
#   x = cluster6_piumet,
#   file = file.path("PIUMet","cluster6_piumet.txt"),
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t"
# )



cluster7_piumet <-
  cluster7 %>%
  dplyr::select(variable_id, mz, fdr) %>%
  dplyr::mutate(fdr = -log(fdr, 10)) %>%
  dplyr::mutate(polarity = case_when(
    stringr::str_detect(variable_id, "POS") ~ "positive",
    stringr::str_detect(variable_id, "NEG") ~ "negative"
  )) %>%
  dplyr::select(mz, polarity, fdr) %>%
  dplyr::arrange(desc(fdr))

colnames(cluster7_piumet) <- NULL

# write.table(
#   x = cluster7_piumet,
#   file = file.path("PIUMet","cluster7_piumet.txt"),
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t"
# )

cluster1_info <- read.table("PIUMet/cluster1_piumet.txt", sep = "\t")
cluster2_info <- read.table("PIUMet/cluster2_piumet.txt", sep = "\t")
cluster3_info <- read.table("PIUMet/cluster3_piumet.txt", sep = "\t")
cluster4_info <- read.table("PIUMet/cluster4_piumet.txt", sep = "\t")
cluster5_info <- read.table("PIUMet/cluster5_piumet.txt", sep = "\t")
cluster6_info <- read.table("PIUMet/cluster6_piumet.txt", sep = "\t")
cluster7_info <- read.table("PIUMet/cluster7_piumet.txt", sep = "\t")

colnames(cluster1_info) <- 
  colnames(cluster2_info) <- 
  colnames(cluster3_info) <- 
  colnames(cluster4_info) <- 
  colnames(cluster5_info) <- 
  colnames(cluster6_info) <- 
  colnames(cluster7_info) <- 
  c("mz", "polarity", "fdr")

# ######cluster1
# readPIUMet(
#   path = "PIUMet/piumet_output_cluster1/",
#   marker = cluster1[, c("name", "mz", "rt")],
#   text = FALSE,
#   layout = "kk",
#   size_range = c(2, 8)
# )
# 
# 
# ######cluster2
# readPIUMet(
#   path = "PIUMet/piumet_output_cluster2/",
#   marker = cluster2[, c("name", "mz", "rt")],
#   text = FALSE,
#   layout = "kk",
#   size_range = c(2, 8)
# )
# 
# # ######cluster3
# # readPIUMet(
# #   path = "PIUMet/piumet_output_cluster3/",
# #   marker = cluster3[, c("name", "mz", "rt")],
# #   text = FALSE,
# #   layout = "kk",
# #   size_range = c(2, 8)
# # )
# 
# ######cluster4
# readPIUMet(
#   path = "PIUMet/piumet_output_cluster4/",
#   marker = cluster4[, c("name", "mz", "rt")],
#   text = FALSE,
#   layout = "kk",
#   size_range = c(2, 8)
# )
# 
# 
# ######cluster5
# readPIUMet(
#   path = "PIUMet/piumet_output_cluster5/",
#   marker = cluster5[, c("name", "mz", "rt")],
#   text = FALSE,
#   layout = "kk",
#   size_range = c(2, 8)
# )
# 
# ######cluster6
# readPIUMet(
#   path = "PIUMet/piumet_output_cluster6/",
#   marker = cluster6[, c("name", "mz", "rt")],
#   text = FALSE,
#   layout = "kk",
#   size_range = c(2, 8)
# )
# 
# 
# ######cluster7
# readPIUMet(
#   path = "PIUMet/piumet_output_cluster7/",
#   marker = cluster7[, c("name", "mz", "rt")],
#   text = FALSE,
#   layout = "kk",
#   size_range = c(2, 8)
# )


# ###cluster1
# load("PIUMet/piumet_output_cluster1/Result/node_data")
# 
# hmdb_id <-
#   node_data %>%
#   dplyr::filter(!is.null(HMDB_ID)) %>%
#   dplyr::filter(HMDB_ID != " " & HMDB_ID != "NULL") %>%
#   dplyr::pull(HMDB_ID) %>%
#   unique()
# 
# kegg_id <-
#   lapply(hmdb_id, function(x){
#     metflow2::transID(query = x,
#                       from = "Human cytokine Database",
#                       to = "KEGG", top = 1)
#   }) %>%
#   do.call(rbind, .)
# 
# kegg_id
# 
# load("PIUMet/hmdbMS1Database0.0.1")
# 
# hmdb_data <- hmdbMS1Database0.0.1@spectra.info
# 
# kegg_id2 <-
#   hmdb_data$KEGG.ID[match(kegg_id$`Human cytokine Database`, hmdb_data$HMDB.ID)]
# 
# kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)
# 
# colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")
# 
# KEGG_ID <-
#   apply(kegg_id, 1, function(x){
#     x <- as.character(x)
#     if(is.na(x[2]) & is.na(x[3])){
#       return(NA)
#     }
# 
#     if(!is.na(x[2])){
#       return(x[2])
#     }
# 
#     if(!is.na(x[3])){
#       return(x[3])
#     }
# 
# 
#   })
# 
# kegg_id <-
#   kegg_id %>%
#   dplyr::mutate(KEGG_ID = KEGG_ID) %>%
#   dplyr::select(-c(KEGG_ID1, KEGG_ID2))
# 
# kegg_id
# 
# node_data <-
#   node_data %>%
#   dplyr::left_join(kegg_id, by = c("HMDB_ID"))
# 
# annotation_result <- node_data
# 
# save(annotation_result, file = "PIUMet/piumet_output_cluster1/Result/annotation_result")
# 
# 
# ###cluster2
# load("PIUMet/piumet_output_cluster2/Result/node_data")
# 
# hmdb_id <-
#   node_data %>%
#   dplyr::filter(!is.null(HMDB_ID)) %>%
#   dplyr::filter(HMDB_ID != " " & HMDB_ID != "NULL") %>%
#   dplyr::pull(HMDB_ID) %>%
#   unique()
# 
# kegg_id <-
#   lapply(hmdb_id, function(x){
#     metflow2::transID(query = x,
#                       from = "Human cytokine Database",
#                       to = "KEGG", top = 1)
#   }) %>%
#   do.call(rbind, .)
# 
# kegg_id
# 
# # load("PIUMet/hmdbMS1Database0.0.1")
# 
# hmdb_data <- hmdbMS1Database0.0.1@spectra.info
# 
# kegg_id2 <-
#   hmdb_data$KEGG.ID[match(kegg_id$`Human cytokine Database`, hmdb_data$HMDB.ID)]
# 
# kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)
# 
# colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")
# 
# KEGG_ID <-
#   apply(kegg_id, 1, function(x){
#     x <- as.character(x)
#     if(is.na(x[2]) & is.na(x[3])){
#       return(NA)
#     }
# 
#     if(!is.na(x[2])){
#       return(x[2])
#     }
# 
#     if(!is.na(x[3])){
#       return(x[3])
#     }
# 
# 
#   })
# 
# kegg_id <-
#   kegg_id %>%
#   dplyr::mutate(KEGG_ID = KEGG_ID) %>%
#   dplyr::select(-c(KEGG_ID1, KEGG_ID2))
# 
# kegg_id
# 
# node_data <-
#   node_data %>%
#   dplyr::left_join(kegg_id, by = c("HMDB_ID"))
# 
# annotation_result <- node_data
# 
# save(annotation_result, file = "PIUMet/piumet_output_cluster2/Result/annotation_result")
# 
# 
# ###cluster3
# load("PIUMet/piumet_output_cluster3/Result/node_data")
# hmdb_id <-
#   node_data %>%
#   dplyr::filter(!is.null(HMDB_ID)) %>%
#   dplyr::filter(HMDB_ID != " " & HMDB_ID != "NULL") %>%
#   dplyr::pull(HMDB_ID) %>%
#   unique()
# 
# kegg_id <-
#   lapply(hmdb_id, function(x){
#     metflow2::transID(query = x,
#                       from = "Human cytokine Database",
#                       to = "KEGG", top = 1)
#   }) %>%
#   do.call(rbind, .)
# 
# kegg_id
# 
# # load("PIUMet/hmdbMS1Database0.0.1")
# 
# hmdb_data <- hmdbMS1Database0.0.1@spectra.info
# 
# kegg_id2 <-
#   hmdb_data$KEGG.ID[match(kegg_id$`Human cytokine Database`, hmdb_data$HMDB.ID)]
# 
# kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)
# 
# colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")
# 
# KEGG_ID <-
#   apply(kegg_id, 1, function(x){
#     x <- as.character(x)
#     if(is.na(x[2]) & is.na(x[3])){
#       return(NA)
#     }
# 
#     if(!is.na(x[2])){
#       return(x[2])
#     }
# 
#     if(!is.na(x[3])){
#       return(x[3])
#     }
# 
# 
#   })
# 
# kegg_id <-
#   kegg_id %>%
#   dplyr::mutate(KEGG_ID = KEGG_ID) %>%
#   dplyr::select(-c(KEGG_ID1, KEGG_ID2))
# 
# kegg_id
# 
# node_data <-
#   node_data %>%
#   dplyr::left_join(kegg_id, by = c("HMDB_ID"))
# 
# annotation_result <- node_data
# 
# save(annotation_result, file = "PIUMet/piumet_output_cluster3/Result/annotation_result")
# 
# 
# 
# ###cluster4
# load("PIUMet/piumet_output_cluster4/Result/node_data")
# hmdb_id <-
#   node_data %>%
#   dplyr::filter(!is.null(HMDB_ID)) %>%
#   dplyr::filter(HMDB_ID != " " & HMDB_ID != "NULL") %>%
#   dplyr::pull(HMDB_ID) %>%
#   unique()
# 
# kegg_id <-
#   lapply(hmdb_id, function(x){
#     metflow2::transID(query = x,
#                       from = "Human cytokine Database",
#                       to = "KEGG", top = 1)
#   }) %>%
#   do.call(rbind, .)
# 
# kegg_id
# 
# # load("PIUMet/hmdbMS1Database0.0.1")
# 
# hmdb_data <- hmdbMS1Database0.0.1@spectra.info
# 
# kegg_id2 <-
#   hmdb_data$KEGG.ID[match(kegg_id$`Human cytokine Database`, hmdb_data$HMDB.ID)]
# 
# kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)
# 
# colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")
# 
# KEGG_ID <-
#   apply(kegg_id, 1, function(x){
#     x <- as.character(x)
#     if(is.na(x[2]) & is.na(x[3])){
#       return(NA)
#     }
# 
#     if(!is.na(x[2])){
#       return(x[2])
#     }
# 
#     if(!is.na(x[3])){
#       return(x[3])
#     }
# 
# 
#   })
# 
# kegg_id <-
#   kegg_id %>%
#   dplyr::mutate(KEGG_ID = KEGG_ID) %>%
#   dplyr::select(-c(KEGG_ID1, KEGG_ID2))
# 
# kegg_id
# 
# node_data <-
#   node_data %>%
#   dplyr::left_join(kegg_id, by = c("HMDB_ID"))
# 
# annotation_result <- node_data
# 
# save(annotation_result, file = "PIUMet/piumet_output_cluster4/Result/annotation_result")
# 
# 
# ###cluster5
# load("PIUMet/piumet_output_cluster5/Result/node_data")
# 
# hmdb_id <-
#   node_data %>%
#   dplyr::filter(!is.null(HMDB_ID)) %>%
#   dplyr::filter(HMDB_ID != " " & HMDB_ID != "NULL") %>%
#   dplyr::pull(HMDB_ID) %>%
#   unique()
# 
# kegg_id <-
#   lapply(hmdb_id, function(x){
#     metflow2::transID(query = x,
#                       from = "Human cytokine Database",
#                       to = "KEGG", top = 1)
#   }) %>%
#   do.call(rbind, .)
# 
# kegg_id
# 
# # load("PIUMet/hmdbMS1Database0.0.1")
# 
# hmdb_data <- hmdbMS1Database0.0.1@spectra.info
# 
# kegg_id2 <-
#   hmdb_data$KEGG.ID[match(kegg_id$`Human cytokine Database`, hmdb_data$HMDB.ID)]
# 
# kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)
# 
# colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")
# 
# KEGG_ID <-
#   apply(kegg_id, 1, function(x){
#     x <- as.character(x)
#     if(is.na(x[2]) & is.na(x[3])){
#       return(NA)
#     }
# 
#     if(!is.na(x[2])){
#       return(x[2])
#     }
# 
#     if(!is.na(x[3])){
#       return(x[3])
#     }
# 
# 
#   })
# 
# kegg_id <-
#   kegg_id %>%
#   dplyr::mutate(KEGG_ID = KEGG_ID) %>%
#   dplyr::select(-c(KEGG_ID1, KEGG_ID2))
# 
# kegg_id
# 
# node_data <-
#   node_data %>%
#   dplyr::left_join(kegg_id, by = c("HMDB_ID"))
# 
# annotation_result <- node_data
# 
# save(annotation_result, file = "PIUMet/piumet_output_cluster5/Result/annotation_result")
# 
# ###cluster6
# load("PIUMet/piumet_output_cluster6/Result/node_data")
# 
# hmdb_id <-
#   node_data %>%
#   dplyr::filter(!is.null(HMDB_ID)) %>%
#   dplyr::filter(HMDB_ID != " " & HMDB_ID != "NULL") %>%
#   dplyr::pull(HMDB_ID) %>%
#   unique()
# 
# kegg_id <-
#   lapply(hmdb_id, function(x){
#     metflow2::transID(query = x,
#                       from = "Human cytokine Database",
#                       to = "KEGG", top = 1)
#   }) %>%
#   do.call(rbind, .)
# 
# kegg_id
# 
# # load("PIUMet/hmdbMS1Database0.0.1")
# 
# hmdb_data <- hmdbMS1Database0.0.1@spectra.info
# 
# kegg_id2 <-
#   hmdb_data$KEGG.ID[match(kegg_id$`Human cytokine Database`, hmdb_data$HMDB.ID)]
# 
# kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)
# 
# colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")
# 
# KEGG_ID <-
#   apply(kegg_id, 1, function(x){
#     x <- as.character(x)
#     if(is.na(x[2]) & is.na(x[3])){
#       return(NA)
#     }
# 
#     if(!is.na(x[2])){
#       return(x[2])
#     }
# 
#     if(!is.na(x[3])){
#       return(x[3])
#     }
# 
# 
#   })
# 
# kegg_id <-
#   kegg_id %>%
#   dplyr::mutate(KEGG_ID = KEGG_ID) %>%
#   dplyr::select(-c(KEGG_ID1, KEGG_ID2))
# 
# kegg_id
# 
# node_data <-
#   node_data %>%
#   dplyr::left_join(kegg_id, by = c("HMDB_ID"))
# 
# annotation_result <- node_data
# 
# save(annotation_result, file = "PIUMet/piumet_output_cluster6/Result/annotation_result")
# 
# 
# ###cluster7
# load("PIUMet/piumet_output_cluster7/Result/node_data")
# hmdb_id <-
#   node_data %>%
#   dplyr::filter(!is.null(HMDB_ID)) %>%
#   dplyr::filter(HMDB_ID != " " & HMDB_ID != "NULL") %>%
#   dplyr::pull(HMDB_ID) %>%
#   unique()
# 
# kegg_id <-
#   lapply(hmdb_id, function(x){
#     metflow2::transID(query = x,
#                       from = "Human cytokine Database",
#                       to = "KEGG", top = 1)
#   }) %>%
#   do.call(rbind, .)
# 
# kegg_id
# 
# # load("PIUMet/hmdbMS1Database0.0.1")
# 
# hmdb_data <- hmdbMS1Database0.0.1@spectra.info
# 
# kegg_id2 <-
#   hmdb_data$KEGG.ID[match(kegg_id$`Human cytokine Database`, hmdb_data$HMDB.ID)]
# 
# kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)
# 
# colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")
# 
# KEGG_ID <-
#   apply(kegg_id, 1, function(x){
#     x <- as.character(x)
#     if(is.na(x[2]) & is.na(x[3])){
#       return(NA)
#     }
# 
#     if(!is.na(x[2])){
#       return(x[2])
#     }
# 
#     if(!is.na(x[3])){
#       return(x[3])
#     }
# 
# 
#   })
# 
# kegg_id <-
#   kegg_id %>%
#   dplyr::mutate(KEGG_ID = KEGG_ID) %>%
#   dplyr::select(-c(KEGG_ID1, KEGG_ID2))
# 
# kegg_id
# 
# node_data <-
#   node_data %>%
#   dplyr::left_join(kegg_id, by = c("HMDB_ID"))
# 
# annotation_result <- node_data
# 
# save(annotation_result, file = "PIUMet/piumet_output_cluster7/Result/annotation_result")


######pathway enrichment
######
######
######
##cluster1
load("PIUMet/piumet_output_cluster1/Result/annotation_result")
load("PIUMet/piumet_output_cluster1/Result/edge_data")

edge_data <-
  edge_data %>%
  dplyr::filter(stringr::str_detect(from, "POS") |
                  stringr::str_detect(from, "NEG") |
                  stringr::str_detect(to, "POS") |
                  stringr::str_detect(to, "NEG")
  ) %>%
  dplyr::select(from, to) %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(.f = function(x){
    if(stringr::str_detect(x[1], "POS") | stringr::str_detect(x[1], "NEG")){
      return(x)
    }else{
      return(rev(x))
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::arrange(V1) %>%
  dplyr::distinct()

colnames(edge_data) <- c('peak', "metabolite")

annotation_result <-
  annotation_result %>%
  dplyr::filter(node_class == "Metabolite") %>%
  dplyr::select(node, HMDB_ID, KEGG_ID, super.class)

annotation_result_cluster1 <-
  edge_data %>%
  dplyr::left_join(annotation_result, by = c("metabolite" = "node"))

# save(annotation_result_cluster1,
#      file = "PIUMet/piumet_output_cluster1/Result/annotation_result_cluster1")
load("PIUMet/piumet_output_cluster1/Result/annotation_result_cluster1")


##cluster2
load("PIUMet/piumet_output_cluster2/Result/annotation_result")
load("PIUMet/piumet_output_cluster2/Result/edge_data")

edge_data <-
  edge_data %>%
  dplyr::filter(stringr::str_detect(from, "POS") |
                  stringr::str_detect(from, "NEG") |
                  stringr::str_detect(to, "POS") |
                  stringr::str_detect(to, "NEG")
  ) %>%
  dplyr::select(from, to) %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(.f = function(x){
    if(stringr::str_detect(x[1], "POS") | stringr::str_detect(x[1], "NEG")){
      return(x)
    }else{
      return(rev(x))
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::arrange(V1) %>%
  dplyr::distinct()

colnames(edge_data) <- c('peak', "metabolite")

annotation_result <-
  annotation_result %>%
  dplyr::filter(node_class == "Metabolite") %>%
  dplyr::select(node, HMDB_ID, KEGG_ID, super.class)

annotation_result_cluster2 <-
  edge_data %>%
  dplyr::left_join(annotation_result, by = c("metabolite" = "node"))

# save(annotation_result_cluster2,
#      file = "PIUMet/piumet_output_cluster2/Result/annotation_result_cluster2")
load("PIUMet/piumet_output_cluster2/Result/annotation_result_cluster2")

##cluster3
load("PIUMet/piumet_output_cluster3/Result/annotation_result")
load("PIUMet/piumet_output_cluster3/Result/edge_data")

edge_data <-
  edge_data %>%
  dplyr::filter(stringr::str_detect(from, "POS") |
                  stringr::str_detect(from, "NEG") |
                  stringr::str_detect(to, "POS") |
                  stringr::str_detect(to, "NEG")
  ) %>%
  dplyr::select(from, to) %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(.f = function(x){
    if(stringr::str_detect(x[1], "POS") | stringr::str_detect(x[1], "NEG")){
      return(x)
    }else{
      return(rev(x))
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::arrange(V1) %>%
  dplyr::distinct()

colnames(edge_data) <- c('peak', "metabolite")

annotation_result <-
  annotation_result %>%
  dplyr::filter(node_class == "Metabolite") %>%
  dplyr::select(node, HMDB_ID, KEGG_ID, super.class)

annotation_result_cluster3 <-
  edge_data %>%
  dplyr::left_join(annotation_result, by = c("metabolite" = "node"))

# save(annotation_result_cluster3,
#      file = "PIUMet/piumet_output_cluster3/Result/annotation_result_cluster3")
load("PIUMet/piumet_output_cluster3/Result/annotation_result_cluster3")


##cluster4
load("PIUMet/piumet_output_cluster4/Result/annotation_result")
load("PIUMet/piumet_output_cluster4/Result/edge_data")

edge_data <-
  edge_data %>%
  dplyr::filter(stringr::str_detect(from, "POS") |
                  stringr::str_detect(from, "NEG") |
                  stringr::str_detect(to, "POS") |
                  stringr::str_detect(to, "NEG")
  ) %>%
  dplyr::select(from, to) %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(.f = function(x){
    if(stringr::str_detect(x[1], "POS") | stringr::str_detect(x[1], "NEG")){
      return(x)
    }else{
      return(rev(x))
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::arrange(V1) %>%
  dplyr::distinct()

colnames(edge_data) <- c('peak', "metabolite")

annotation_result <-
  annotation_result %>%
  dplyr::filter(node_class == "Metabolite") %>%
  dplyr::select(node, HMDB_ID, KEGG_ID, super.class)

annotation_result_cluster4 <-
  edge_data %>%
  dplyr::left_join(annotation_result, by = c("metabolite" = "node"))

# save(annotation_result_cluster4,
#      file = "PIUMet/piumet_output_cluster4/Result/annotation_result_cluster4")
load("PIUMet/piumet_output_cluster4/Result/annotation_result_cluster4")

##cluster5
load("PIUMet/piumet_output_cluster5/Result/annotation_result")
load("PIUMet/piumet_output_cluster5/Result/edge_data")

edge_data <-
  edge_data %>%
  dplyr::filter(stringr::str_detect(from, "POS") |
                  stringr::str_detect(from, "NEG") |
                  stringr::str_detect(to, "POS") |
                  stringr::str_detect(to, "NEG")
  ) %>%
  dplyr::select(from, to) %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(.f = function(x){
    if(stringr::str_detect(x[1], "POS") | stringr::str_detect(x[1], "NEG")){
      return(x)
    }else{
      return(rev(x))
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::arrange(V1) %>%
  dplyr::distinct()

colnames(edge_data) <- c('peak', "metabolite")

annotation_result <-
  annotation_result %>%
  dplyr::filter(node_class == "Metabolite") %>%
  dplyr::select(node, HMDB_ID, KEGG_ID, super.class)

annotation_result_cluster5 <-
  edge_data %>%
  dplyr::left_join(annotation_result, by = c("metabolite" = "node"))

# save(annotation_result_cluster5,
#      file = "PIUMet/piumet_output_cluster5/Result/annotation_result_cluster5")
load("PIUMet/piumet_output_cluster5/Result/annotation_result_cluster5")

##cluster6
load("PIUMet/piumet_output_cluster6/Result/annotation_result")
load("PIUMet/piumet_output_cluster6/Result/edge_data")

edge_data <-
  edge_data %>%
  dplyr::filter(stringr::str_detect(from, "POS") |
                  stringr::str_detect(from, "NEG") |
                  stringr::str_detect(to, "POS") |
                  stringr::str_detect(to, "NEG")
  ) %>%
  dplyr::select(from, to) %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(.f = function(x){
    if(stringr::str_detect(x[1], "POS") | stringr::str_detect(x[1], "NEG")){
      return(x)
    }else{
      return(rev(x))
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::arrange(V1) %>%
  dplyr::distinct()

colnames(edge_data) <- c('peak', "metabolite")

annotation_result <-
  annotation_result %>%
  dplyr::filter(node_class == "Metabolite") %>%
  dplyr::select(node, HMDB_ID, KEGG_ID, super.class)

annotation_result_cluster6 <-
  edge_data %>%
  dplyr::left_join(annotation_result, by = c("metabolite" = "node"))

# save(annotation_result_cluster6,
#      file = "PIUMet/piumet_output_cluster6/Result/annotation_result_cluster6")

load("PIUMet/piumet_output_cluster6/Result/annotation_result_cluster6")

##cluster7
load("PIUMet/piumet_output_cluster7/Result/annotation_result")
load("PIUMet/piumet_output_cluster7/Result/edge_data")

edge_data <-
  edge_data %>%
  dplyr::filter(stringr::str_detect(from, "POS") |
                  stringr::str_detect(from, "NEG") |
                  stringr::str_detect(to, "POS") |
                  stringr::str_detect(to, "NEG")
  ) %>%
  dplyr::select(from, to) %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(.f = function(x){
    if(stringr::str_detect(x[1], "POS") | stringr::str_detect(x[1], "NEG")){
      return(x)
    }else{
      return(rev(x))
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::arrange(V1) %>%
  dplyr::distinct()

colnames(edge_data) <- c('peak', "metabolite")

annotation_result <-
  annotation_result %>%
  dplyr::filter(node_class == "Metabolite") %>%
  dplyr::select(node, HMDB_ID, KEGG_ID, super.class)

annotation_result_cluster7 <-
  edge_data %>%
  dplyr::left_join(annotation_result, by = c("metabolite" = "node"))

# save(annotation_result_cluster7,
#      file = "PIUMet/piumet_output_cluster7/Result/annotation_result_cluster7")
load("PIUMet/piumet_output_cluster7/Result/annotation_result_cluster7")

########functional annotation
load("PIUMet/hsa_pathway")

sxtTools::setwd_project()
source("R/2020_02_13/cytokine/pathway_enrichment.R")
setwd("data_analysis_2020_02_13/cytokine_analysis/k_means_clustring")

##cluster1
load("PIUMet/piumet_output_cluster1/Result/annotation_result_cluster1")

kegg_id <-
  annotation_result_cluster1$KEGG_ID

kegg_id <-
  kegg_id[!is.na(kegg_id)]

enrichment_kegg_cluster1 <-
  enrich_pathway(id = kegg_id,
                 pathway_database = hsa_pathway)

# save(enrichment_kegg_cluster1, file = "PIUMet/piumet_output_cluster1/Result/enrichment_kegg_cluster1")
load("PIUMet/piumet_output_cluster1/Result/enrichment_kegg_cluster1")

##cluster2
load("PIUMet/piumet_output_cluster2/Result/annotation_result_cluster2")
# hmdb_id <-
#   annotation_result_cluster2$HMDB_ID
#
# hmdb_id <-
#   hmdb_id[!is.na(hmdb_id)]

kegg_id <-
  annotation_result_cluster2$KEGG_ID

kegg_id <-
  kegg_id[!is.na(kegg_id)]

enrichment_kegg_cluster2 <-
  enrich_pathway(id = kegg_id,
                 pathway_database = hsa_pathway)


# save(enrichment_kegg_cluster2, file = "PIUMet/piumet_output_cluster2/Result/enrichment_kegg_cluster2")
load("PIUMet/piumet_output_cluster2/Result/enrichment_kegg_cluster2")

##cluster3
load("PIUMet/piumet_output_cluster3/Result/annotation_result_cluster3")
# hmdb_id <-
#   annotation_result_cluster3$HMDB_ID
#
# hmdb_id <-
#   hmdb_id[!is.na(hmdb_id)]
# 
# kegg_id <-
#   annotation_result_cluster3$KEGG_ID
# 
# kegg_id <-
#   kegg_id[!is.na(kegg_id)]
# 
# enrichment_kegg_cluster3 <-
#   enrich_pathway(id = kegg_id,
#                  pathway_database = hsa_pathway)
# 
# save(enrichment_kegg_cluster3, file = "PIUMet/piumet_output_cluster3/Result/enrichment_kegg_cluster3")

load("PIUMet/piumet_output_cluster3/Result/enrichment_kegg_cluster3")

##cluster4
load("PIUMet/piumet_output_cluster4/Result/annotation_result_cluster4")
# hmdb_id <-
#   annotation_result_cluster4$HMDB_ID
#
# hmdb_id <-
#   hmdb_id[!is.na(hmdb_id)] 
#             
# kegg_id <-
#   annotation_result_cluster4$KEGG_ID
# 
# kegg_id <-
#   kegg_id[!is.na(kegg_id)]
# 
# enrichment_kegg_cluster4 <-
#   enrich_pathway(id = kegg_id,
#                  pathway_database = hsa_pathway)
# 
# save(enrichment_kegg_cluster4, file = "PIUMet/piumet_output_cluster4/Result/enrichment_kegg_cluster4")
load("PIUMet/piumet_output_cluster4/Result/enrichment_kegg_cluster4")

##cluster5
load("PIUMet/piumet_output_cluster5/Result/annotation_result_cluster5")
# hmdb_id <-
#   annotation_result_cluster5$HMDB_ID
#
# hmdb_id <-
#   hmdb_id[!is.na(hmdb_id)]
# 
# kegg_id <-
#   annotation_result_cluster5$KEGG_ID
# 
# kegg_id <-
#   kegg_id[!is.na(kegg_id)]
# 
# enrichment_kegg_cluster5 <-
#   enrich_pathway(id = kegg_id,
#                  pathway_database = hsa_pathway)
# 
# 
# save(enrichment_kegg_cluster5, file = "PIUMet/piumet_output_cluster5/Result/enrichment_kegg_cluster5")
load("PIUMet/piumet_output_cluster5/Result/enrichment_kegg_cluster5")


##cluster6
load("PIUMet/piumet_output_cluster6/Result/annotation_result_cluster6")
# hmdb_id <-
#   annotation_result_cluster6$HMDB_ID
#
# hmdb_id <-
#   hmdb_id[!is.na(hmdb_id)]
# 
# kegg_id <-
#   annotation_result_cluster6$KEGG_ID
# 
# kegg_id <-
#   kegg_id[!is.na(kegg_id)]
# 
# enrichment_kegg_cluster6 <-
#   enrich_pathway(id = kegg_id,
#                  pathway_database = hsa_pathway)
# 
# save(enrichment_kegg_cluster6, file = "PIUMet/piumet_output_cluster6/Result/enrichment_kegg_cluster6")
load("PIUMet/piumet_output_cluster6/Result/enrichment_kegg_cluster6")

##cluster7
load("PIUMet/piumet_output_cluster7/Result/annotation_result_cluster7")
# hmdb_id <-
#   annotation_result_cluster7$HMDB_ID
#
# hmdb_id <-
#   hmdb_id[!is.na(hmdb_id)]
# 
# kegg_id <-
#   annotation_result_cluster7$KEGG_ID
# 
# kegg_id <-
#   kegg_id[!is.na(kegg_id)]
# 
# enrichment_kegg_cluster7 <-
#   enrich_pathway(id = kegg_id,
#                  pathway_database = hsa_pathway)
# 
# 
# save(enrichment_kegg_cluster7, file = "PIUMet/piumet_output_cluster7/Result/enrichment_kegg_cluster7")
load("PIUMet/piumet_output_cluster7/Result/enrichment_kegg_cluster7")


load("PIUMet/piumet_output_cluster1/Result/enrichment_kegg_cluster1")
load("PIUMet/piumet_output_cluster2/Result/enrichment_kegg_cluster2")
load("PIUMet/piumet_output_cluster3/Result/enrichment_kegg_cluster3")
load("PIUMet/piumet_output_cluster4/Result/enrichment_kegg_cluster4")
load("PIUMet/piumet_output_cluster5/Result/enrichment_kegg_cluster5")
load("PIUMet/piumet_output_cluster6/Result/enrichment_kegg_cluster6")   
load("PIUMet/piumet_output_cluster7/Result/enrichment_kegg_cluster7")

enrichment_kegg_cluster1 <-
  enrichment_kegg_cluster1 %>%
  dplyr::arrange(p.value.fdr)
# dplyr::filter(Overlap >= 3) %>%
# head(10)

enrichment_kegg_cluster2 <-
  enrichment_kegg_cluster2 %>%
  dplyr::arrange(p.value.fdr)

enrichment_kegg_cluster3 <-
  enrichment_kegg_cluster3 %>%
  dplyr::arrange(p.value.fdr)

enrichment_kegg_cluster4 <-
  enrichment_kegg_cluster4 %>%
  dplyr::arrange(p.value.fdr)

enrichment_kegg_cluster5 <-
  enrichment_kegg_cluster5 %>%
  dplyr::arrange(p.value.fdr)

enrichment_kegg_cluster6 <-
  enrichment_kegg_cluster6 %>%
  dplyr::arrange(p.value.fdr)

enrichment_kegg_cluster7 <-
  enrichment_kegg_cluster7 %>%
  dplyr::arrange(p.value.fdr)

intersect(enrichment_kegg_cluster1$Pathway.name, 
          enrichment_kegg_cluster2$Pathway.name)

result_cluster <-
  list(
    enrichment_kegg_cluster1,
    enrichment_kegg_cluster2,
    enrichment_kegg_cluster3,
    enrichment_kegg_cluster4,
    enrichment_kegg_cluster5,
    enrichment_kegg_cluster6,
    enrichment_kegg_cluster7
  )

unique_name <- 
  unique(c(result_cluster[[1]]$Pathway.name,
           result_cluster[[2]]$Pathway.name,
           result_cluster[[3]]$Pathway.name,
           result_cluster[[4]]$Pathway.name,
           result_cluster[[5]]$Pathway.name,
           result_cluster[[6]]$Pathway.name,
           result_cluster[[7]]$Pathway.name)) %>% 
  sort()



for(x in unique_name){
  cat(x, "\n")
  idx1 <- which(result_cluster[[1]]$Pathway.name == x)
  idx2 <- which(result_cluster[[2]]$Pathway.name == x)
  idx3 <- which(result_cluster[[3]]$Pathway.name == x)
  idx4 <- which(result_cluster[[4]]$Pathway.name == x)
  idx5 <- which(result_cluster[[5]]$Pathway.name == x)
  idx6 <- which(result_cluster[[6]]$Pathway.name == x)
  idx7 <- which(result_cluster[[7]]$Pathway.name == x)
  
  temp_data <-
    rbind(
      data.frame(result_cluster[[1]][idx1, c(1:7)], class = rep("no", length(idx1)),cluster = rep(1, length(idx1))),
      data.frame(result_cluster[[2]][idx2, c(1:7)], class = rep("increase", length(idx2)),cluster = rep(2, length(idx2))),
      data.frame(result_cluster[[3]][idx3, c(1:7)], class = rep("increase", length(idx3)),cluster = rep(3, length(idx3))),
      data.frame(result_cluster[[4]][idx4, c(1:7)], class = rep("decrease", length(idx4)),cluster = rep(4, length(idx4))),
      data.frame(result_cluster[[5]][idx5, c(1:7)], class = rep("no", length(idx5)),cluster = rep(5, length(idx5))),
      data.frame(result_cluster[[6]][idx6, c(1:7)], class = rep("no", length(idx6)),cluster = rep(6, length(idx6))),
      data.frame(result_cluster[[7]][idx7, c(1:7)], class = rep("decrease", length(idx7)),cluster = rep(7, length(idx7)))
    ) %>%
    dplyr::arrange(p.value.bh) %>% 
    dplyr::filter(class != "no")
  
  if(nrow(temp_data) == 1){
    next()
  }
  
  if(all(c("decrease", "increase") %in% unique(temp_data$class))){
    temp_data <- 
      temp_data %>% 
      dplyr::filter(class != temp_data$class[1])
    
    ##remove wrong pathways in other cluster
    for(i in 1:nrow(temp_data)){
      # cat(i, "")
      temp_cluster <- temp_data$cluster[i]
      temp_id <- temp_data$Pathway.ID[i]
      result_cluster[[temp_cluster]] <-
        result_cluster[[temp_cluster]] %>%
        dplyr::filter(Pathway.ID != temp_id)
    }
  }
}







library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")

purrr::walk(
  .x = 1:length(result_cluster),
  .f = function(i) {
    addWorksheet(wb, sheetName = paste("Cluster", i, sep = ""))
    freezePane(wb = wb, sheet = i, firstRow = TRUE, firstCol = FALSE) 
    writeDataTable(wb, sheet = i, x = result_cluster[[i]],
                   colNames = TRUE, rowNames = FALSE)
  }
)

saveWorkbook(wb, "cluster_annotation.xlsx", overwrite = TRUE)





save(enrichment_kegg_cluster1, file = "PIUMet/piumet_output_cluster1/Result/enrichment_kegg_cluster1")
save(enrichment_kegg_cluster2, file = "PIUMet/piumet_output_cluster1/Result/enrichment_kegg_cluster2")
save(enrichment_kegg_cluster4, file = "PIUMet/piumet_output_cluster1/Result/enrichment_kegg_cluster3")
save(enrichment_kegg_cluster5, file = "PIUMet/piumet_output_cluster1/Result/enrichment_kegg_cluster4")
save(enrichment_kegg_cluster6, file = "PIUMet/piumet_output_cluster1/Result/enrichment_kegg_cluster6")
save(enrichment_kegg_cluster8, file = "PIUMet/piumet_output_cluster1/Result/enrichment_kegg_cluster8")
save(enrichment_kegg_cluster9, file = "PIUMet/piumet_output_cluster1/Result/enrichment_kegg_cluster9")


enrichment_kegg_cluster1 %>%
  dplyr::mutate(fdr = -log(p.value.fdr, 10)) 


plot <- 
  enrichment_kegg_cluster1 %>% 
  dplyr::mutate(Pathway.name = stringr::str_replace(Pathway.name, " - Homo sapiens \\(human\\)", "")) %>% 
  dplyr::mutate(class = case_when(
    p.value.fdr < 0.05 ~ "< 0.05",
    TRUE ~ "> 0.05"
  )) %>% 
  dplyr::mutate(fdr = -log(p.value.fdr, 10)) %>% 
  dplyr::arrange(fdr) %>% 
  dplyr::mutate(Pathway.ID = factor(Pathway.ID, levels = Pathway.ID)) %>% 
  tail(10) %>% 
  ggplot(aes(x = fdr, y = Pathway.ID)) +
  geom_point(aes(size = Overlap, fill = class), shape = 21) +
  geom_segment(aes(x = 0, xend = fdr, y = Pathway.ID, yend = Pathway.ID)) +
  geom_text(aes(x = fdr, y = Pathway.ID, label = Pathway.name)) +
  scale_size_continuous(range = c(2,10)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "-log10(FDR adjusted P-values)", y = "") + 
  scale_fill_manual(values = c(
    "< 0.05" = ggsci::pal_aaas()(n=10)[2],
    "> 0.05" = "black"
  )) +
  guides(fill = guide_legend(title = "FDR adjusted P value", 
                             override.aes = list(size = 7)),
         size = guide_legend(title = "Metabolite number")) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

ggsave(plot,
       filename = "PIUMet/piumet_output_cluster1/Result/cluster1_annotation.pdf",
       width = 7,
       height = 7)




plot <- 
  enrichment_kegg_cluster2 %>% 
  dplyr::mutate(Pathway.name = stringr::str_replace(Pathway.name, " - Homo sapiens \\(human\\)", "")) %>% 
  dplyr::mutate(class = case_when(
    p.value.fdr < 0.05 ~ "< 0.05",
    TRUE ~ "> 0.05"
  )) %>% 
  dplyr::mutate(fdr = -log(p.value.fdr, 10)) %>% 
  dplyr::arrange(fdr) %>% 
  dplyr::mutate(Pathway.ID = factor(Pathway.ID, levels = Pathway.ID)) %>% 
  tail(10) %>% 
  ggplot(aes(x = fdr, y = Pathway.ID)) +
  geom_point(aes(size = Overlap, fill = class), shape = 21) +
  geom_segment(aes(x = 0, xend = fdr, y = Pathway.ID, yend = Pathway.ID)) +
  geom_text(aes(x = fdr, y = Pathway.ID, label = Pathway.name)) +
  scale_size_continuous(range = c(2,10)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  labs(x = "-log10(FDR adjusted P-values)", y = "") + 
  scale_fill_manual(values = c(
    "< 0.05" = ggsci::pal_aaas()(n=10)[2],
    "> 0.05" = "black"
  )) +
  guides(fill = guide_legend(title = "FDR adjusted P value", 
                             override.aes = list(size = 7)),
         size = guide_legend(title = "Metabolite number")) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

ggsave(plot,
       filename = "PIUMet/piumet_output_cluster2/Result/cluster2_annotation.pdf",
       width = 7,
       height = 7)



plot <- 
  enrichment_kegg_cluster4 %>% 
  dplyr::mutate(Pathway.name = stringr::str_replace(Pathway.name, " - Homo sapiens \\(human\\)", "")) %>% 
  dplyr::mutate(class = case_when(
    p.value.fdr < 0.05 ~ "< 0.05",
    TRUE ~ "> 0.05"
  )) %>% 
  dplyr::mutate(fdr = -log(p.value.fdr, 10)) %>% 
  dplyr::arrange(fdr) %>% 
  dplyr::mutate(Pathway.ID = factor(Pathway.ID, levels = Pathway.ID)) %>% 
  tail(10) %>% 
  ggplot(aes(x = fdr, y = Pathway.ID)) +
  geom_point(aes(size = Overlap, fill = class), shape = 21) +
  geom_segment(aes(x = 0, xend = fdr, y = Pathway.ID, yend = Pathway.ID)) +
  geom_text(aes(x = fdr, y = Pathway.ID, label = Pathway.name)) +
  scale_size_continuous(range = c(2,10)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "-log10(FDR adjusted P-values)", y = "") + 
  scale_fill_manual(values = c(
    "< 0.05" = ggsci::pal_aaas()(n=10)[2],
    "> 0.05" = "black"
  )) +
  guides(fill = guide_legend(title = "FDR adjusted P value", 
                             override.aes = list(size = 7)),
         size = guide_legend(title = "Metabolite number", nrow = 1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

ggsave(plot,
       filename = "PIUMet/piumet_output_cluster4/Result/cluster4_annotation.pdf",
       width = 7,
       height = 7)


plot <- 
  enrichment_kegg_cluster5 %>% 
  dplyr::mutate(Pathway.name = stringr::str_replace(Pathway.name, " - Homo sapiens \\(human\\)", "")) %>% 
  dplyr::mutate(class = case_when(
    p.value.fdr < 0.05 ~ "< 0.05",
    TRUE ~ "> 0.05"
  )) %>% 
  dplyr::mutate(fdr = -log(p.value.fdr, 10)) %>% 
  dplyr::arrange(fdr) %>% 
  dplyr::mutate(Pathway.ID = factor(Pathway.ID, levels = Pathway.ID)) %>% 
  tail() %>% 
  ggplot(aes(x = fdr, y = Pathway.ID)) +
  geom_point(aes(size = Overlap, fill = class), shape = 21) +
  geom_segment(aes(x = 0, xend = fdr, y = Pathway.ID, yend = Pathway.ID)) +
  geom_text(aes(x = fdr, y = Pathway.ID, label = Pathway.name)) +
  scale_size_continuous(range = c(2,10)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "-log10(FDR adjusted P-values)", y = "") + 
  scale_fill_manual(values = c(
    "< 0.05" = ggsci::pal_aaas()(n=10)[2],
    "> 0.05" = "black"
  )) +
  guides(fill = guide_legend(title = "FDR adjusted P value", 
                             override.aes = list(size = 7)),
         size = guide_legend(title = "Metabolite number")) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

ggsave(plot,
       filename = "PIUMet/piumet_output_cluster5/Result/cluster5_annotation.pdf",
       width = 7,
       height = 7)


plot <- 
  enrichment_kegg_cluster6 %>% 
  dplyr::mutate(Pathway.name = stringr::str_replace(Pathway.name, " - Homo sapiens \\(human\\)", "")) %>% 
  dplyr::mutate(class = case_when(
    p.value.fdr < 0.05 ~ "< 0.05",
    TRUE ~ "> 0.05"
  )) %>% 
  dplyr::mutate(fdr = -log(p.value.fdr, 10)) %>% 
  dplyr::arrange(fdr) %>% 
  dplyr::mutate(Pathway.ID = factor(Pathway.ID, levels = Pathway.ID)) %>% 
  tail() %>% 
  ggplot(aes(x = fdr, y = Pathway.ID)) +
  geom_point(aes(size = Overlap, fill = class), shape = 21) +
  geom_segment(aes(x = 0, xend = fdr, y = Pathway.ID, yend = Pathway.ID)) +
  geom_text(aes(x = fdr, y = Pathway.ID, label = Pathway.name)) +
  scale_size_continuous(range = c(2,10)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "-log10(FDR adjusted P-values)", y = "") + 
  scale_fill_manual(values = c(
    "< 0.05" = ggsci::pal_aaas()(n=10)[2],
    "> 0.05" = "black"
  )) +
  guides(fill = guide_legend(title = "FDR adjusted P value", 
                             override.aes = list(size = 7)),
         size = guide_legend(title = "Metabolite number", nrow = 1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

ggsave(plot,
       filename = "PIUMet/piumet_output_cluster6/Result/cluster6_annotation.pdf",
       width = 7,
       height = 7)













