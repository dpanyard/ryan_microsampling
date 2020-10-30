#to avoid source
no_exist_function()
sxtTools::setwd_project()
rm(list = ls())
library(tidyverse)


###load data
###metabolome
load("data/shake_study/metabolome_data_analysis/invidual_variation/distance_individual")
metabolome_distance_individual <- distance_individual
load("data/shake_study/metabolome_data_analysis/invidual_variation/metabolome_icc")
load("data/shake_study/metabolome_data_analysis/invidual_variation/metabolome_iqr")
load("data/shake_study/metabolome_data_analysis/data_preparation/expression_data")
load("data/shake_study/metabolome_data_analysis/data_preparation/sample_info")
load("data/shake_study/metabolome_data_analysis/data_preparation/variable_info")
metabolome_expression_data <-
  log(expression_data + 1, 2)
metabolome_sample_info <- sample_info
metabolome_variable_info <- variable_info

###lipidomics
load("data/shake_study/lipidomics_data_anlaysis/invidual_variation/distance_individual")
lipidomics_distance_individual <- distance_individual
load("data/shake_study/lipidomics_data_anlaysis/invidual_variation/lipidomics_icc")
load("data/shake_study/lipidomics_data_anlaysis/invidual_variation/lipidomics_iqr")
load("data/shake_study/lipidomics_data_anlaysis/data_preparation/expression_data")
load("data/shake_study/lipidomics_data_anlaysis/data_preparation/sample_info")
load("data/shake_study/lipidomics_data_anlaysis/data_preparation/variable_info")
lipidomics_expression_data <-
  log(expression_data + 1, 2)
lipidomics_sample_info <- sample_info
lipidomics_variable_info <- variable_info

###cytokine
load("data/shake_study/cytokine_data_analysis/invidual_variation/distance_individual")
cytokine_distance_individual <- distance_individual
load("data/shake_study/cytokine_data_analysis/invidual_variation/cytokine_icc")
load("data/shake_study/cytokine_data_analysis/invidual_variation/cytokine_iqr")
load("data/shake_study/cytokine_data_analysis/data_preparation/expression_data")
load("data/shake_study/cytokine_data_analysis/data_preparation/sample_info")
load("data/shake_study/cytokine_data_analysis/data_preparation/variable_info")
cytokine_expression_data <-
  log(expression_data + 1, 2)
cytokine_sample_info <- sample_info
cytokine_variable_info <- variable_info




setwd("data/shake_study/individual_variation")


####inter distance
metabolome_distance_individual <-
  metabolome_distance_individual %>%
  do.call(cbind, .) %>%
  data.frame(time = c("30-0",
                      "60-30",
                      "120-60",
                      "240-120")) %>%
  tidyr::pivot_longer(cols = -time,
                      names_to = "subject_id",
                      values_to = "value") %>%
  dplyr::mutate(time = factor(time, levels = c("30-0",
                                               "60-30",
                                               "120-60",
                                               "240-120"))) %>% 
  dplyr::group_by(subject_id) %>%
  dplyr::summarise(mean_distance = mean(value, na.rm = TRUE)) %>%
  dplyr::ungroup()


lipidomics_distance_individual <-
  lipidomics_distance_individual %>%
  do.call(cbind, .) %>%
  data.frame(time = c("30-0",
                      "60-30",
                      "120-60",
                      "240-120")) %>%
  tidyr::pivot_longer(cols = -time,
                      names_to = "subject_id",
                      values_to = "value") %>%
  dplyr::mutate(time = factor(time, levels = c("30-0",
                                               "60-30",
                                               "120-60",
                                               "240-120"))) %>% 
  dplyr::group_by(subject_id) %>%
  dplyr::summarise(mean_distance = mean(value, na.rm = TRUE)) %>%
  dplyr::ungroup()


cytokine_distance_individual <-
  cytokine_distance_individual %>%
  do.call(cbind, .) %>%
  data.frame(time = c("30-0",
                      "60-30",
                      "120-60",
                      "240-120")) %>%
  tidyr::pivot_longer(cols = -time,
                      names_to = "subject_id",
                      values_to = "value") %>%
  dplyr::mutate(time = factor(time, levels = c("30-0",
                                               "60-30",
                                               "120-60",
                                               "240-120"))) %>% 
  dplyr::group_by(subject_id) %>%
  dplyr::summarise(mean_distance = mean(value, na.rm = TRUE)) %>%
  dplyr::ungroup()


###find the top 10 variable individuals
distance_individual <-
  metabolome_distance_individual %>%
  dplyr::rename(metabolome_distance = mean_distance) %>%
  dplyr::left_join(
    lipidomics_distance_individual %>%
      dplyr::rename(lipidomics_distance = mean_distance),
    by = c("subject_id")
  ) %>%
  dplyr::left_join(
    cytokine_distance_individual %>%
      dplyr::rename(cytokine_distance = mean_distance),
    by = c("subject_id")
  )

mean_distance_individual <-
  apply(distance_individual, 1, function(x) {
    mean(as.numeric(x[-1]), na.rm = TRUE)
  })

names(mean_distance_individual) <- distance_individual$subject_id
top10_id <-
  sort(mean_distance_individual) %>% tail(10) %>% names()

subject_color <-
  c(ggsci::pal_d3()(n = 10)[c(1:7, 9:10)],
    ggsci::pal_aaas()(n = 10)[4], "grey")

names(subject_color) <- c(top10_id, "Other")


distance_plot <-
  distance_individual %>%
  tidyr::pivot_longer(cols = -subject_id,
                      names_to = "data",
                      values_to = "value") %>%
  dplyr::mutate(data = stringr::str_replace(data, "_distance", "")) %>%
  dplyr::mutate(data = factor(
    data,
    levels = c(
      "metabolome",
      "lipidomics",
      "cytokine"
    )
  )) %>%
  dplyr::mutate(top10 = case_when(subject_id %in% top10_id ~ subject_id,
                                  TRUE ~ "Other")) %>%
  ggplot(aes(data, value, group = subject_id)) +
  geom_point(aes(fill = top10), shape = 21, size = 3) +
  geom_line(aes(color = top10, group = subject_id), alpha = 0.7) +
  scale_fill_manual(values = subject_color) +
  scale_color_manual(values = subject_color) +
  theme_bw() +
  labs(x = "", y = "Scaled Euclidean Distance") +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(
      size = 12,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "right",
    # legend.justification = c(0,1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 15),
    panel.grid.minor = element_blank()
  ) +
  ggforce::facet_zoom(ylim = c(1, 7.5), zoom.size = 1)


distance_plot


ggsave(distance_plot, filename = "distance_plot.pdf", width = 14, height = 7)


####ICC
icc_data <-
  rbind(
    metabolome_icc %>%
      data.frame(
        value = .,
        data = "Metabolome",
        stringsAsFactors = FALSE
      ) %>%
      tibble::rownames_to_column(var = "variable_id"),
    lipidomics_icc %>%
      data.frame(
        value = .,
        data = "Lipidomics",
        stringsAsFactors = FALSE
      ) %>%
      tibble::rownames_to_column(var = "variable_id"),
    cytokine_icc %>%
      data.frame(
        value = .,
        data = "Cytokine",
        stringsAsFactors = FALSE
      ) %>%
      tibble::rownames_to_column(var = "variable_id")
  ) %>%
  dplyr::mutate(data = factor(data, levels = unique(data)))


value <-
  c(
    "Transcriptome" = ggsci::pal_aaas()(10)[1],
    "Proteome" = ggsci::pal_aaas()(10)[2],
    "Metabolome" = ggsci::pal_aaas()(10)[3],
    "Cytokine" = ggsci::pal_aaas()(10)[4],
    "Vaginal microbiome" = ggsci::pal_aaas()(10)[5],
    "Clinical" = ggsci::pal_aaas()(10)[6],
    "Lipidomics" = ggsci::pal_aaas()(10)[8]
    # "Phenotye" = ggsci::pal_aaas()(10)[7]
  )

plot <-
  icc_data %>%
  as.data.frame() %>%
  ggplot(aes(x = data, y = value)) +
  geom_jitter(aes(size = data), alpha = 0.5, shape = 16) +
  geom_boxplot(fill = "transparent", aes(color = data), outlier.shape = NA) +
  scale_color_manual(values = value) +
  scale_size_manual(
    values = c(
      "Transcriptome" = 1,
      "Proteome" = 3,
      "Metabolome" = 1,
      "Cytokine" = 3,
      "Lipidomics" = 3,
      "Gut microbiome" = 2
    )
  ) +
  theme_bw() +
  labs(x = "", y = "Scaled Euclidean Distance") +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(
      size = 12,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "none",
    # legend.justification = c(0,1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 15),
    panel.grid.minor = element_blank()
  )

plot

ggsave(plot, filename = "ICC.pdf", width = 7, height = 7)

icc_data

write.csv(icc_data, "icc_data.csv", row.names = FALSE)


##for each person
######IQR
##metabolome
metabolome_head10_id <- head(metabolome_iqr %>% sort(decreasing = TRUE), 10)
metabolome_tail10_id <- tail(metabolome_iqr %>% sort(decreasing = TRUE), 10)

plot <- 
  metabolome_expression_data[names(metabolome_head10_id),] %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "variable_id") %>% 
  tidyr::pivot_longer(cols = -variable_id, names_to = "sample_id", values_to = "value") %>% 
  dplyr::mutate(variable_id = factor(variable_id, levels = names(metabolome_head10_id))) %>% 
  dplyr::left_join(metabolome_variable_info[,c("variable_id", "mol_name")], by = "variable_id") %>% 
  ggplot(aes(variable_id, value)) +
  geom_jitter(shape = 16, alpha = 0.8, size = 2) +
  geom_boxplot(fill = "transparent", color = ggsci::pal_aaas()(n=10)[2], outlier.shape = NA) +
  theme_bw() +
  labs(x = "", y = "log2(Intensity)") +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(
      size = 12,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "right",
    # legend.justification = c(0,1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 15),
    panel.grid.minor = element_blank()
  )

plot

ggsave(plot, filename = "metabolome_head10.pdf", width = 8, height = 7)

plot <- 
  metabolome_expression_data[names(metabolome_tail10_id),] %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "variable_id") %>% 
  tidyr::pivot_longer(cols = -variable_id, names_to = "sample_id", values_to = "value") %>% 
  dplyr::mutate(variable_id = factor(variable_id, levels = names(metabolome_tail10_id))) %>% 
  dplyr::left_join(metabolome_variable_info[,c("variable_id", "mol_name")], by = "variable_id") %>% 
  ggplot(aes(variable_id, value)) +
  geom_jitter(shape = 16, alpha = 0.8, size = 2) +
  geom_boxplot(fill = "transparent", color = ggsci::pal_aaas()(n=10)[1], outlier.shape = NA) +
  theme_bw() +
  labs(x = "", y = "log2(Intensity)") +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(
      size = 12,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "right",
    # legend.justification = c(0,1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 15),
    panel.grid.minor = element_blank()
  ) 
# scale_y_continuous(limits = c(0, 0.05)) 

plot

ggsave(plot, filename = "metabolome_tail10.pdf", width = 8, height = 7)


metabolome_iqr <- 
  metabolome_iqr %>% 
  data.frame(IQR = ., stringsAsFactors = FALSE) %>% 
  tibble::rownames_to_column(var = "variable_id") %>% 
  dplyr::left_join(metabolome_variable_info, 
                   by = c("variable_id"))

write.csv(metabolome_iqr, "metabolome_iqr.csv", row.names = FALSE)





##lipidomics
lipidomics_head10_id <- head(lipidomics_iqr %>% sort(decreasing = TRUE), 10)
lipidomics_tail10_id <- tail(lipidomics_iqr %>% sort(decreasing = TRUE), 10)

plot <- 
  lipidomics_expression_data[names(lipidomics_head10_id),] %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "variable_id") %>% 
  tidyr::pivot_longer(cols = -variable_id, names_to = "sample_id", values_to = "value") %>% 
  dplyr::mutate(variable_id = factor(variable_id, levels = names(lipidomics_head10_id))) %>% 
  dplyr::left_join(lipidomics_variable_info[,c("variable_id", "mol_name")], by = "variable_id") %>% 
  ggplot(aes(mol_name, value)) +
  geom_jitter(shape = 16, alpha = 0.8, size = 2) +
  geom_boxplot(fill = "transparent", color = ggsci::pal_aaas()(n=10)[2], outlier.shape = NA) +
  theme_bw() +
  labs(x = "", y = "log2(Intensity)") +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(
      size = 12,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "right",
    # legend.justification = c(0,1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 15),
    panel.grid.minor = element_blank()
  )

plot

ggsave(plot, filename = "lipidomics_head10.pdf", width = 8, height = 7)

plot <- 
  lipidomics_expression_data[names(lipidomics_tail10_id),] %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "variable_id") %>% 
  tidyr::pivot_longer(cols = -variable_id, names_to = "sample_id", values_to = "value") %>% 
  dplyr::mutate(variable_id = factor(variable_id, levels = names(lipidomics_tail10_id))) %>% 
  dplyr::left_join(lipidomics_variable_info[,c("variable_id", "mol_name")], by = "variable_id") %>% 
  ggplot(aes(mol_name, value)) +
  geom_jitter(shape = 16, alpha = 0.8, size = 2) +
  geom_boxplot(fill = "transparent", color = ggsci::pal_aaas()(n=10)[1], outlier.shape = NA) +
  theme_bw() +
  labs(x = "", y = "log2(Intensity)") +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(
      size = 12,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "right",
    # legend.justification = c(0,1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 15),
    panel.grid.minor = element_blank()
  ) 
# scale_y_continuous(limits = c(0, 0.05)) 

plot

ggsave(plot, filename = "lipidomics_tail10.pdf", width = 8, height = 7)


lipidomics_iqr <- 
  lipidomics_iqr %>% 
  data.frame(IQR = ., stringsAsFactors = FALSE) %>% 
  tibble::rownames_to_column(var = "variable_id") %>% 
  dplyr::left_join(lipidomics_variable_info, 
                   by = c("variable_id"))

write.csv(lipidomics_iqr, "lipidomics_iqr.csv", row.names = FALSE)


##cytokine
cytokine_head10_id <- head(cytokine_iqr %>% sort(decreasing = TRUE), 10)
cytokine_tail10_id <- tail(cytokine_iqr %>% sort(decreasing = TRUE), 10)

plot <- 
  cytokine_expression_data[names(cytokine_head10_id),] %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "variable_id") %>% 
  tidyr::pivot_longer(cols = -variable_id, names_to = "sample_id", values_to = "value") %>% 
  dplyr::mutate(variable_id = factor(variable_id, levels = names(cytokine_head10_id))) %>% 
  dplyr::left_join(cytokine_variable_info[,c("variable_id", "mol_name")], by = "variable_id") %>% 
  ggplot(aes(mol_name, value)) +
  geom_jitter(shape = 16, alpha = 0.8, size = 2) +
  geom_boxplot(fill = "transparent", color = ggsci::pal_aaas()(n=10)[2], outlier.shape = NA) +
  theme_bw() +
  labs(x = "", y = "log2(Intensity)") +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(
      size = 12,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "right",
    # legend.justification = c(0,1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 15),
    panel.grid.minor = element_blank()
  )

plot

ggsave(plot, filename = "cytokine_head10.pdf", width = 8, height = 7)

plot <- 
  cytokine_expression_data[names(cytokine_tail10_id),] %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "variable_id") %>% 
  tidyr::pivot_longer(cols = -variable_id, names_to = "sample_id", values_to = "value") %>% 
  dplyr::mutate(variable_id = factor(variable_id, levels = names(cytokine_tail10_id))) %>% 
  dplyr::left_join(cytokine_variable_info[,c("variable_id", "mol_name")], by = "variable_id") %>% 
  ggplot(aes(mol_name, value)) +
  geom_jitter(shape = 16, alpha = 0.8, size = 2) +
  geom_boxplot(fill = "transparent", color = ggsci::pal_aaas()(n=10)[1], outlier.shape = NA) +
  theme_bw() +
  labs(x = "", y = "log2(Intensity)") +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(
      size = 12,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "right",
    # legend.justification = c(0,1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 15),
    panel.grid.minor = element_blank()
  )

plot

ggsave(plot, filename = "cytokine_tail10.pdf", width = 8, height = 7)


cytokine_iqr1 <- 
  cytokine_iqr %>% 
  data.frame(IQR = ., stringsAsFactors = FALSE) %>% 
  tibble::rownames_to_column(var = "variable_id") %>% 
  dplyr::left_join(cytokine_variable_info, 
                   by = c("variable_id"))

write.csv(cytokine_iqr1, "cytokine_iqr.csv", row.names = FALSE)



#####z-score for each variable and each subject
library(plyr)
metabolome_rsd <-
  metabolome_expression_data %>%
  t() %>%
  as.data.frame() %>%
  dplyr::mutate(subject_id = metabolome_sample_info$subject_id) %>%
  plyr::dlply(.variables = .(subject_id)) %>% 
  purrr::map(function(x){
    rsd <- 
      apply(x %>% dplyr::select(-subject_id), 2, function(y){
        sd(y)/mean(y)
      })
    rsd[is.na(rsd)] <- 0
    sum(rsd > 0.3)
  }) %>% 
  unlist() %>% 
  `/`(nrow(metabolome_expression_data))


lipidomics_rsd <-
  lipidomics_expression_data %>%
  t() %>%
  as.data.frame() %>%
  dplyr::mutate(subject_id = lipidomics_sample_info$subject_id) %>%
  plyr::dlply(.variables = .(subject_id)) %>% 
  purrr::map(function(x){
    rsd <- 
      apply(x %>% dplyr::select(-subject_id), 2, function(y){
        sd(y)/mean(y)
      })
    rsd[is.na(rsd)] <- 0
    sum(rsd > 0.3)
  }) %>% 
  unlist() %>% 
  `/`(nrow(lipidomics_expression_data))

cytokine_rsd <-
  cytokine_expression_data %>%
  t() %>%
  as.data.frame() %>%
  dplyr::mutate(subject_id = cytokine_sample_info$subject_id) %>%
  plyr::dlply(.variables = .(subject_id)) %>% 
  purrr::map(function(x){
    rsd <- 
      apply(x %>% dplyr::select(-subject_id), 2, function(y){
        sd(y)/mean(y)
      })
    rsd[is.na(rsd)] <- 0
    sum(rsd > 0.3)
  }) %>% 
  unlist() %>% 
  `/`(nrow(cytokine_expression_data))

rsd_data <- 
  metabolome_rsd %>% 
  data.frame(metabolome = ., stringsAsFactors = FALSE) %>% 
  tibble::rownames_to_column(var = "variable_id") %>% 
  dplyr::left_join(
    lipidomics_rsd %>% 
      data.frame(lipidomics = ., stringsAsFactors = FALSE) %>% 
      tibble::rownames_to_column(var = "variable_id"),
    by = "variable_id"
  )  %>% 
  dplyr::left_join(
    cytokine_rsd %>% 
      data.frame(cytokine = ., stringsAsFactors = FALSE) %>% 
      tibble::rownames_to_column(var = "variable_id"),
    by = "variable_id"
  ) 

temp_data1 <- 
  rsd_data %>%
  apply(1, function(x) {
    sum(as.numeric(x[-1]), na.rm = TRUE)
  }) %>%
  `/`(6) %>%
  data.frame(value = .) %>%
  dplyr::mutate(subject_id = rsd_data$variable_id) %>%
  dplyr::arrange(value) %>%
  dplyr::mutate(subject_id = factor(subject_id, subject_id))

plot1 <-
  temp_data1 %>%
  ggplot(aes(x =  value * 100, y = subject_id)) +
  geom_bar(stat = "identity", fill = ggsci::pal_aaas()(n=10)[2],
           color = "black") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_bw() +
  labs(x = "Varying features percentage (%)", y = "") +
  theme(
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.position = "right",
    # legend.justification = c(0,1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 15),
    panel.grid.minor = element_blank()
  )

plot1

temp_data2 <-
  rsd_data %>% 
  tidyr::pivot_longer(cols = -variable_id, names_to = "data", values_to = "value") %>% 
  plyr::dlply(.variables = .(variable_id)) %>% 
  purrr::map(function(x){
    x$value[is.na(x$value)] <- 0
    x$value <- x$value/sum(x$value)
    x
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  dplyr::mutate(subject_id = 
                  factor(variable_id, levels = temp_data1$subject_id)) %>% 
  dplyr::mutate(
    data =
      case_when(
        data == "transcriptome" ~ "Transcriptome",
        data == "proteome" ~ "Proteome",
        data == "metabolome" ~ "Metabolome",
        data == "cytokine" ~ "Cytokine",
        data == "microbiome" ~ "Vaginal microbiome",
        data == "lipidomics" ~ "Lipidomics"
        # TRUE ~ NA
      )
  ) %>% 
  dplyr::mutate(data = factor(data, levels = c("Metabolome",
                                               "Lipidomics",
                                               "Cytokine") %>% rev()))

plot2 <-
  temp_data2 %>%
  ggplot(aes(x = value, y = subject_id)) +
  geom_bar(position = "stack", stat = "identity", aes(fill = data),
           show.legend = FALSE) +
  scale_fill_manual(values = value) +
  scale_x_continuous(expand = expansion(mult = c(0,0))) +
  scale_y_discrete(expand = expansion(mult = c(0,0))) +
  theme_bw() +
  labs(x = "Varying features percentage (%)", y = "") +
  theme(
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.position = "right",
    # legend.justification = c(0,1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 15),
    panel.grid.minor = element_blank()
  )

plot2

plot3 <-
  temp_data2 %>%
  dplyr::mutate(class = case_when(value == 0 ~ "no",
                                  TRUE ~ "yes")) %>%
  ggplot(aes(x = data, y = subject_id)) +
  geom_point(aes(fill = data,
                 shape = class),
             size = 4) +
  scale_shape_manual(values = c("no" = NA,
                                "yes" = 21)) +
  scale_fill_manual(values = value)  +
  guides(fill = guide_legend(title = "", override.aes = list(size = 4, shape = 21))) +
  theme_bw() +
  labs(x = "Varying features percentage (%)", y = "") +
  theme(
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.position = "right",
    # legend.justification = c(0,1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 15),
    panel.grid.minor = element_blank()
  )

plot3

library(patchwork)

plot <- 
  plot1 + plot2 + plot3 + 
  patchwork::plot_layout(nrow = 1, widths = c(1, 4, 2))

plot

ggsave(plot, filename = "percentage_plot.pdf", width = 8, height = 10)

