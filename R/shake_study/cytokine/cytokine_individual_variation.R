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
setwd("data/shake_study/cytokine_data_analysis/invidual_variation")

dim(expression_data)
dim(sample_info)
dim(variable_info)


###calculate distance for each person
##for cytokine, we use the euclidean distance
plot(density(as.numeric(expression_data[1,])))

expression_data1 <- log(expression_data + 1, 2)
plot(density(as.numeric(expression_data1[1,])))

distance_individual <-
  purrr::map(
    unique(sample_info$subject_id),
    .f = function(x) {
      temp_sample_info <- sample_info %>%
        dplyr::filter(subject_id == x) %>%
        dplyr::arrange(TP)
      
      temp_data <- expression_data1 %>%
        dplyr::select(temp_sample_info$sample_id)
      
      tp_list <-
        c(0, 30, 60, 120, 240)
      
      distance <-
        purrr::map(
          2:length(tp_list),
          .f = function(idx) {
            temp_tp <- tp_list[c(idx - 1, idx)]
            temp_idx <-
            which(temp_sample_info$TP %in% temp_tp)
            if (length(temp_idx) != 2) {
              return(NA)
            } else{
              dist(x = t(as.matrix(temp_data[, temp_idx])), method = "euclidean") %>%
                as.numeric()
            }
          }
        ) %>%
        unlist()
      
      distance
      
    }
  )

names(distance_individual) <- unique(sample_info$subject_id)

save(distance_individual, file = "distance_individual")
load("distance_individual")

tp_list <- c(0, 30, 60, 120, 240)

distance_individual <-
  distance_individual %>%
  do.call(cbind, .) %>%
  data.frame(
    time = c(
      "30-0",
      "60-30",
      "120-60",
      "240-120"
    )
  ) %>%
  tidyr::pivot_longer(cols = -time,
                      names_to = "subject_id",
                      values_to = "value") %>% 
  dplyr::mutate(time = factor(time, levels = c(
    "30-0",
    "60-30",
    "120-60",
    "240-120"
  )))

length(unique(sample_info$subject_id))

subject_col <-
  colorRampPalette(colors = RColorBrewer::brewer.pal(11, name = "Accent"))(n = 28)

names(subject_col) <-
  stringr::str_sort(unique(sample_info$subject_id), numeric = TRUE)

distance_individual %>% 
ggplot(aes(time, value, group = subject_id)) +
  geom_point(aes(fill = subject_id, group = subject_id), shape = 21) +
  geom_line(aes(color = subject_id, group = subject_id)) +
  scale_fill_manual(values = subject_col) +
  scale_color_manual(values = subject_col) +
  theme_bw() +
  labs(x = "Time", y = "Euclidean Distance") +
  theme_bw() +
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
    legend.position = "bottom",
    # legend.justification = c(0,1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 15),
    panel.grid.minor = element_blank()
  )


plot1 <-
  distance_individual %>%
  dplyr::mutate(subject_id = factor(subject_id, levels = stringr::str_sort(unique(subject_id), numeric = TRUE, decreasing = TRUE))) %>% 
  dplyr::mutate(value = value / max(value, na.rm = TRUE)) %>%
  ggplot(aes(time, value, group = subject_id)) +
  geom_tile(aes(x = time, y = subject_id, fill = value), color = "white") +
  labs(x = "Time", y = "") +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  geom_text(aes(
    x = time,
    y = subject_id,
    label = round(value, 2)
  ), color = "white") +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    # axis.text.x = element_text(
    #   size = 12,
    #   angle = 45,
    #   hjust = 1,
    #   vjust = 1
    # ),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "left",
    # legend.justification = c(0,1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 15),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    )
  )

plot1


plot2 <-
  distance_individual %>%
  dplyr::mutate(subject_id = factor(subject_id, levels = stringr::str_sort(unique(subject_id), numeric = TRUE, decreasing = TRUE))) %>% 
  dplyr::mutate(value = value / max(value, na.rm = TRUE)) %>%
  dplyr::group_by(subject_id) %>%
  dplyr::summarise(mean = mean(value, na.rm = TRUE)) %>%
  ggplot(aes(x = mean, y = subject_id)) +
  geom_segment(aes(
    x = 0,
    y = subject_id,
    xend = mean,
    yend = subject_id
  )) +
  geom_point(aes(fill = mean), shape = 21, size = 3) +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_classic() +
  labs(x = "Average Scaled Euclidean Distance", y = "") +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    axis.text.y = element_blank(),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "none",
    # legend.justification = c(0,1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 15),
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    )
  )

plot2

library(patchwork)

plot <-
  plot1 + plot2 + plot_layout(ncol = 2, widths = c(4, 1))

plot

ggsave(plot, filename = "cytokine_inter_vist_distance.pdf", width = 9, height = 7)




####intra-class cor- relation (ICC) analysis 
dim(expression_data1)

temp_data <- data.frame(
  subject_id = sample_info$subject_id,
  t(expression_data1),
  stringsAsFactors = TRUE
)

variable_id <- variable_info$variable_id

colnames(temp_data)[-1] <-
  variable_id

library(future)
plan(multiprocess)
library(ICC)
cytokine_icc <-
  furrr::future_map(.x = variable_id, function(x) {
    icc <-
      ICCest(
        x = subject_id,
        y = x,
        data = temp_data,
        alpha = 0.05,
        CI.type = "S"
      )
    icc$ICC
  }) %>%
  unlist()

names(cytokine_icc) <- variable_id
save(cytokine_icc, file = "cytokine_icc")

load("cytokine_icc")

##IQR for each metabolie
cytokine_iqr <-
  furrr::future_map(.x = variable_id, function(x) {
  x <- temp_data %>%
      pull(x)
  IQR(x)
  }) %>%
  unlist()

names(cytokine_iqr) <- variable_id
save(cytokine_iqr, file = "cytokine_iqr")

load("cytokine_iqr")

plot(cytokine_icc, cytokine_iqr)
