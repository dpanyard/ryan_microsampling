##
no_function()
sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())

###load data
load("data/shake_study/metabolome_data_analysis/data_preparation/expression_data")
load("data/shake_study/metabolome_data_analysis/data_preparation/sample_info")
load("data/shake_study/metabolome_data_analysis/data_preparation/variable_info")

sxtTools::setwd_project()
setwd("data/shake_study/metabolome_data_analysis/DEG")

dim(expression_data)
dim(sample_info)
dim(variable_info)

rownames(expression_data) == variable_info$variable_id

colnames(expression_data) == sample_info$sample_id

###find marker which are change according to pregnancy
expression_data[5,] %>% 
  as.numeric() %>% 
  density() %>%
  plot()

subject_data <- 
  expression_data %>% 
  dplyr::select(one_of(sample_info$sample_id))

##log transformation
subject_data <-
  log(subject_data + 1, 2)

### SAM analysi
# sam_test <-
#   samr::SAM(
#     x = subject_data,
#     y = sample_info$TP,
#     resp.type = "Quantitative",
#     geneid = rownames(subject_data),
#     genenames = rownames(subject_data),
#     return.x = FALSE,
#     fdr.output = 0.05,
#     regression.method = "ranks",
#     random.seed = "123",
#     nperms = 1000
#   )
# 
# save(sam_test, file = "sam_test")

load("sam_test")

metabolite_up <- sam_test$siggenes.table$genes.up

metabolite_down <- sam_test$siggenes.table$genes.lo

metabolite_marker <- c(metabolite_up[,1], metabolite_down[,1])

plot(sample_info$TP, subject_data[metabolite_marker[1],])

plot(sample_info$sample_id, expression_data[metabolite_down[7,1],])


metabolite_marker <-
  rbind(metabolite_up, metabolite_down) %>%
  as.data.frame() %>%
  dplyr::rename(score = `Score(d)`) %>%
  dplyr::mutate(score = as.numeric(score)) %>%
  dplyr::arrange(score) %>%
  dplyr::rename(variable_id = `Gene ID`,
                variable_name = `Gene Name`) %>%
  dplyr::mutate(variable_name = variable_id) %>%
  dplyr::left_join(variable_info, by = c("variable_id"))

save(metabolite_marker, file = "metabolite_marker")

load("metabolite_marker")
load("sam_test")



#####ga range
#to avoid source
no_exist_function()

sxtTools::setwd_project()
rm(list = ls())
library(tidyverse)

##load data
load("data/shake_study/metabolome_data_analysis/data_preparation/expression_data")
load("data/shake_study/metabolome_data_analysis/data_preparation/sample_info")
load("data/shake_study/metabolome_data_analysis/data_preparation/variable_info")

sxtTools::setwd_project()
setwd("data/shake_study/metabolome_data_analysis/DEG")

## for each person, just combine the samples in the same ga range
###combine different samples in one time together
library(plyr)

subject_data <- 
  expression_data %>% 
  dplyr::select(one_of(sample_info$sample_id))

##log transformation
subject_data <- 
  log(subject_data + 1, 2)

subject_data <- 
  apply(subject_data, 1, function(x){
    (x) / sd(x)
  })

subject_data2 <- 
  subject_data %>% 
  data.frame(., ga_range = sample_info$ga_range, stringsAsFactors = FALSE) %>% 
  mutate(ga_range = factor(ga_range,levels = sample_info$ga_range %>% unique() %>% stringr::str_sort(numeric = TRUE))) %>% 
  plyr::dlply(.variables = .(ga_range))

# subject_data_mean <-
#   lapply(subject_data2, function(x){
#     x <-
#       x %>%
#       dplyr::select(-ga_range)
#     apply(x, 2, mean)
#   }) %>%
#   do.call(cbind, . )
# 
# subject_data_sd <-
#   lapply(subject_data2, function(x){
#     x <-
#       x %>%
#       dplyr::select(-ga_range)
#     apply(x, 2, sd)
#   }) %>%
#   do.call(cbind, . )
# 
# 
# subject_data_sem <-
#   lapply(subject_data2, function(x){
#     x <-
#       x %>%
#       dplyr::select(-ga_range)
#     apply(x, 2, function(y){
#       sd(y)/sqrt(nrow(x))
#     })
#   }) %>%
#   do.call(cbind, . )
# 
# save(subject_data_mean, file = "subject_data_mean")
# save(subject_data_sd, file = "subject_data_sd")
# save(subject_data_sem, file = "subject_data_sem")

load("subject_data_mean")
load("subject_data_sd")
load("subject_data_sem")

subject_data2 <-
  lapply(subject_data2, function(x) {
    x <-
      x %>%
      dplyr::select(-ga_range)
  })

load("metabolite_marker")

plot(subject_data_mean[metabolite_marker$variable_id[2],])

# #find all the peaks in different time points
# fc_p_value <-
# pbapply::pblapply(subject_data2[-1], function(x){
#   p_value <- lapply(1:ncol(x), function(idx){
#     t.test(x[,idx], subject_data2[[1]][,idx])$p.value
#   }) %>%
#     unlist() %>%
#     p.adjust(method = "fdr")
# 
#   fc <- lapply(1:ncol(x), function(idx){
#     mean(x[,idx]) /mean(subject_data2[[1]][,idx])
#   }) %>%
#     unlist()
# 
#   fc[is.infinite(fc)] <- max(fc[!is.infinite(fc)])
# 
#   data.frame(p_value, fc, stringsAsFactors = FALSE)
# })
# 
# save(fc_p_value, file = "fc_p_value")

load("fc_p_value")

names(fc_p_value)

fc_p_value$`(14,16]`

dir.create("marker_in_different_points")

# for(idx in 1:length(fc_p_value)){
#   cat(idx, " ")
#   plot <-
#     volcano_plot(fc = fc_p_value[[idx]][,2],
#                  p_value = fc_p_value[[idx]][,1],
#                  p.cutoff = 0.05, fc.cutoff = 1,
#                  theme = "light")
# 
#   plot <-
#     plot +
#     labs(title = paste(names(fc_p_value)[idx], "(4,12]", sep = "/"))
# 
#   plot
# 
#   ggsave(
#     plot,
#     file = file.path(
#       "marker_in_different_points",
#       paste(names(fc_p_value)[idx], "(4,12]_light.pdf", sep = "_")
#     ),
#     width = 7,
#     height = 7,
#     bg = "transparent"
#   )
# 
#   ggsave(
#     plot,
#     file = file.path(
#       "marker_in_different_points",
#       paste(names(fc_p_value)[idx], "(4,12]_light.png", sep = "_")
#     ),
#     width = 7,
#     height = 7,
#     bg = "transparent"
#   )
# 
# 
#   plot <-
#     volcano_plot(fc = fc_p_value[[idx]][,2],
#                  p_value = fc_p_value[[idx]][,1],
#                  p.cutoff = 0.05, fc.cutoff = 1,
#                  theme = "dark")
# 
#   plot <-
#     plot +
#     labs(title = paste(names(fc_p_value)[idx], "(4,12]", sep = "/"))
# 
#   # plot
# 
#   ggsave(
#     plot,
#     file = file.path(
#       "marker_in_different_points",
#       paste(names(fc_p_value)[idx], "(4,12]_dark.pdf", sep = "_")
#     ),
#     width = 7,
#     height = 7,
#     bg = "transparent"
#   )
# 
# 
#   ggsave(
#     plot,
#     file = file.path(
#       "marker_in_different_points",
#       paste(names(fc_p_value)[idx], "(4,12]_dark.png", sep = "_")
#     ),
#     width = 7,
#     height = 7,
#     bg = "transparent"
#   )
# 
# }


##find markers for each time points
marker_each_point <- 
  lapply(fc_p_value, function(x){
    idx1 <- which(x$p_value < 0.05 & x$fc > 1)
    idx2 <- which(x$p_value < 0.05 & x$fc < 1)
    
    gene1 <- 
      try(
        data.frame(variable_id = variable_info$variable_id[idx1],
                   x[idx1,],
                   class = "increase",
                   stringsAsFactors = FALSE
        ),silent = TRUE 
      )
    
    if(class(gene1) == "try-error"){
      gene1 <- NULL
    }
    
    gene2 <- 
      try(
        data.frame(variable_id = variable_info$variable_id[idx2],
                   x[idx2,],
                   class = "decrease",
                   stringsAsFactors = FALSE
        ),silent = TRUE
      )
    
    if(class(gene2) == "try-error"){
      gene2 <- NULL
    }
    
    rbind(gene1, gene2)
  })

marker_each_point[[10]]

names(marker_each_point)

save(marker_each_point, file = "marker_each_point")

#####a sankey 
marker_each_point %>% 
  lapply(nrow) %>% 
  unlist()

all_marker_name <- 
  lapply(marker_each_point, function(x){
    x$variable_id
  }) %>% 
  unlist() %>% 
  unique()

length(all_marker_name)

library(ggalluvial)

temp_data <- 
  lapply(marker_each_point, function(x){
    if(is.null(x)){
      return(NULL)
    }
    x <- 
      data.frame(variable_id = all_marker_name,
                 stringsAsFactors = FALSE) %>% 
      left_join(x, by = "variable_id") %>% 
      dplyr::select(variable_id, class)
    
    x$class[is.na(x$class)] <- "no"
    x$freq <- 1
    x
    
  })

temp_data <-
  purrr::map2(.x = temp_data, .y = names(temp_data), .f = function(x,y){
    if(is.null(x)){
      return(NULL)
    }
    data.frame(x, point = y, stringsAsFactors = FALSE)
  })

temp_data <- 
  do.call(rbind, temp_data)

temp_data$point <- 
  factor(temp_data$point, levels = unique(temp_data$point))

RColorBrewer::display.brewer.all()

plot1 <- 
  ggplot(temp_data,
         aes(x = point, 
             y = freq,
             stratum = class, 
             alluvium = variable_id,
             fill = class, 
             label = class)) +
  scale_x_discrete(expand = c(.1, .1)) +
  ggalluvial::geom_flow() +
  labs(x = "", y = "") +
  scale_fill_manual(values = c(
    "increase" = ggsci::pal_aaas()(10)[2],
    "decrease" = ggsci::pal_aaas()(10)[1],
    "no" = "azure2"
  )) +
  ggalluvial::geom_stratum(alpha = 1, color = NA) +
  # geom_text(stat = "stratum", size = 3) +
  theme_bw() +
  theme(legend.position = "top", 
        panel.border = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 2),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
  )

plot1

# ggsave(
#   plot1,
#   file = file.path("marker_in_different_points", "gene_sankey_light.pdf"),
#   width = 14,
#   height = 7,
#   bg = "transparent"
# )

load("metabolite_marker")

overlap_data <-
  lapply(marker_each_point, function(x) {
    if (is.null(x)) {
      return(c(0, 0, 0))
    } else{
      up <-
        length(intersect(x$variable_id[x$class == "increase"],
                         metabolite_marker$variable_id[metabolite_marker$class == "up"])) * 100 / sum(x$class == "increase")
      down <-
        length(intersect(x$variable_id[x$class == "decrease"],
                         metabolite_marker$variable_id[metabolite_marker$class == "down"])) * 100 / sum(x$class == "decrease")
      total <-
        length(intersect(x$variable_id,
                         metabolite_marker$variable_id)) * 100 / nrow(x)
      c(up, down, total)
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::rename("up" = 'V1',
                "down" = 'V2',
                "total" = "V3")

plot3 <- 
  overlap_data %>% 
  tibble::rownames_to_column(var = "point") %>% 
  tidyr::pivot_longer(cols = -point, names_to = "class", values_to = "value") %>% 
  # dplyr::filter(value > 0) %>% 
  dplyr::mutate(class = factor(class, levels = c("total", "up", "down"))) %>% 
  ggplot(aes(point, class)) +
  geom_point(aes(fill = class, 
                 size = value), 
             shape = 21,
             show.legend = TRUE) +
  scale_size_continuous(range = c(1, 10)) +
  scale_fill_manual(values = c(
    "up" = ggsci::pal_aaas()(10)[2],
    "down" = ggsci::pal_aaas()(10)[1],
    "total" = "#FDB462"
  )) +
  guides(size = guide_legend(override.aes = list(color = "black"))) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_bw() +
  labs(x = "", y = "") +
  theme(legend.position = "bottom", 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 13), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  coord_fixed(ratio = 0.6)

plot3

# ggsave(
#   plot3,
#   file = file.path(
#     "marker_in_different_points",
#     "gene_marker_overlap_light.pdf"
#   ),
#   width = 14,
#   height = 7,
#   bg = "transparent"
# )





##if the DES are consistent 
increase_marker <- 
  lapply(marker_each_point, function(x){
    if(is.null(x)){
      return(NULL)
    }
    x$variable_id[x$class == "increase"]
  }) %>% 
  unlist() %>% 
  unique()


decrease_marker <- 
  lapply(marker_each_point, function(x){
    if(is.null(x)){
      return(NULL)
    }
    x$variable_id[x$class == "decrease"]
  }) %>% 
  unlist() %>% 
  unique()

library(VennDiagram)

up_marker <- 
  metabolite_marker$variable_id[metabolite_marker$class == "up"]

down_marker <- 
  metabolite_marker$variable_id[metabolite_marker$class == "down"]

# library(VennDiagram)
# 
# plot <- venn.diagram(x = list(sam = up_marker, each = increase_marker), filename = NULL)
# 
# grid.draw(plot)
# 
# plot <- venn.diagram(x = list(sam = down_marker, each = decrease_marker), filename = NULL)
# 
# grid.draw(plot)
# 
# example_metabolite1 <- 
#   metabolite_marker %>% 
#   dplyr::filter(class == "up",
#                 score == max(score, na.rm = TRUE)) %>% 
#   dplyr::pull(variable_id)
# 
# example_metabolite2 <- 
#   metabolite_marker %>% 
#   dplyr::filter(class == "down")
# # dplyr::pull(variable_id)
# 
# example_metabolite2 <- example_metabolite2$variable_id[1]
# 
# 
##log transformation
subject_data <- 
  expression_data %>% 
  dplyr::select(one_of(sample_info$sample_id))

subject_data <-
  log(subject_data + 1, 2)

subject_data <-
  apply(subject_data, 1, function(x){
    (x - mean(x)) / sd(x)
  })

# plot <- 
#   data.frame(ga = sample_info$g_stage, 
#              value = as.numeric(subject_data[,example_metabolite1])) %>% 
#   dplyr::mutate(
#     class = case_when(
#       ga == 50 ~ "After birth",
#       TRUE ~ "During pregnancy"
#     )
#   ) %>% 
#   ggplot(aes(ga, value)) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_point(shape = 21, aes(fill = class), size = 3) +
#   scale_color_manual(values = c("After birth" = "#FB8072",
#                                 "During pregnancy" = "black")) +
#   geom_smooth(method = "loess") +
#   labs(x = "Gestational age (GA, week)", y = "Scaled intensity") +
#   theme_bw() +
#   theme(axis.title = element_text(size = 13),
#         axis.text = element_text(size = 12), 
#         legend.title = element_text(size = 13),
#         legend.text = element_text(size = 12),
#         legend.position = c(0,1), 
#         legend.justification = c(0,1),
#         panel.background = element_rect(fill = "transparent", color = NA),
#         plot.background = element_rect(fill = "transparent", color = NA),
#         legend.background = element_rect(fill = "transparent", color = NA),
#         strip.background = element_rect(fill = "#0099B47F"),
#         strip.text = element_text(color = "white", size = 15),
#         panel.grid.major = element_blank()) 
# plot
# ggsave(plot, filename = "example_metabolite_up.pdf", height = 7, width = 7)
# 
# plot <- 
#   data.frame(ga = sample_info$g_stage, 
#              value = as.numeric(subject_data[,example_metabolite2])) %>% 
#   dplyr::mutate(
#     class = case_when(
#       ga == 50 ~ "After birth",
#       TRUE ~ "During pregnancy"
#     )
#   ) %>% 
#   ggplot(aes(ga, value)) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_point(shape = 21, aes(fill = class), size = 3) +
#   scale_color_manual(values = c("After birth" = "#FB8072",
#                                 "During pregnancy" = "black")) +
#   geom_smooth(method = "loess") +
#   labs(x = "Gestational age (GA, week)", y = "Scaled intensity") +
#   theme_bw() +
#   theme(axis.title = element_text(size = 13),
#         axis.text = element_text(size = 12), 
#         legend.title = element_text(size = 13),
#         legend.text = element_text(size = 12),
#         legend.position = c(0,1), 
#         legend.justification = c(0,1),
#         panel.background = element_rect(fill = "transparent", color = NA),
#         plot.background = element_rect(fill = "transparent", color = NA),
#         legend.background = element_rect(fill = "transparent", color = NA),
#         strip.background = element_rect(fill = "#0099B47F"),
#         strip.text = element_text(color = "white", size = 15),
#         panel.grid.major = element_blank()) 
# 
# plot
# 
# ggsave(plot, filename = "example_metabolite_down.pdf", height = 7, width = 7)


# ##get the very important metabolites
# consistent_up_info <- vector(mode = "list", length = length(marker_each_point) - 2)
# consistent_down_info <- vector(mode = "list", length = length(marker_each_point) - 2)
# 
# names(consistent_up_info) <- 
#   names(consistent_down_info) <- 
#   names(marker_each_point)[1:length(consistent_up_info)]
# 
# for(i in 1:length(consistent_up_info)){
#   cat(i, " ")
# 
#   if(is.null(marker_each_point[[i]])){
#     consistent_up_info[[i]] <- NA
#     consistent_down_info[[i]] <- NA
#     next()
#   }
# 
#   up_marker <-
#     marker_each_point[[i]] %>%
#     dplyr::filter(class == "increase") %>%
#     dplyr::pull(variable_id)
# 
#   down_marker <-
#     marker_each_point[[i]] %>%
#     dplyr::filter(class == "decrease") %>%
#     dplyr::pull(variable_id)
# 
#   if(i == 1){
#     previous_up_marker <- NULL
#     previous_down_marker <- NULL
#   }else{
#     previous_up_marker <-
#       marker_each_point[1:(i-1)] %>%
#       lapply(function(z){
#         if(is.null(z)){
#           return(NULL)
#         }else{
#           z %>%
#             dplyr::filter(class == "increase") %>%
#             dplyr::pull(variable_id)
#         }
#       }) %>%
#       unlist() %>%
#       unique()
# 
#     previous_down_marker <-
#       marker_each_point[1:(i-1)] %>%
#       lapply(function(z){
#         if(is.null(z)){
#           return(NULL)
#         }else{
#           z %>%
#             dplyr::filter(class == "decrease") %>%
#             dplyr::pull(variable_id)
#         }
#       }) %>%
#       unlist() %>%
#       unique()
#   }
# 
#   up_marker <-
#     setdiff(up_marker, previous_up_marker)
# 
#   down_marker <-
#     setdiff(down_marker, previous_down_marker)
# 
#   if(length(up_marker) == 0){
#     up_info <- NA
#   }else{
#     up_info <-
#       purrr::map(up_marker, function(x){
#         lapply(marker_each_point[-c(1:i)], function(y){
#           class <-
#             dplyr::filter(y, variable_id == x) %>%
#             dplyr::pull(class)
# 
#           if(length(class) == 0){
#             class <- "no"
#           }
#           class
#         }) %>%
#           unlist()
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
# 
#     rownames(up_info) <- up_marker
#   }
# 
#   if(length(down_marker) == 0){
#     down_info <- NA
#   }else{
#     down_info <-
#       purrr::map(down_marker, function(x){
#         lapply(marker_each_point[-c(1:i)], function(y){
#           class <-
#             dplyr::filter(y, variable_id == x) %>%
#             dplyr::pull(class)
# 
#           if(length(class) == 0){
#             class <- "no"
#           }
#           class
#         }) %>%
#           unlist()
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
# 
#     rownames(down_info) <- down_marker
#   }
#   consistent_up_info[[i]] <- up_info
#   consistent_down_info[[i]] <- down_info
# }
# 
# 
# save(consistent_up_info, file = "consistent_up_info")
# save(consistent_down_info, file = "consistent_down_info")

load("consistent_up_info")
load("consistent_down_info")

consistent_up_info2 <- 
  purrr::map2(.x = consistent_up_info, 
              .y = names(consistent_up_info), 
              .f = function(x, y){
                if(is.na(x)){
                  return(NULL)
                }
                
                number <- 
                  apply(x[,-1, drop = FALSE], 1, function(x){
                    sum(as.character(x) == "increase")
                  })
                
                x <- 
                  data.frame(x, number = number, 
                             check.names = FALSE,
                             stringsAsFactors = FALSE) %>% 
                  tibble::rownames_to_column(var = "variable_id") %>% 
                  dplyr::arrange(number) %>% 
                  tibble::column_to_rownames(var = "variable_id") %>% 
                  dplyr::select(-c(number))
                
                data.frame(
                  stage = y,
                  x,
                  stringsAsFactors = FALSE,
                  check.names = FALSE
                ) %>% 
                  tibble::rownames_to_column(var = "variable_id") %>%
                  tidyr::pivot_longer(cols = -c(variable_id, stage),
                                      names_to = "ga_range",
                                      values_to = "class")
              }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()


consistent_down_info2 <- 
  purrr::map2(.x = consistent_down_info, 
              .y = names(consistent_down_info), 
              .f = function(x, y){
                if(is.na(x)){
                  return(NULL)
                }
                
                number <- 
                  apply(x[,-1, drop = FALSE], 1, function(x){
                    sum(as.character(x) == "decrease")
                  })
                
                x <- 
                  data.frame(x, number = number, 
                             check.names = FALSE,
                             stringsAsFactors = FALSE) %>% 
                  tibble::rownames_to_column(var = "variable_id") %>% 
                  dplyr::arrange(number) %>% 
                  tibble::column_to_rownames(var = "variable_id") %>% 
                  dplyr::select(-c(number))
                
                data.frame(
                  stage = y,
                  x,
                  stringsAsFactors = FALSE,
                  check.names = FALSE
                ) %>% 
                  tibble::rownames_to_column(var = "variable_id") %>%
                  tidyr::pivot_longer(cols = -c(variable_id, stage),
                                      names_to = "ga_range",
                                      values_to = "class")
              }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()


plot_up1 <- 
  consistent_up_info2 %>% 
  dplyr::filter(!stage %in% c("(38,42]", "PP")) %>%
  dplyr::filter(!ga_range %in% c("PP")) %>% 
  dplyr::mutate(variable_id = factor(variable_id, 
                                     levels = unique(variable_id))) %>% 
  ggplot(aes(x = ga_range, y = variable_id)) +
  geom_tile(aes(fill = class), color = "white") +
  scale_fill_manual(values = c(
    "increase" = ggsci::pal_aaas()(n = 10)[2],
    "decrease" = ggsci::pal_aaas()(n = 10)[1],
    "no" = "white"
  )) +
  labs(y = "", x = "") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        plot.margin = margin(0,0,0,0),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey")
  )

plot_up1

consistent_up_info3 <-
  lapply(consistent_up_info, function(x) {
    if (is.na(x)) {
      return(x)
    } else{
      x <-
        x %>%
        dplyr::select(-PP)
      apply(x, 1, function(y) {
        sum(as.character(y) == "increase")
      }) %>%
        data.frame(count = .) %>%
        tibble::rownames_to_column(var = "variable_id")
    }
  })

consistent_up_info3 <-
  purrr::map2(
    .x = consistent_up_info3,
    .y = 12:1,
    .f = function(x, y) {
      if (is.na(x)) {
        return(x)
      }
      data.frame(x,
                 total = y,
                 stringsAsFactors = FALSE) %>%
        dplyr::mutate(per = count / total) %>%
        dplyr::arrange(per)
    }
  ) %>%
  do.call(rbind, .) %>%
  dplyr::filter(!is.na(per))

plot_up2 <-
  consistent_up_info3 %>%
  as.data.frame() %>%
  dplyr::mutate(variable_id = factor(variable_id, levels = variable_id)) %>%
  ggplot(aes(y = per, x = variable_id)) +
  geom_hline(yintercept = 2/3) +
  geom_line(aes(group = 1)) +
  # geom_point() +
  labs(y = "", x = "") +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    plot.margin = margin(0,0,0,0),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "grey")
  ) +
  coord_flip()

plot_up2  

library(patchwork)

plot_up <- 
  plot_up1 + plot_up2 + patchwork::plot_layout(nrow = 1, widths = c(3,1))

plot_up

# ggsave(plot_up, filename = "plot_up.pdf", width = 8, height = 7)

plot_down1 <-
  consistent_down_info2 %>%
  dplyr::filter(!stage %in% c("(38,42]", "PP")) %>%
  dplyr::filter(!ga_range %in% c("PP")) %>%
  dplyr::mutate(variable_id = factor(variable_id,
                                     levels = unique(variable_id))) %>%
  ggplot(aes(x = ga_range, y = variable_id)) +
  geom_tile(aes(fill = class)) +
  scale_fill_manual(
    values = c(
      "increase" = ggsci::pal_aaas()(n = 10)[2],
      "decrease" = ggsci::pal_aaas()(n = 10)[1],
      "no" = "white"
    )
  ) +
  labs(y = "", x = "") +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    plot.margin = margin(0, 0, 0, 0),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "grey")
  )

plot_down1

consistent_down_info3 <-
  lapply(consistent_down_info, function(x) {
    if (is.na(x)) {
      return(x)
    } else{
      x <-
        x %>%
        dplyr::select(-PP)
      apply(x, 1, function(y) {
        sum(as.character(y) == "decrease")
      }) %>%
        data.frame(count = .) %>%
        tibble::rownames_to_column(var = "variable_id")
    }
  })

consistent_down_info3 <-
  purrr::map2(
    .x = consistent_down_info3,
    .y = 12:1,
    .f = function(x, y) {
      if (is.na(x)) {
        return(x)
      }
      data.frame(x,
                 total = y,
                 stringsAsFactors = FALSE) %>%
        dplyr::mutate(per = count / total) %>%
        dplyr::arrange(per)
    }
  ) %>%
  do.call(rbind, .) %>%
  dplyr::filter(!is.na(per))


plot_down2 <-
  consistent_down_info3 %>%
  as.data.frame() %>%
  dplyr::mutate(variable_id = factor(variable_id, levels = variable_id)) %>%
  ggplot(aes(y = per, x = variable_id)) +
  geom_hline(yintercept = 2 / 3) +
  geom_line(aes(group = 1)) +
  # geom_point() +
  labs(y = "", x = "") +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    plot.margin = margin(0, 0, 0, 0),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "grey")
  ) +
  coord_flip()

plot_down2  

library(patchwork)

plot_down <- 
  plot_down1 + plot_down2 + patchwork::plot_layout(nrow = 1, widths = c(3,1))

plot_down

# ggsave(plot_down, filename = "plot_down.pdf", width = 8, height = 7)

consistent_up_info3$variable_id[consistent_up_info3$per > 2/3]

consistent_down_info3$variable_id[consistent_down_info3$per > 2/3]

##remove the metabolite which are not resume after birth
importance_up_metabolite <- 
  consistent_up_info3$variable_id[consistent_up_info3$per > 2/3]

importance_down_metabolite <- 
  consistent_down_info3$variable_id[consistent_down_info3$per > 2/3]

importance_up_metabolite <- 
  purrr::map(.x = importance_up_metabolite, .f = function(x){
    fc <- 
      mean(subject_data2$PP[,x])/mean(subject_data2$`(38,42]`[,x])
    p <- t.test(subject_data2$PP[,x], subject_data2$`(38,42]`[,x])$p.value
    c(fc, p)
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  dplyr::rename(fc = V1, p = V2) %>% 
  dplyr::mutate(fdr = p.adjust(p, method = "fdr")) %>% 
  data.frame(variable_id = importance_up_metabolite, ., stringsAsFactors = FALSE) %>% 
  dplyr::filter(fc < 1 & fdr < 0.05) %>% 
  dplyr::pull(variable_id)


importance_down_metabolite <- 
  purrr::map(.x = importance_down_metabolite, .f = function(x){
    fc <- 
      mean(subject_data2$PP[,x])/mean(subject_data2$`(38,42]`[,x])
    p <- t.test(subject_data2$PP[,x], subject_data2$`(38,42]`[,x])$p.value
    c(fc, p)
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  dplyr::rename(fc = V1, p = V2) %>% 
  dplyr::mutate(fdr = p.adjust(p, method = "fdr")) %>% 
  data.frame(variable_id = importance_down_metabolite, ., stringsAsFactors = FALSE) %>% 
  dplyr::filter(fc > 1 & fdr < 0.05) %>% 
  dplyr::pull(variable_id)

# library(VennDiagram)
# 
# plot <- 
#   venn.diagram(x = list(
#     sam = metabolite_marker$variable_id[metabolite_marker$class == "up"],
#     each = importance_up_metabolite
#   ), filename = NULL)
# 
# grid.draw(plot)
# 
# plot <- 
#   venn.diagram(x = list(
#     sam = metabolite_marker$variable_id[metabolite_marker$class == "down"],
#     each = importance_down_metabolite
#   ), filename = NULL)
# 
# grid.draw(plot)

importance_up_metabolite <- 
  intersect(
    metabolite_marker$variable_id[metabolite_marker$class == "up"],
    importance_up_metabolite
  )

importance_up_metabolite <-   
  dplyr::filter(metabolite_marker, variable_id %in% importance_up_metabolite)

importance_down_metabolite <- 
  intersect(
    metabolite_marker$variable_id[metabolite_marker$class == "down"],
    importance_down_metabolite
  )

importance_down_metabolite <-   
  dplyr::filter(metabolite_marker, variable_id %in% importance_down_metabolite)

# xlsx::write.xlsx(importance_up_metabolite,
#                  "importance_metabolites.xlsx",
#                  row.names = FALSE, sheetName = "UP",
#                  append = FALSE)
# xlsx::write.xlsx(importance_down_metabolite,
#                  "importance_metabolites.xlsx",
#                  row.names = FALSE, sheetName = "DOWN", append = TRUE)
# save(importance_up_metabolite, file = "importance_up_metabolite")
# save(importance_down_metabolite, file = "importance_down_metabolite")
load("importance_up_metabolite")
load("importance_down_metabolite")


####line plot for peak
temp_data <- subject_data_mean[importance_up_metabolite$variable_id,] %>% 
  apply(1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

plot <- 
  temp_data %>%
  tibble::rownames_to_column(var = "metabolite") %>%
  tidyr::pivot_longer(cols = -metabolite,
                      names_to = "ga",
                      values_to = "value") %>%
  dplyr::mutate(ga = factor(ga, levels = unique(ga))) %>%
  ggplot(aes(x  = ga, y = value)) +
  geom_line(aes(group = metabolite), color = ggsci::pal_aaas()(n = 10)[2]) +
  geom_boxplot(aes(x = ga, y = value), outlier.shape = NA) +
  geom_rect(
    aes(
      # xmin = "(38,42]",
      xmin = 14.5,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = ggsci::pal_aaas()(n = 10)[5],
    alpha = 0.5,
    data = data.frame(),
    inherit.aes = FALSE
  ) +
  theme_bw() +
  labs(x = "", y = "Z-score") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 12
    ),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.grid.minor = element_blank()
  )

plot

# ggsave(plot, filename = "important_up_peak_lineplot.pdf", width = 8, height = 7)

temp_data <- subject_data_mean[importance_down_metabolite$variable_id,] %>% 
  apply(1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

plot <- 
  temp_data %>%
  tibble::rownames_to_column(var = "metabolite") %>%
  tidyr::pivot_longer(cols = -metabolite,
                      names_to = "ga",
                      values_to = "value") %>%
  dplyr::mutate(ga = factor(ga, levels = unique(ga))) %>%
  ggplot(aes(x  = ga, y = value)) +
  geom_line(aes(group = metabolite), color = ggsci::pal_aaas()(n = 10)[1]) +
  geom_boxplot(aes(x = ga, y = value), outlier.shape = NA) +
  geom_rect(
    aes(
      # xmin = "(38,42]",
      xmin = 14.5,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = ggsci::pal_aaas()(n = 10)[5],
    alpha = 0.5,
    data = data.frame(),
    inherit.aes = FALSE
  ) +
  theme_bw() +
  labs(x = "", y = "Z-score") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 12
    ),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.grid.minor = element_blank()
  )

plot

# ggsave(plot, filename = "important_down_peak_lineplot.pdf", width = 8, height = 7)


importance_up_metabolite %>% 
  dplyr::filter(!is.na(Level)) %>% 
  dplyr::filter(Level == 1 | Level == 2) %>% 
  dplyr::pull(Compound.name)

importance_down_metabolite %>% 
  dplyr::filter(!is.na(Level)) %>% 
  dplyr::filter(Level == 1 | Level == 2) %>% 
  dplyr::pull(Compound.name)

##line plot for metabolite
temp_data <- subject_data_mean[importance_up_metabolite %>% 
                                 dplyr::filter(!is.na(Level)) %>% 
                                 dplyr::filter(Level == 1 | Level == 2) %>% 
                                 dplyr::pull(variable_id),] %>% 
  apply(1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

plot <- 
  temp_data %>%
  tibble::rownames_to_column(var = "metabolite") %>%
  tidyr::pivot_longer(cols = -metabolite,
                      names_to = "ga",
                      values_to = "value") %>%
  dplyr::mutate(ga = factor(ga, levels = unique(ga))) %>%
  ggplot(aes(x  = ga, y = value)) +
  geom_line(aes(group = metabolite), color = ggsci::pal_aaas()(n = 10)[2]) +
  geom_boxplot(aes(x = ga, y = value), outlier.shape = NA) +
  geom_rect(
    aes(
      # xmin = "(38,42]",
      xmin = 14.5,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = ggsci::pal_aaas()(n = 10)[5],
    alpha = 0.5,
    data = data.frame(),
    inherit.aes = FALSE
  ) +
  theme_bw() +
  labs(x = "", y = "Z-score") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 12
    ),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.grid.minor = element_blank()
  )

plot

# ggsave(plot, filename = "important_up_metabolite_lineplot.pdf", width = 8, height = 7)



temp_data <- subject_data_mean[importance_down_metabolite %>% 
                                 dplyr::filter(!is.na(Level)) %>% 
                                 dplyr::filter(Level == 1 | Level == 2) %>% 
                                 dplyr::pull(variable_id),] %>% 
  apply(1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

plot <- 
  temp_data %>%
  tibble::rownames_to_column(var = "metabolite") %>%
  tidyr::pivot_longer(cols = -metabolite,
                      names_to = "ga",
                      values_to = "value") %>%
  dplyr::mutate(ga = factor(ga, levels = unique(ga))) %>%
  ggplot(aes(x  = ga, y = value)) +
  geom_line(aes(group = metabolite), color = ggsci::pal_aaas()(n = 10)[1]) +
  geom_boxplot(aes(x = ga, y = value), outlier.shape = NA) +
  geom_rect(
    aes(
      # xmin = "(38,42]",
      xmin = 14.5,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = ggsci::pal_aaas()(n = 10)[5],
    alpha = 0.5,
    data = data.frame(),
    inherit.aes = FALSE
  ) +
  theme_bw() +
  labs(x = "", y = "Z-score") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 12
    ),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.grid.minor = element_blank()
  )

plot

# ggsave(plot, filename = "important_down_metabolite_lineplot.pdf", width = 8, height = 7)


###module analysis for each time points
# dir.create("PIUMet")
# piumet_file_up <-
#   importance_up_metabolite %>%
#   dplyr::select(name, mz, fdr) %>%
#   dplyr::mutate(fdr = -log(fdr, 10)) %>%
#   dplyr::mutate(polarity = case_when(
#     stringr::str_detect(name, "POS") ~ "positive",
#     stringr::str_detect(name, "NEG") ~ "negative"
#   )) %>%
#   dplyr::select(mz, polarity, fdr)
# 
# colnames(piumet_file_up) <- NULL
# 
# write.table(
#   piumet_file_up,
#   "PIUMet/piumet_file_up.txt",
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE
# )
# 
# 
# 
# piumet_file_down <-
#   importance_down_metabolite %>%
#   dplyr::select(name, mz, fdr) %>%
#   dplyr::mutate(fdr = -log(fdr, 10)) %>%
#   dplyr::mutate(polarity = case_when(
#     stringr::str_detect(name, "POS") ~ "positive",
#     stringr::str_detect(name, "NEG") ~ "negative"
#   )) %>%
#   dplyr::select(mz, polarity, fdr)
# 
# colnames(piumet_file_down) <- NULL
# 
# write.table(
#   piumet_file_down,
#   "PIUMet/piumet_file_down.txt",
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE
# )
# 
# 
# 
# ######up PIUMet
# readPIUMet(
#   path = "PIUMet/piumet_output_up/",
#   marker = importance_up_metabolite[, c("name", "mz", "rt")],
#   text = TRUE,
#   layout = "kk",
#   size_range = c(2, 8)
# )
# 
# 
# 
# load("PIUMet/piumet_output_up/Result/node_data")
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
#                       from = "Human Metabolome Database",
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
#   hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]
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
# load("PIUMet/piumet_output_up/Result/edge_data")
# 
# edge_data <-
# edge_data %>%
#   dplyr::filter(stringr::str_detect(from, "POS") |
#                   stringr::str_detect(from, "NEG") |
#                   stringr::str_detect(to, "POS") |
#                   stringr::str_detect(to, "NEG")
#                   ) %>%
#   dplyr::select(from, to) %>%
#   t() %>%
#   as.data.frame() %>%
#   purrr::map(.f = function(x){
#     if(stringr::str_detect(x[1], "POS") | stringr::str_detect(x[1], "NEG")){
#       return(x)
#     }else{
#       return(rev(x))
#     }
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame() %>%
#   dplyr::arrange(V1) %>%
#   dplyr::distinct()
# 
# colnames(edge_data) <- c('peak', "metabolite")
# 
# annotation_result <-
# annotation_result %>%
#   dplyr::filter(node_class == "Metabolite") %>%
#   dplyr::select(node, HMDB_ID, KEGG_ID, super.class)
# 
# annotation_result_up <-
# edge_data %>%
#   dplyr::left_join(annotation_result, by = c("metabolite" = "node"))
# 
# save(annotation_result_up,
#      file = "PIUMet/piumet_output_up/Result/annotation_result_up")
# 
# 
# 
# 
# ######down PIUMet
# readPIUMet(
#   path = "PIUMet/piumet_output_down/",
#   marker = importance_up_metabolite[, c("name", "mz", "rt")],
#   text = TRUE,
#   layout = "kk",
#   size_range = c(2, 8)
# )
# 
# 
# load("PIUMet/piumet_output_down/Result/node_data")
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
#                       from = "Human Metabolome Database",
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
#   hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]
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
# load("PIUMet/piumet_output_down/Result/edge_data")
# 
# edge_data <-
#   edge_data %>%
#   dplyr::filter(stringr::str_detect(from, "POS") |
#                   stringr::str_detect(from, "NEG") |
#                   stringr::str_detect(to, "POS") |
#                   stringr::str_detect(to, "NEG")
#   ) %>%
#   dplyr::select(from, to) %>%
#   t() %>%
#   as.data.frame() %>%
#   purrr::map(.f = function(x){
#     if(stringr::str_detect(x[1], "POS") | stringr::str_detect(x[1], "NEG")){
#       return(x)
#     }else{
#       return(rev(x))
#     }
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame() %>%
#   dplyr::arrange(V1) %>%
#   dplyr::distinct()
# 
# colnames(edge_data) <- c('peak', "metabolite")
# 
# annotation_result <-
#   annotation_result %>%
#   dplyr::filter(node_class == "Metabolite") %>%
#   dplyr::select(node, HMDB_ID, KEGG_ID, super.class)
# 
# annotation_result_down <-
#   edge_data %>%
#   dplyr::left_join(annotation_result, by = c("metabolite" = "node"))
# 
# save(annotation_result_down,
#      file = "PIUMet/piumet_output_down/Result/annotation_result_down")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ########functional annotation
# load("PIUMet/hsa_pathway")
# 
# sxtTools::setwd_project()
# source("R/2020_02_13/metabolome/pathway_enrichment.R")
# setwd("data_analysis_2020_02_13/metabolome_analysis/metabolome_DEG_analysis")
# 
# ##up cluster
# load("PIUMet/piumet_output_up/Result/annotation_result_up")
# 
# kegg_id <-
#   annotation_result_up$KEGG_ID
# 
# kegg_id <-
#   kegg_id[!is.na(kegg_id)]
# 
# 
# enrichment_kegg_up <-
#   enrich_pathway(id = kegg_id,
#                  pathway_database = hsa_pathway)
# 
# save(enrichment_kegg_up, file = "PIUMet/piumet_output_up/Result/enrichment_kegg_up")
# 
# 
# 
# 
# 
# 
# 
# 
# ##down cluster
# load("PIUMet/piumet_output_down/Result/annotation_result_down")
# 
# kegg_id <-
#   annotation_result_down$KEGG_ID
# 
# kegg_id <-
#   kegg_id[!is.na(kegg_id)]
# 
# 
# enrichment_kegg_down <-
#   enrich_pathway(id = kegg_id,
#                  pathway_database = hsa_pathway)
# 
# save(enrichment_kegg_down, file = "PIUMet/piumet_output_down/Result/enrichment_kegg_down")
# 
# 
# 
# 
# name <- 
#   stringr::str_split(names(hsa_pathway), pattern = ";") %>% 
#   lapply(function(x){
#     x[1]
#   }) %>% 
#   unlist()
# 
# temp_up <- 
# enrichment_kegg_up %>% 
#   dplyr::filter(p.value.fdr < 0.05) %>% 
#   dplyr::pull(Pathway.name) %>% 
#   purrr::map(.f = function(x){
#     id <- hsa_pathway[[match(x, name)]]
#     id <- annotation_result_up$metabolite[match(id, annotation_result_up$KEGG_ID)]
#     id[!is.na(id)]
#   })
# 
# names(temp_up) <- 
#   enrichment_kegg_up %>% 
#   dplyr::filter(p.value.fdr < 0.05) %>% 
#   dplyr::pull(Pathway.name)
# 
# 
# 
# 
# 
# readPIUMet(
#   path = "PIUMet/piumet_output_up/",
#   marker = importance_up_metabolite[, c("name", "mz", "rt")],
#   text = TRUE,
#   layout = "auto",
#   size_range = c(2, 8),
#   label.name = temp_up %>% unlist() %>% unique() %>% sort()
# )




###output all the metabolites with level 1 and level 2 annotation
####important up metabolites
temp_up_metabolite <- 
  importance_up_metabolite %>% 
  dplyr::filter(!is.na(Level)) %>% 
  dplyr::filter(Level == 1 | Level == 2)

dir.create("top100_up_metabolites")

subject_col <- colorRampPalette(colors = RColorBrewer::brewer.pal(11, name = "BrBG"))(n = 30)
names(subject_col) <- stringr::str_sort(unique(sample_info$subject_id), numeric = TRUE)

# for(idx in 1:nrow(temp_up_metabolite)){
#   cat(idx, " ")
# metabolite_id <- temp_up_metabolite$variable_id[idx]
# 
# score <- temp_up_metabolite$score[idx]
#   
# plot <-
#     data.frame(
#       subject_id = sample_info$subject_id,
#       ga = sample_info$g_stage,
#                value = as.numeric(subject_data[,temp_up_metabolite$variable_id[idx]])) %>%
#     dplyr::mutate(
#       class = case_when(
#         ga == 50 ~ "After birth",
#         TRUE ~ "During pregnancy"
#       )
#     ) %>%
#     ggplot(aes(ga, value)) +
#     geom_hline(yintercept = 0, linetype = 2) +
#     geom_point(shape = 21, aes(fill = subject_id), size = 3, alpha = 0.7) +
#     scale_fill_manual(values = subject_col) +
#     geom_smooth(aes(group = subject_id, color = subject_id),
#                 method = "loess",
#                 # color = "red",
#                 se = FALSE) +
#     guides(fill = guide_legend(title = "", nrow = 3),
#            color = guide_legend(title = "", nrow = 3)) +
#     scale_color_manual(values = subject_col) +
#     scale_x_continuous(breaks = c(10, 20, 30, 40, 50),
#                        labels = c(10, 20, 30, 40, "PP")) +
#     labs(x = "Gestational age (GA, week)", y = "Scaled intensity") +
#     theme_bw() +
#     theme(axis.title = element_text(size = 13),
#           axis.text = element_text(size = 12),
#           legend.title = element_text(size = 13),
#           legend.text = element_text(size = 12),
#           legend.position = "bottom",
#           # legend.justification = c(0,1),
#           panel.background = element_rect(fill = "transparent", color = NA),
#           plot.background = element_rect(fill = "transparent", color = NA),
#           legend.background = element_rect(fill = "transparent", color = NA),
#           strip.background = element_rect(fill = "#0099B47F"),
#           strip.text = element_text(color = "white", size = 15),
#           panel.grid.major = element_blank())
#   # plot
#   ggsave(
#     plot,
#     filename = file.path("top100_up_metabolites",
#                          paste(score, "_", metabolite_id, ".pdf", sep = "")),
#     height = 7.2,
#     width = 9.6
#   )
# 
#   ggsave(
#     plot,
#     filename = file.path("top100_up_metabolites",
#                          paste(score, "_", metabolite_id, ".png", sep = "")),
#     height = 7.2,
#     width = 9.6
#   )
# 
# }



dir.create("top100_up_metabolites_boxplot")

subject_col <- colorRampPalette(colors = RColorBrewer::brewer.pal(11, name = "BrBG"))(n = 30)
names(subject_col) <- stringr::str_sort(unique(sample_info$subject_id), numeric = TRUE)

# for(idx in 1:nrow(temp_up_metabolite)){
#   cat(idx, " ")
#   metabolite_id <- temp_up_metabolite$variable_id[idx]
#   
#   score <- temp_up_metabolite$score[idx]
#   
#   temp_data1 <- 
#     data.frame(
#       subject_id = sample_info$subject_id,
#       ga = sample_info$ga_range,
#       value = as.numeric(subject_data[,temp_up_metabolite$variable_id[idx]])
#       ) %>%
#     dplyr::mutate(
#       class = case_when(
#         ga == 50 ~ "After birth",
#         TRUE ~ "During pregnancy"
#       )
#     ) %>%
#     dplyr::mutate(ga = factor(ga, levels = sample_info$ga_range %>% 
#                                 unique() %>% 
#                                 stringr::str_sort(numeric = TRUE)))
#   
#   library(plyr)
#   temp_data2 <- 
#     temp_data1 %>%  
#     dlply(.variables = .(subject_id, ga)) %>% 
#     purrr::map(.f = function(x){
#       x$value <- mean(x$value)
#       x <- 
#         x %>% 
#         dplyr::distinct()
#       x
#     }) %>% 
#     do.call(rbind, .) %>% 
#     as.data.frame()
#     
#   
#   plot <-
#     temp_data1 %>% 
#     ggplot(aes(ga, value)) +
#     geom_hline(yintercept = 0, linetype = 2) +
#     geom_line(aes(x = ga, 
#                   y = value, 
#                   color = subject_id, 
#                   group = subject_id), 
#               data = temp_data2) +
#     geom_boxplot(aes(x = ga, y = value)) +
#     scale_fill_manual(values = subject_col) +
#     guides(fill = guide_legend(title = "", nrow = 3),
#            color = guide_legend(title = "", nrow = 3)) +
#     scale_color_manual(values = subject_col) +
#     labs(x = "Gestational age (GA, week)", y = "Scaled intensity") +
#     theme_bw() +
#     theme(axis.title = element_text(size = 13),
#           axis.text = element_text(size = 12),
#           legend.title = element_text(size = 13),
#           legend.text = element_text(size = 12),
#           legend.position = "bottom",
#           axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
#           # legend.justification = c(0,1),
#           panel.background = element_rect(fill = "transparent", color = NA),
#           plot.background = element_rect(fill = "transparent", color = NA),
#           legend.background = element_rect(fill = "transparent", color = NA),
#           strip.background = element_rect(fill = "#0099B47F"),
#           strip.text = element_text(color = "white", size = 15),
#           panel.grid.major = element_blank()
#           )
#   # plot
#   ggsave(
#     plot,
#     filename = file.path("top100_up_metabolites_boxplot",
#                          paste(score, "_", metabolite_id, ".pdf", sep = "")),
#     height = 7.2,
#     width = 9.6
#   )
#   
#   ggsave(
#     plot,
#     filename = file.path("top100_up_metabolites_boxplot",
#                          paste(score, "_", metabolite_id, ".png", sep = "")),
#     height = 7.2,
#     width = 9.6
#   )
#   
# }




####important down metabolites
temp_down_metabolite <-
  importance_down_metabolite %>%
  dplyr::filter(!is.na(Level)) %>%
  dplyr::filter(Level == 1 | Level == 2)
# 
# 
# dir.create("top100_down_metabolites")
# 
# subject_col <- colorRampPalette(colors = RColorBrewer::brewer.pal(11, name = "BrBG"))(n = 30)
# names(subject_col) <- stringr::str_sort(unique(sample_info$subject_id), numeric = TRUE)
# 
# for(idx in 1:nrow(temp_down_metabolite)){
#   cat(idx, " ")
#   metabolite_id <- temp_down_metabolite$variable_id[idx]
#   
#   score <- temp_down_metabolite$score[idx]
#   
#   plot <-
#     data.frame(
#       subject_id = sample_info$subject_id,
#       ga = sample_info$g_stage,
#       value = as.numeric(subject_data[,temp_down_metabolite$variable_id[idx]])) %>%
#     dplyr::mutate(
#       class = case_when(
#         ga == 50 ~ "After birth",
#         TRUE ~ "During pregnancy"
#       )
#     ) %>%
#     ggplot(aes(ga, value)) +
#     geom_hline(yintercept = 0, linetype = 2) +
#     geom_point(shape = 21, aes(fill = subject_id), size = 3, alpha = 0.7) +
#     scale_fill_manual(values = subject_col) +
#     geom_smooth(aes(group = subject_id, color = subject_id),
#                 method = "loess",
#                 # color = "red",
#                 se = FALSE) +
#     guides(fill = guide_legend(title = "", nrow = 3),
#            color = guide_legend(title = "", nrow = 3)) +
#     scale_color_manual(values = subject_col) +
#     scale_x_continuous(breaks = c(10, 20, 30, 40, 50),
#                        labels = c(10, 20, 30, 40, "PP")) +
#     labs(x = "Gestational age (GA, week)", y = "Scaled intensity") +
#     theme_bw() +
#     theme(axis.title = element_text(size = 13),
#           axis.text = element_text(size = 12),
#           legend.title = element_text(size = 13),
#           legend.text = element_text(size = 12),
#           legend.position = "bottom",
#           # legend.justification = c(0,1),
#           panel.background = element_rect(fill = "transparent", color = NA),
#           plot.background = element_rect(fill = "transparent", color = NA),
#           legend.background = element_rect(fill = "transparent", color = NA),
#           strip.background = element_rect(fill = "#0099B47F"),
#           strip.text = element_text(color = "white", size = 15),
#           panel.grid.major = element_blank())
#   # plot
#   ggsave(
#     plot,
#     filename = file.path("top100_down_metabolites",
#                          paste(score, "_", metabolite_id, ".pdf", sep = "")),
#     height = 7.2,
#     width = 9.6
#   )
#   
#   ggsave(
#     plot,
#     filename = file.path("top100_down_metabolites",
#                          paste(score, "_", metabolite_id, ".png", sep = "")),
#     height = 7.2,
#     width = 9.6
#   )
#   
# }
# 
# dir.create("top100_down_metabolites_boxplot")
# 
# subject_col <- colorRampPalette(colors = RColorBrewer::brewer.pal(11, name = "BrBG"))(n = 30)
# names(subject_col) <- stringr::str_sort(unique(sample_info$subject_id), numeric = TRUE)
# 
# for(idx in 1:nrow(temp_down_metabolite)){
#   cat(idx, " ")
#   metabolite_id <- temp_down_metabolite$variable_id[idx]
#   
#   score <- temp_down_metabolite$score[idx]
#   
#   temp_data1 <- 
#     data.frame(
#       subject_id = sample_info$subject_id,
#       ga = sample_info$ga_range,
#       value = as.numeric(subject_data[,temp_down_metabolite$variable_id[idx]])
#     ) %>%
#     dplyr::mutate(
#       class = case_when(
#         ga == 50 ~ "After birth",
#         TRUE ~ "During pregnancy"
#       )
#     ) %>%
#     dplyr::mutate(ga = factor(ga, levels = sample_info$ga_range %>% 
#                                 unique() %>% 
#                                 stringr::str_sort(numeric = TRUE)))
#   
#   library(plyr)
#   temp_data2 <- 
#     temp_data1 %>%  
#     dlply(.variables = .(subject_id, ga)) %>% 
#     purrr::map(.f = function(x){
#       x$value <- mean(x$value)
#       x <- 
#         x %>% 
#         dplyr::distinct()
#       x
#     }) %>% 
#     do.call(rbind, .) %>% 
#     as.data.frame()
#   
#   
#   plot <-
#     temp_data1 %>% 
#     ggplot(aes(ga, value)) +
#     geom_hline(yintercept = 0, linetype = 2) +
#     geom_line(aes(x = ga, 
#                   y = value, 
#                   color = subject_id, 
#                   group = subject_id), data = temp_data2) +
#     geom_boxplot(aes(x = ga, y = value)) +
#     scale_fill_manual(values = subject_col) +
#     guides(fill = guide_legend(title = "", nrow = 3),
#            color = guide_legend(title = "", nrow = 3)) +
#     scale_color_manual(values = subject_col) +
#     labs(x = "Gestational age (GA, week)", y = "Scaled intensity") +
#     theme_bw() +
#     theme(axis.title = element_text(size = 13),
#           axis.text = element_text(size = 12),
#           legend.title = element_text(size = 13),
#           legend.text = element_text(size = 12),
#           legend.position = "bottom",
#           axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
#           # legend.justification = c(0,1),
#           panel.background = element_rect(fill = "transparent", color = NA),
#           plot.background = element_rect(fill = "transparent", color = NA),
#           legend.background = element_rect(fill = "transparent", color = NA),
#           strip.background = element_rect(fill = "#0099B47F"),
#           strip.text = element_text(color = "white", size = 15),
#           panel.grid.major = element_blank()
#     )
#   # plot
#   ggsave(
#     plot,
#     filename = file.path("top100_down_metabolites_boxplot",
#                          paste(score, "_", metabolite_id, ".pdf", sep = "")),
#     height = 7.2,
#     width = 9.6
#   )
#   
#   ggsave(
#     plot,
#     filename = file.path("top100_down_metabolites_boxplot",
#                          paste(score, "_", metabolite_id, ".png", sep = "")),
#     height = 7.2,
#     width = 9.6
#   )
#   
# }





dir.create("top100_up_metabolites_whole")

subject_col <- colorRampPalette(colors = RColorBrewer::brewer.pal(11, name = "BrBG"))(n = 30)
names(subject_col) <- stringr::str_sort(unique(sample_info$subject_id), numeric = TRUE)

for(idx in 1:nrow(temp_up_metabolite)){
  cat(idx, " ")
  metabolite_id <- temp_up_metabolite$variable_id[idx]
  
  score <- temp_up_metabolite$score[idx]
  
  temp_data <-
    data.frame(
      ga = stringr::str_sort(unique(sample_info$ga_range), numeric = TRUE),
      mean = as.numeric(subject_data_mean[temp_up_metabolite$variable_id[idx],]),
      sd = as.numeric(subject_data_sd[temp_up_metabolite$variable_id[idx],]),
      sem = as.numeric(subject_data_sem[temp_up_metabolite$variable_id[idx],])
    ) %>%
    dplyr::mutate(ga = factor(ga, levels = ga))
  
  
  plot <-
    temp_data %>%
    ggplot(aes(ga, mean)) +
    # geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = ga, y = mean), shape = 21, 
               size = 3,
               fill = ggsci::pal_aaas()(n=10)[2]) +
    geom_errorbar(aes(x  = ga, 
                      y = mean, ymin = mean - sem, 
                      ymax = mean + sem), 
                  width= 0,
                  color = ggsci::pal_aaas()(n=10)[2]) +
    geom_line(aes(x = ga,
                  y = mean, group = 1),
              color = ggsci::pal_aaas()(n=10)[2]) +
    geom_rect(
      aes(
        xmin = 14.5,
        xmax = Inf,
        ymin = -Inf,
        ymax = Inf
      ),
      fill = ggsci::pal_aaas()(n = 10)[5],
      alpha = 0.5,
      data = data.frame(),
      inherit.aes = FALSE
    ) +
    labs(x = "Gestational age (GA, week)", y = "Scaled intensity") +
    theme_bw() +
    theme(axis.title = element_text(size = 13),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
          # legend.justification = c(0,1),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          legend.background = element_rect(fill = "transparent", color = NA),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15),
          panel.grid.minor = element_blank()
    )
  # plot
  ggsave(
    plot,
    filename = file.path("top100_up_metabolites_whole",
                         paste(score, "_", metabolite_id, ".pdf", sep = "")),
    height = 7.2,
    width = 9.6
  )
  
  ggsave(
    plot,
    filename = file.path("top100_up_metabolites_whole",
                         paste(score, "_", metabolite_id, ".png", sep = "")),
    height = 7.2,
    width = 9.6
  )
  
}








dir.create("top100_down_metabolites_whole")

subject_col <- colorRampPalette(colors = RColorBrewer::brewer.pal(11, name = "BrBG"))(n = 30)
names(subject_col) <- stringr::str_sort(unique(sample_info$subject_id), numeric = TRUE)

for(idx in 1:nrow(temp_down_metabolite)){
  cat(idx, " ")
  metabolite_id <- temp_down_metabolite$variable_id[idx]
  
  score <- temp_down_metabolite$score[idx]
  
  temp_data <-
    data.frame(
      ga = stringr::str_sort(unique(sample_info$ga_range), numeric = TRUE),
      mean = as.numeric(subject_data_mean[temp_down_metabolite$variable_id[idx],]),
      sd = as.numeric(subject_data_sd[temp_down_metabolite$variable_id[idx],]),
      sem = as.numeric(subject_data_sem[temp_down_metabolite$variable_id[idx],])
    ) %>%
    dplyr::mutate(ga = factor(ga, levels = ga))
  
  
  plot <-
    temp_data %>%
    ggplot(aes(ga, mean)) +
    # geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = ga, y = mean), shape = 21, 
               size = 3,
               fill = ggsci::pal_aaas()(n=10)[1]) +
    geom_errorbar(aes(x  = ga, 
                      y = mean, ymin = mean - sem, 
                      ymax = mean + sem), 
                  width= 0,
                  color = ggsci::pal_aaas()(n=10)[1]) +
    geom_line(aes(x = ga,
                  y = mean, group = 1),
              color = ggsci::pal_aaas()(n=10)[1]) +
    geom_rect(
      aes(
        xmin = 14.5,
        xmax = Inf,
        ymin = -Inf,
        ymax = Inf
      ),
      fill = ggsci::pal_aaas()(n = 10)[5],
      alpha = 0.5,
      data = data.frame(),
      inherit.aes = FALSE
    ) +
    labs(x = "Gestational age (GA, week)", y = "Scaled intensity") +
    theme_bw() +
    theme(axis.title = element_text(size = 13),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
          # legend.justification = c(0,1),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          legend.background = element_rect(fill = "transparent", color = NA),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15),
          panel.grid.minor = element_blank()
    )
  # plot
  ggsave(
    plot,
    filename = file.path("top100_down_metabolites_whole",
                         paste(score, "_", metabolite_id, ".pdf", sep = "")),
    height = 7.2,
    width = 9.6
  )
  
  ggsave(
    plot,
    filename = file.path("top100_down_metabolites_whole",
                         paste(score, "_", metabolite_id, ".png", sep = "")),
    height = 7.2,
    width = 9.6
  )
  
}









