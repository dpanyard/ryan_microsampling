sxtTools::setwd_project()
setwd("data/")
library(tidyverse)
load("shake_omics.RData")

cytokine_data <- 
  omics %>% 
  dplyr::filter(Class == "Cytokine")

lipid_data <- 
  omics %>% 
  dplyr::filter(Class == "Lipid")

met_data <- 
  omics %>% 
  dplyr::filter(Class == "Metabolite")

##metabolomics data
  met_sample_info <- 
    met_data %>% 
    dplyr::ungroup() %>% 
    dplyr::select(sample_id = SampleID,
                  subject_id = PID, 
                  TP = TP) %>% 
    dplyr::distinct(.keep_all = TRUE)
  
  met_variable_info <- 
    met_data %>% 
    dplyr::ungroup() %>% 
    dplyr::select(mol_name = MolName,
                  pathway = Pathway,
                  subclass = Subclass,
                  p.value) %>% 
    dplyr::distinct(.keep_all = TRUE)

  met_expression_data <- 
    met_data %>% 
    dplyr::ungroup() %>% 
    dplyr::select(sample_id = SampleID,
                  mol_name = MolName,
                  intensity = Intensity
                  ) %>% 
    tidyr::pivot_wider(names_from = sample_id, 
                       values_from = intensity)
  
    

  
  
  
  
  
  
  