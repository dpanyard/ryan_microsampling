##
no_function()

sxtTools::setwd_project()
setwd("data/7_24_mike/")

library(tidyverse)
rm(list = ls())

load("MOmics_01.RData")

time_anno = read.csv("TimeAnno.csv")
time_anno = time_anno %>% select(-time, -date) %>% rename(DT = date_time)
events = gather(time_anno,"Type","Details",5:6)
events = events[!is.na(events$Details),]

library(lubridate)

start = as_datetime("2019-04-29 03:00:00", tz="America/Los_Angeles")
events$DT = as.POSIXct(events$DT, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")
events = events %>% mutate(day = date(DT), hour_of_day=hour(DT), tod = hms::as.hms(DT, tz="America/Los_Angeles"))



####lipidomics
lipids = DF247_lipids2[,c(1,2,4,5,10,13)]
lipids$SampleID = as.numeric(gsub("24_7_prepsample","",lipids$SampleID))
lipids = lipids[!is.na(lipids$SampleID),]
lipids$MolClass = "Lipid"
lipids = rename(lipids, MolName = LipidSpecies, MolSubclass = LipidClass, Intensity = log_sample_nmol_per_g_concentration, DT = CollectionTime_format)
lipids = select(lipids, -collector)
lipids$SampleID = as.character(lipids$SampleID)


sample_info <-
  lipids %>%
  dplyr::ungroup() %>%
  dplyr::select(sample_id = SampleID,
                # subject_id = PID,
                DT = DT) %>%
  dplyr::distinct(.keep_all = TRUE) %>% 
  as.data.frame()

variable_info <-
  lipids %>%
  dplyr::ungroup() %>%
  dplyr::select(mol_name = MolName,
                mol_class = MolClass,
                mol_subclass = MolSubclass) %>%
  dplyr::distinct(.keep_all = TRUE) %>% 
  as.data.frame() %>% 
  dplyr::mutate(mol_name = as.character(mol_name))

expression_data <-
  lipids %>%
  dplyr::ungroup() %>%
  dplyr::select(sample_id = SampleID,
                mol_name = MolName,
                intensity = Intensity) %>%
  tidyr::pivot_wider(names_from = sample_id,
                     values_from = intensity) %>% 
  as.data.frame()

expression_data$mol_name == variable_info$mol_name

variable_info <- 
  variable_info %>% 
  dplyr::mutate(variable_id = paste("lipid", 1:nrow(variable_info), sep = "_")) %>% 
  dplyr::select(variable_id, everything())

rownames(expression_data) <- variable_info$variable_id

expression_data <- 
  expression_data %>% 
  dplyr::select(-mol_name) %>% 
  as.data.frame()

dim(expression_data)
dim(variable_info)
dim(sample_info)

colnames(expression_data) == sample_info$sample_id

variable_info

save(sample_info, file = "lipidomics/data_preparation/sample_info")
save(variable_info, file = "lipidomics/data_preparation/variable_info")
save(expression_data, file = "lipidomics/data_preparation/expression_data")



###### proteins - plasma
library(data.table)
proteins <- data.table(fread("allproteins_overtime_annotated_Mike247plasma.csv"), sep='\t')
uniprot = read.csv("uniprot_match.csv")
uniprot$Entry_name = gsub("_HUMAN","",uniprot$Entry_name)
proteinsL = gather(proteins[,c(1:313,321,325)], Sample, Intensity, 1:313)
proteinsL = left_join(proteinsL,select(uniprot,Entry,Entry_name), by=c("Sample"="Entry"))
proteinsL = subset(proteinsL,!is.na(Entry_name))
proteinsL = unite(data = proteinsL, col = "MolName", Sample, Entry_name)
proteinsL = proteinsL %>% rename(DT = date_time, SampleID = SampleIndex)
proteinsL$MolClass = "Protein"
proteinsL$MolSubclass = ""
proteinsL$SampleID = as.character(proteinsL$SampleID)
proteinsL$DT = as.POSIXct(proteinsL$DT, tz = "America/Los_Angeles", format = "%Y-%m-%d %H:%M:%S")

sample_info <-
  proteinsL %>%
  dplyr::ungroup() %>%
  dplyr::select(sample_id = SampleID,
                # subject_id = PID,
                DT = DT) %>%
  dplyr::distinct(.keep_all = TRUE) %>% 
  as.data.frame()

variable_info <-
  proteinsL %>%
  dplyr::ungroup() %>%
  dplyr::select(mol_name = MolName,
                mol_class = MolClass,
                mol_subclass = MolSubclass) %>%
  dplyr::distinct(.keep_all = TRUE) %>% 
  as.data.frame() %>% 
  dplyr::mutate(mol_name = as.character(mol_name))

expression_data <-
  proteinsL %>%
  dplyr::ungroup() %>%
  dplyr::select(sample_id = SampleID,
                mol_name = MolName,
                intensity = Intensity) %>%
  tidyr::pivot_wider(names_from = sample_id,
                     values_from = intensity) %>% 
  as.data.frame()

expression_data$mol_name == variable_info$mol_name

variable_info <- 
  variable_info %>% 
  dplyr::mutate(variable_id = paste("protein", 1:nrow(variable_info), sep = "_")) %>% 
  dplyr::select(variable_id, everything())

rownames(expression_data) <- variable_info$variable_id

expression_data <- 
  expression_data %>% 
  dplyr::select(-mol_name) %>% 
  as.data.frame()

dim(expression_data)
dim(variable_info)
dim(sample_info)

colnames(expression_data) == sample_info$sample_id

variable_info <-
  variable_info %>% 
  dplyr::mutate(protein_name2 = 
                  stringr::str_split(mol_name, "_", n = 2) %>% 
                  purrr::map(.f = function(x)x[2]) %>% 
                  unlist())

variable_info <- 
variable_info %>% 
  dplyr::left_join(uniprot, by = c("protein_name2" = "Entry_name"))


save(sample_info, file = "proteomics/data_preparation/sample_info")
save(variable_info, file = "proteomics/data_preparation/variable_info")
save(expression_data, file = "proteomics/data_preparation/expression_data")











############## metabolomics----------------------------------------------------
metab_all = read.csv("raw_data_from_box/Microsampling_247_combined_all.csv")

### metID
pR = read.csv("metabolomics/metID_Microsampling_247/RPLC/pos/identification.table.new_pRPLC.csv")
nR = read.csv("metabolomics/metID_Microsampling_247/RPLC/neg/identification.table.new_nRPLC.csv")
pH = read.csv("metabolomics/metID_Microsampling_247/HILIC/pos/identification.table.new_pHILIC.csv")
nH = read.csv("metabolomics/metID_Microsampling_247/HILIC/neg/identification.table.new_nHILIC.csv")

pR$mode="pR"
nR$mode="nR"
pH$mode="pH"
nH$mode="nH"
all_metID = rbind(pR,nR,pH,nH)

mets_join = left_join(metab_all,all_metID,by=c("Compounds_ID"="name"))
write.table(mets_join,file="metabolomics/24-7_joined_110820b.csv",sep=",",row.names=F)

mets_join %>% ungroup %>% group_by(mode) %>% summarise(count = n())

metab_all = subset(metab_all, !Metabolite_val %in% 
                     c("0", "", "NA")) %>% 
  select(Metabolite_val:CAS_val, Super.pathway:Sub.pathway, X1:X106)

j = inner_join(
  metab_all,
  select(all_metID, name, Compound.name:Database),
  by = c("Compounds_ID" = "name")
)

j=select(j,Compound.name:Database,everything())

j=j[!duplicated(j$Compounds_ID),]

write.table(j,file="24-7_joined.csv",sep=",",row.names=F)

inj_order = read_excel("/Users/ryankellogg/Box/Microsampling/metabolites/MS_247/Data/metabolomics/M-Sampling-Brittany.xlsx","247 ome")
inj_order = rename(inj_order, SampleID = "PrepIndex (SampleName)")
inj_order$AcqOrder = paste0("X", inj_order$AcqOrder)

metab = gather(metab_all,AcqOrder,Intensity,-(Metabolite_val:Sub.pathway))
metab = rename(metab,MolName = Metabolite_val, HMDB_ID = HMDB_val)
metab = left_join(metab,select(inj_order,AcqOrder,SampleID)) %>% select(-AcqOrder)

# median norm
metab = metab %>% group_by(SampleID) %>% mutate(Intensity = log10(Intensity/median(Intensity)))
#metab$SampleID = as.numeric(gsub("X24_7_prepsample","",metab$SampleID))
metab = metab[!is.na(metab$SampleID),]
metab$SampleID = as.character(metab$SampleID)

anno = read.csv("/Users/ryankellogg/Box/Microsampling/Analysis/annotations/24_7_Stability_MasterIDReferenz_DH4_KE2_rk.csv")
metab = left_join(metab,anno[,c("PrepID","CollectionTime", "collector")], by=c("SampleID" = "PrepID"))
metab$DT = as.POSIXct(metab$CollectionTime, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")
metab = select(metab, -collector, -CollectionTime)
metab$MolClass = "Metabolite"
metab$MolSubclass = NA
metab = subset(metab, !is.na(DT))
#metab$MolName = make.unique(as.character(metab$MolName))
metab = metab %>% group_by(SampleID) %>% mutate(MolName = make.unique(as.character(MolName)))

# weird sample #50
metab = metab[metab$SampleID != 50,]
metab=subset(metab, !(DT=="2019-05-05 09:11:00" & metab$MolName=="Caffeine"))



# metabolic protein panel
MetaData=read.csv("/Users/ryankellogg/Box/Microsampling/Analysis/annotations/Copy of 24_7_Stability_MasterIDReferenz_DH4_KE2.csv") 

MetaData = MetaData[MetaData$Study == "F0" & MetaData$Original_SampleID %in% 1:97,]
MetaData$Name = paste("DH-",MetaData$Original_SampleID,sep="")

MP = read.csv("/Users/ryankellogg/Box/Microsampling/cytokines/24-7 omics-MetabolicHormone-Ryan Metabolic Plate-1.csv")
MP = MP[2:97,]

MP = left_join(MP,MetaData[,c("Name","CollectionTime","PrepIndex")])
MP$DT = as.POSIXct(MP$CollectionTime, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")
MPL = gather(MP,"Cytokine","MFI",5:21)
MPL$Intensity = log10(as.numeric(MPL$MFI))
MPL = rename(MPL, MolName = Cytokine, SampleID=PrepIndex)
MPL$MolClass = "MetabolicPanel"
MPL$MolSubclass = NA
MPL = select(MPL, SampleID, DT, MolName, Intensity, MolClass, MolSubclass)
MPL[MPL$MolName=="TNFA",]$MolName = "TNFA_MP"
MPL[MPL$MolName=="MCP1",]$MolName = "MCP1_MP"
MPL[MPL$MolName=="CHEX1",]$MolName = "CHEX1_MP"
MPL[MPL$MolName=="CHEX2",]$MolName = "CHEX2_MP"
MPL[MPL$MolName=="CHEX3",]$MolName = "CHEX3_MP"
MPL[MPL$MolName=="CHEX4",]$MolName = "CHEX4_MP"
MPL[MPL$MolName=="IL6",]$MolName = "IL6_MP"




# cortisol 
MetaData=read.csv("raw_data_from_box/Copy of 24_7_Stability_MasterIDReferenz_DH4_KE2.csv") 

MetaData = MetaData[MetaData$Study == "F0" & MetaData$Original_SampleID %in% 1:97,]
MetaData$Name = paste("DH-",MetaData$Original_SampleID,sep="")

cort = read.csv("raw_data_from_box/RYAN CORTISOL 8-28-19.csv")
cort = cort[10:102,]

cort = left_join(cort,MetaData[,c("Name","CollectionTime","PrepIndex")])
cort$DT = as.POSIXct(cort$CollectionTime, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")
cort = rename(cort, SampleID=PrepIndex)

cortL = gather(cort,"Cytokine","MFI",5)
cortL$Intensity = log10(as.numeric(cortL$MFI))
cortL$MolClass = "CortisolEnzymatic"
cortL$MolSubclass = NA
cortL$MolName = "Cortisol"
cortL = select(cortL, SampleID, DT, MolName, Intensity, MolClass, MolSubclass)

#41-plex
cyto = read.csv("raw_data_from_box/24-7 omics-HumanLuminexMAG42plex-Mitra H41 Plex-1.csv")
cyto = cyto[10:105,]

cyto = left_join(cyto,MetaData[,c("Name","CollectionTime","PrepIndex")])
cyto$DT = as.POSIXct(cyto$CollectionTime, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")

cytoL = gather(cyto,"Cytokine","MFI",5:49)
cytoL$Intensity = log10(as.numeric(cytoL$MFI))
cytoL$MolClass = "Cytokine_41Plex"
cytoL$MolSubclass = NA
cytoL = rename(cytoL, SampleID = PrepIndex, MolName = Cytokine)
cytoL = select(cytoL, SampleID, DT, MolName, Intensity, MolClass, MolSubclass)

all_omes = bind_rows(metab, lipids, proteinsL, MPL, cortL, cytoL)

# add MolNum
temp = tibble(MolName = unique(all_omes$MolName)) #, MolIndex = 1:n(MolName))
temp$MolNum = 1:nrow(temp)
all_omes = left_join(all_omes, temp)

all_omes_wear = bind_rows(all_omes, wearables)
all_omes_wear_food = bind_rows(all_omes, wearables, FoodTimeL_sum)

all_omes = all_omes %>% mutate(day = date(DT), hour_of_day=hour(DT), tod = hms::as.hms(DT, tz="America/Los_Angeles"))
all_omes = all_omes[all_omes$DT < as_date("2019-05-07"),]

all_omes_wear = all_omes_wear %>% mutate(day = date(DT), hour_of_day=hour(DT), tod = hms::as.hms(DT, tz="America/Los_Angeles"))
all_omes_wear_food = all_omes_wear_food %>% mutate(day = date(DT), hour_of_day=hour(DT), tod = hms::as.hms(DT, tz="America/Los_Angeles"))


grid = seq(as.POSIXct(floor_date(min(all_omes$DT),"hour")),as.POSIXct(max(all_omes$DT)),by="1 hour")

all_omes_hr = all_omes %>% select(-KEGG_val, -HMDB_ID, -CAS_val, -Super.pathway, -Sub.pathway) %>% group_by(MolName,MolClass,MolSubclass,MolNum) %>% mutate(MI = mean(Intensity), SDI = sd(Intensity), Intensity_z = (Intensity - mean(Intensity,na.rm=T))/sd(Intensity,na.rm=T)) %>% nest() %>%
  mutate(hr = map(data, function(x) as_tibble(approx(x$DT, x$Intensity_z, xout=grid, rule=1)))) %>%
  unnest(hr) %>% rename(DT = x, Intensity = y)

all_omes_hr = all_omes_hr[!is.na(all_omes_hr$Intensity),]
all_omes_hr = all_omes_hr %>% group_by(MolName) %>% mutate(Hr = as.factor(as.numeric(difftime(DT, first(DT)), units="hours")), Hr_day = hour(DT), AM = am(DT))
all_omes_hr = all_omes_hr[all_omes_hr$DT < as_date("2019-05-07"),]
all_omes_hr = all_omes_hr %>% mutate(day = date(DT), hour_of_day=hour(DT), tod = hms::as.hms(DT, tz="America/Los_Angeles"))


all_omes_wear = all_omes_wear %>% select(-KEGG_val, -HMDB_ID, -CAS_val, -Super.pathway, -Sub.pathway) %>% group_by(MolName,MolClass,MolSubclass,MolNum) %>% mutate(MI = mean(Intensity), SDI = sd(Intensity), Intensity_z = (Intensity - mean(Intensity,na.rm=T))/sd(Intensity,na.rm=T)) %>% nest() 
all_omes_wear = all_omes_wear %>%
  mutate(hr = map(data, function(x) as_tibble(approx(x$DT, x$Intensity_z, xout=grid, rule=1)))) %>%
  unnest(hr) %>% rename(DT = x, Intensity = y)
all_omes_wear = all_omes_wear[!is.na(all_omes_wear$Intensity),]
all_omes_wear = all_omes_wear %>% group_by(MolName) %>% mutate(Hr = as.factor(as.numeric(difftime(DT, first(DT)), units="hours")), Hr_day = hour(DT), AM = am(DT))
all_omes_wear = all_omes_wear[all_omes_wear$DT > "2019-04-29 03:00:00" & all_omes_wear$DT < as_date("2019-05-07"),]
all_omes_wear = all_omes_wear %>% mutate(day = date(DT), hour_of_day=hour(DT), tod = hms::as.hms(DT, tz="America/Los_Angeles"))







