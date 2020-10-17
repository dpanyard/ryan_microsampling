library(metID)
library(devtools)
library(tidyverse)
library(metIdentify)
library(tinyTools)

#######
#combine peak table

##HILIC
##positive
setwd("/Users/brittanylee/Desktop/metID_Microsampling_247/HILIC/pos/")

rm(list = ls())

load("result.pHILIC")

library(metID)

identification.table1 <- getIdentificationTable(result.pHILIC[[1]],
                                                candidate.num = 1,
                                                type = "old")

identification.table1 <-
  identification.table1 %>%
  dplyr::filter(!is.na(Identification)) %>%
  dplyr::mutate(Level = 1)

dim(identification.table1)

identification.table2 <-
  getIdentificationTable(
    result.pHILIC[[2]],
    result.pHILIC[[3]],
    result.pHILIC[[4]],
    result.pHILIC[[5]],
    candidate.num = 1,
    type = "old"
  )

identification.table2 <-
  identification.table2 %>%
  dplyr::filter(!is.na(Identification), !is.element(name, identification.table1$name)) %>%
  dplyr::mutate(Level = 2)

dim(identification.table2)


identification.table3 <-
  getIdentificationTable2(result.pHILIC[[6]],
                          candidate.num = 3,
                          type = "old")


identification.table3 <-
  identification.table3 %>%
  dplyr::filter(!is.na(Identification), !is.element(name, identification.table1$name), !is.element(name, identification.table2$name)) %>%
  dplyr::mutate(Level = 3)

dim(identification.table3)


identification.table <-
  dplyr::full_join(rbind(identification.table1, identification.table2),
                   identification.table3,
                   by = colnames(identification.table3))

write.csv(identification.table, "identification.table_pHILIC.csv", row.names = FALSE)


identification.table.new <- trans2newStyle(identification.table = identification.table)
write.csv(identification.table.new, "identification.table.new_pHILIC.csv", row.names = FALSE)






##negative

rm(list = ls())

setwd("/Users/brittanylee/Desktop/metID_Microsampling_247/HILIC/neg")

load("result.nHILIC")

library(metID)

identification.table1 <- getIdentificationTable(result.nHILIC[[1]],
                                                candidate.num = 1,
                                                type = "old")

identification.table1 <-
  identification.table1 %>%
  dplyr::filter(!is.na(Identification)) %>%
  dplyr::mutate(Level = 1)

dim(identification.table1)

identification.table2 <-
  getIdentificationTable(
    result.nHILIC[[2]],
    result.nHILIC[[3]],
    result.nHILIC[[4]],
    result.nHILIC[[5]],
    candidate.num = 1,
    type = "old"
  )

identification.table2 <-
  identification.table2 %>%
  dplyr::filter(!is.na(Identification), !is.element(name, identification.table1$name)) %>%
  dplyr::mutate(Level = 2)

dim(identification.table2)


identification.table3 <-
  getIdentificationTable2(result.nHILIC[[6]],
                          candidate.num = 3,
                          type = "old")


identification.table3 <-
  identification.table3 %>%
  dplyr::filter(!is.na(Identification), !is.element(name, identification.table1$name), !is.element(name, identification.table2$name)) %>%
  dplyr::mutate(Level = 3)

dim(identification.table3)


identification.table <-
  dplyr::full_join(rbind(identification.table1, identification.table2),
                   identification.table3,
                   by = colnames(identification.table3))

write.csv(identification.table, "identification.table.csv", row.names = FALSE)


identification.table.new <- trans2newStyle(identification.table = identification.table)
write.csv(identification.table.new, "identification.table.new.csv", row.names = FALSE)

#combine peak table
####RPLC positive

rm(list=ls())

setwd("/Users/brittanylee/Desktop/metID_Microsampling_247/RPLC/pos")
load("result.pRPLC")

library(metID)

identification.table1 <- getIdentificationTable(result.pRPLC[[1]],
                                                candidate.num = 1,
                                                type = "old")
identification.table1 <-
  identification.table1 %>%
  dplyr::filter(!is.na(Identification)) %>%
  dplyr::mutate(Level = 1)

dim(identification.table1)

identification.table2 <-
  getIdentificationTable(
    result.pRPLC[[2]],
    result.pRPLC[[3]],
    result.pRPLC[[4]],
    result.pRPLC[[5]],
    candidate.num = 1,
    type = "old"
  )

identification.table2 <-
  identification.table2 %>%
  dplyr::filter(!is.na(Identification), !is.element(name, identification.table1$name)) %>%
  dplyr::mutate(Level = 2)

dim(identification.table2)


identification.table3 <-
  getIdentificationTable2(result.pRPLC[[6]],
                          candidate.num = 3,
                          type = "old")


identification.table3 <-
  identification.table3 %>%
  dplyr::filter(!is.na(Identification), !is.element(name, identification.table1$name), !is.element(name, identification.table2$name)) %>%
  dplyr::mutate(Level = 3)

dim(identification.table3)


identification.table <-
  dplyr::full_join(rbind(identification.table1, identification.table2),
                   identification.table3,
                   by = colnames(identification.table3))

write.csv(identification.table, "identification.table_pRPLC.csv", row.names = FALSE)


identification.table.new <- trans2newStyle(identification.table = identification.table)
write.csv(identification.table.new, "identification.table.new_pRPLC.csv", row.names = FALSE)


####RPLC negative
setwd("/Users/brittanylee/Desktop/metID_Microsampling_247/RPLC/neg")

rm(list = ls())

load("result.nRPLC")

identification.table1 <- getIdentificationTable(result.nRPLC[[1]],
                                                candidate.num = 1,
                                                type = "old")

identification.table1 <-
  identification.table1 %>%
  dplyr::filter(!is.na(Identification)) %>%
  dplyr::mutate(Level = 1)

dim(identification.table1)

identification.table2 <-
  getIdentificationTable(
    result.nRPLC[[2]],
    result.nRPLC[[3]],
    result.nRPLC[[4]],
    result.nRPLC[[5]],
    candidate.num = 1,
    type = "old"
  )

identification.table2 <-
  identification.table2 %>%
  dplyr::filter(!is.na(Identification), !is.element(name, identification.table1$name)) %>%
  dplyr::mutate(Level = 2)

dim(identification.table2)


identification.table3 <-
  getIdentificationTable2(result.nRPLC[[6]],
                          candidate.num = 3,
                          type = "old")

identification.table3 <-
  identification.table3 %>%
  dplyr::filter(!is.na(Identification), !is.element(name, identification.table1$name), !is.element(name, identification.table2$name)) %>%
  dplyr::mutate(Level = 3)

dim(identification.table3)


identification.table <-
  dplyr::full_join(rbind(identification.table1, identification.table2),
                   identification.table3,
                   by = colnames(identification.table3))

write.csv(identification.table, "identification.table_nRPLC.csv", row.names = FALSE)


identification.table.new <- trans2newStyle(identification.table = identification.table)
write.csv(identification.table.new, "identification.table.new_nRPLC.csv", row.names = FALSE)


