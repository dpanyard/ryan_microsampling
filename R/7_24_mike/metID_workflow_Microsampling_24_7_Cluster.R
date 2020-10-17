## Type this into the cluster before starting so the run will continue even if you don't have the session open on your computer
## tmux new -s test_name

library(metID)
library(devtools)
library(tidyverse)
library(metIdentify)
library(tinyTools)

####HILIC
##positive
setwd("/home/leeba/metID_Microsampling_247/HILIC/pos")

library(metID)

parameter1 <- metIdentifyParam(
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 60,
  ms1.match.ppm = 25,
  ms2.match.tol = 0.4,
  rt.match.tol = 90,
  polarity = "positive",
  ce = "all",
  column = "hilic",
  total.score.tol = 0,
  candidate.num = 3,
  database = "msDatabase_hilic0.0.1",
  threads = 10
)

parameter2 <- metIdentifyParam(
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 60,
  ms1.match.ppm = 25,
  ms2.match.tol = 0.4,
  rt.match.tol = 90,
  polarity = "positive",
  ce = "all",
  column = "hilic",
  total.score.tol = 0,
  candidate.num = 3,
  database = "hmdbDatabase0.0.1",
  threads = 10
)

parameter3 <- metIdentifyParam(
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 60,
  ms1.match.ppm = 25,
  ms2.match.tol = 0.4,
  rt.match.tol = 90,
  polarity = "positive",
  ce = "all",
  column = "hilic",
  total.score.tol = 0,
  candidate.num = 3,
  database = "massbankDatabase0.0.1",
  threads = 10
)


parameter4 <- metIdentifyParam(
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 60,
  ms1.match.ppm = 25,
  ms2.match.tol = 0.4,
  rt.match.tol = 90,
  polarity = "positive",
  ce = "all",
  column = "hilic",
  total.score.tol = 0,
  candidate.num = 3,
  database = "monaDatabase0.0.1",
  threads = 10
)



parameter5 <- metIdentifyParam(
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 60,
  ms1.match.ppm = 25,
  ms2.match.tol = 0.4,
  rt.match.tol = 90,
  polarity = "positive",
  ce = "all",
  column = "hilic",
  total.score.tol = 0,
  candidate.num = 3,
  database = "orbitrapDatabase0.0.1",
  threads = 10
)

parameter6 <- mzIdentifyParam(
  ms1.match.ppm = 25,
  polarity = "positive",
  column = "hilic",
  candidate.num = 3,
  database = "HMDB.metabolite.data",
  threads = 10
)
ms1.data <- grep("\\.csv", dir(), value = TRUE)
ms2.data <- grep("mgf", dir(), value = TRUE)
result.pHILIC <- metIdentify4all(
  ms1.data = ms1.data,
  ms2.data = ms2.data,
  parameter.list = c(
    parameter1,
    parameter2,
    parameter3,
    parameter4,
    parameter5,
    parameter6
  ),
  path = "."
)

save(result.pHILIC, file = "result.pHILIC")
print(1)


##negative
setwd("/home/leeba/metID_Microsampling_247/HILIC/neg")

library(metID)

parameter1 <- metIdentifyParam(
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 60,
  ms1.match.ppm = 25,
  ms2.match.tol = 0.4,
  rt.match.tol = 90,
  polarity = "negative",
  ce = "all",
  column = "hilic",
  total.score.tol = 0,
  candidate.num = 3,
  database = "msDatabase_hilic0.0.1",
  threads = 10
)

parameter2 <- metIdentifyParam(
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 60,
  ms1.match.ppm = 25,
  ms2.match.tol = 0.4,
  rt.match.tol = 90,
  polarity = "negative",
  ce = "all",
  column = "hilic",
  total.score.tol = 0,
  candidate.num = 3,
  database = "hmdbDatabase0.0.1",
  threads = 10
)

parameter3 <- metIdentifyParam(
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 60,
  ms1.match.ppm = 25,
  ms2.match.tol = 0.4,
  rt.match.tol = 90,
  polarity = "negative",
  ce = "all",
  column = "hilic",
  total.score.tol = 0,
  candidate.num = 3,
  database = "massbankDatabase0.0.1",
  threads = 10
)


parameter4 <- metIdentifyParam(
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 60,
  ms1.match.ppm = 25,
  ms2.match.tol = 0.4,
  rt.match.tol = 90,
  polarity = "negative",
  ce = "all",
  column = "hilic",
  total.score.tol = 0,
  candidate.num = 3,
  database = "monaDatabase0.0.1",
  threads = 10
)

parameter5 <- metIdentifyParam(
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 60,
  ms1.match.ppm = 25,
  ms2.match.tol = 0.4,
  rt.match.tol = 90,
  polarity = "negative",
  ce = "all",
  column = "hilic",
  total.score.tol = 0,
  candidate.num = 3,
  database = "orbitrapDatabase0.0.1",
  threads = 10
)

parameter6 <- mzIdentifyParam(
  ms1.match.ppm = 25,
  polarity = "negative",
  column = "hilic",
  candidate.num = 3,
  database = "HMDB.metabolite.data",
  threads = 10
)
ms1.data <- grep("\\.csv", dir(), value = TRUE)
ms2.data <- grep("mgf", dir(), value = TRUE)
result.nHILIC <- metIdentify4all(
  ms1.data = ms1.data,
  ms2.data = ms2.data,
  parameter.list = c(
    parameter1,
    parameter2,
    parameter3,
    parameter4,
    parameter5,
    parameter6
  ),
  path = "."
)

save(result.nHILIC, file = "result.nHILIC")
print(1)

####RPLC
##positive
setwd("/home/leeba/metID_Microsampling_247/RPLC/pos")

library(metID)

parameter1 <- metIdentifyParam(
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 60,
  ms1.match.ppm = 25,
  ms2.match.tol = 0.4,
  rt.match.tol = 90,
  polarity = "positive",
  ce = "all",
  column = "rp",
  total.score.tol = 0,
  candidate.num = 3,
  database = "msDatabase_rplc0.0.1",
  threads = 10
)

parameter2 <- metIdentifyParam(
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 60,
  ms1.match.ppm = 25,
  ms2.match.tol = 0.4,
  rt.match.tol = 90,
  polarity = "positive",
  ce = "all",
  column = "rp",
  total.score.tol = 0,
  candidate.num = 3,
  database = "hmdbDatabase0.0.1",
  threads = 10
)

parameter3 <- metIdentifyParam(
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 60,
  ms1.match.ppm = 25,
  ms2.match.tol = 0.4,
  rt.match.tol = 90,
  polarity = "positive",
  ce = "all",
  column = "rp",
  total.score.tol = 0,
  candidate.num = 3,
  database = "massbankDatabase0.0.1",
  threads = 10
)


parameter4 <- metIdentifyParam(
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 60,
  ms1.match.ppm = 25,
  ms2.match.tol = 0.4,
  rt.match.tol = 90,
  polarity = "positive",
  ce = "all",
  column = "rp",
  total.score.tol = 0,
  candidate.num = 3,
  database = "monaDatabase0.0.1",
  threads = 10
)


parameter5 <- metIdentifyParam(
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 60,
  ms1.match.ppm = 25,
  ms2.match.tol = 0.4,
  rt.match.tol = 90,
  polarity = "positive",
  ce = "all",
  column = "rp",
  total.score.tol = 0,
  candidate.num = 3,
  database = "orbitrapDatabase0.0.1",
  threads = 10
)

parameter6 <- mzIdentifyParam(
  ms1.match.ppm = 25,
  polarity = "positive",
  column = "rp",
  candidate.num = 3,
  database = "HMDB.metabolite.data",
  threads = 10
)

ms1.data <- grep("\\.csv", dir(), value = TRUE)
ms2.data <- grep("mgf", dir(), value = TRUE)
result.pRPLC <- metIdentify4all(
  ms1.data = ms1.data,
  ms2.data = ms2.data,
  parameter.list = c(
    parameter1,
    parameter2,
    parameter3,
    parameter4,
    parameter5,
    parameter6
  ),
  path = "."
)

save(result.pRPLC, file = "result.pRPLC")
print(1)

#####RPLC negative
setwd("/home/leeba/metID_Microsampling_247/RPLC/neg")

parameter1 <- metIdentifyParam(
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 60,
  ms1.match.ppm = 25,
  ms2.match.tol = 0.4,
  rt.match.tol = 90,
  polarity = "negative",
  ce = "all",
  column = "rp",
  total.score.tol = 0,
  candidate.num = 3,
  database = "msDatabase_rplc0.0.1",
  threads = 10
)

parameter2 <- metIdentifyParam(
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 60,
  ms1.match.ppm = 25,
  ms2.match.tol = 0.4,
  rt.match.tol = 90,
  polarity = "negative",
  ce = "all",
  column = "rp",
  total.score.tol = 0,
  candidate.num = 3,
  database = "hmdbDatabase0.0.1",
  threads = 10
)

parameter3 <- metIdentifyParam(
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 60,
  ms1.match.ppm = 25,
  ms2.match.tol = 0.4,
  rt.match.tol = 90,
  polarity = "negative",
  ce = "all",
  column = "rp",
  total.score.tol = 0,
  candidate.num = 3,
  database = "massbankDatabase0.0.1",
  threads = 10
)

parameter4 <- metIdentifyParam(
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 60,
  ms1.match.ppm = 25,
  ms2.match.tol = 0.4,
  rt.match.tol = 90,
  polarity = "negative",
  ce = "all",
  column = "rp",
  total.score.tol = 0,
  candidate.num = 3,
  database = "monaDatabase0.0.1",
  threads = 10
)

parameter5 <- metIdentifyParam(
  ms1.ms2.match.mz.tol = 25,
  ms1.ms2.match.rt.tol = 60,
  ms1.match.ppm = 25,
  ms2.match.tol = 0.4,
  rt.match.tol = 90,
  polarity = "negative",
  ce = "all",
  column = "rp",
  total.score.tol = 0,
  candidate.num = 3,
  database = "orbitrapDatabase0.0.1",
  threads = 10
)

parameter6 <- mzIdentifyParam(
  ms1.match.ppm = 25,
  polarity = "negative",
  column = "rp",
  candidate.num = 3,
  database = "HMDB.metabolite.data",
  threads = 10
)

ms1.data <- grep("\\.csv", dir(), value = TRUE)
ms2.data <- grep("mgf", dir(), value = TRUE)
result.nRPLC <- metIdentify4all(
  ms1.data = ms1.data,
  ms2.data = ms2.data,
  parameter.list = c(
    parameter1,
    parameter2,
    parameter3,
    parameter4,
    parameter5,
    parameter6
  ),
  path = "."
)

save(result.nRPLC, file = "result.nRPLC")
print(1)



