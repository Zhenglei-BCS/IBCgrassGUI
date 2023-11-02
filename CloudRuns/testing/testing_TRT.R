#####
# libraries
#####
library(foreach)
library(doParallel)
library(labeling)
library(dplyr)
library(readr)



source("~/Projects/Research/IBCgrassGUI/ExampleAnalyses/runTrt_txt.R")

load("~/Projects/Research/IBCgrassGUI/ExampleAnalyses/TXTfile/HerbicideSettings/SimulationSettings.Rdata") # j

proj_dir <- "~/Projects/Research/IBCgrassGUI/"
# set working directory  ----
setwd(proj_dir)
# Run Original Sensitivity Simulation ----
print("Run Original Simulation")
t0 <- Sys.time()
set.seed(100)
runtrt(HCx=0,EP="Biomass",HT="R",PCommunity="Fieldedge.txt",modelpath="/home/sagemaker-user/Projects/Research/IBCgrassGUI/Model-files/",
       sim_dir="/home/sagemaker-user/Projects/Research/IBCgrassGUI/currentSimulation/",
       sim_subpath="orig",MCruns=1,SaveEnvironment=SaveEnvironment,genHerbFile = FALSE,
       userHerbFile="ExampleAnalyses/Calthion/SimEnvironments/C/HerbicideSettings/HerbFact.txt",
       proj_dir = "~/Projects/Research/IBCgrassGUI/")
t1 <- Sys.time()
t1-t0

print(t1-t0)

# Run Adapated Sensitivity Simulation ----
print("Run Adapated Sensitivity Simulation")
t0 <- Sys.time()
set.seed(100)
runtrt(HCx=0.05,EP="Biomass",HT="R",PCommunity="Fieldedge.txt",
       modelpath="/home/sagemaker-user/Projects/Research/IBCgrassGUI/Model-files/",
       sim_dir="/home/sagemaker-user/Projects/Research/IBCgrassGUI/currentSimulation/",
       sim_subpath="gen",MCruns=1,SaveEnvironment,genHerbFile = TRUE,
       proj_dir = "~/Projects/Research/IBCgrassGUI/",
       useHerbFile = "ExampleAnalyses/Calthion/SimEnvironments/C/HerbicideSettings/HerbFact.txt")
t1 <- Sys.time()
t1-t0

print(t1-t0)

### output from this ##

# > data.table::fread("currentSimulation/gen/HerbicideSettings/Fieldedge1.txt") -> tmp2
# > tmp2$sens
# [1] 0.216919344 0.219075822 0.227498355 0.665467527 0.121318733 0.054908828 0.085731601 0.045253329
# [9] 0.068408890 0.192601006 0.214914352 0.046065066 0.084080641 0.222508023 0.111279679 0.310427702
# [17] 0.184106760 0.209508249 0.003949584 0.101568340 0.166562684 0.237861260 0.220573405 0.018331022
# [25] 0.061024183 0.172091878 0.045354041 0.208009436 0.127217866 0.143545019 0.015563963 0.177642588
# [33] 0.144173044 0.144203149 0.163221211 0.213723735 0.242231218 0.112816293 0.225374764 0.109548844
# [41] 0.050423939 0.107195776 0.036482604 0.199854404 0.631442592 0.100905243 0.241226872 0.029865624
# [49] 0.112078663 0.225921063 0.107047191 0.082238197
# > sum(tmp2$sens>0.25)/length(tmp2$sens)
# [1] 0.05769231
# > data.table::fread("currentSimulation/orig/HerbicideSettings/Fieldedge1.txt") -> tmp2
# > sum(tmp2$sens>0.25)/length(tmp2$sens)
# [1] 0.75

data.table::fread("currentSimulation/gen/1/Grd__type_treatment_Scenario_1_MCrun_1.txt") -> tmp2
head(tmp2)
tail(tmp2)
data.table::fread("currentSimulation/orig/1/Grd__type_treatment_Scenario_1_MCrun_1.txt") -> tmp1
head(tmp1)
tail(tmp1)

data.table::fread("currentSimulation/gen/1/Pt__type_treatment_Scenario_1_MCrun_1.txt") -> tmp2
head(tmp2)
tail(tmp2)
data.table::fread("currentSimulation/orig/1/Pt__type_treatment_Scenario_1_MCrun_1.txt") -> tmp1
head(tmp1)
tail(tmp1)
dim(tmp1)
length(unique(tmp1$Time))
nrow(tmp1)/length(unique(tmp1$Time))
table(tmp1$PFT)
table(tmp2$PFT)
##########################################