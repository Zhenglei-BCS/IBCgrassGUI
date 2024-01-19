source("~/Projects/Research/IBCgrassGUI/ExampleAnalyses/runTrt_txt.R")


#####
# libraries
#####
library(foreach)
library(doParallel)
library(labeling)
library(dplyr)
#####
# Notes
#####
# -------------- Make sure to compile IBC before running this code ----------- #
# ----------------- run Model-files/CompileIBC.bat --------------------------- #
# code is copied/adapted from R-files/SimulationSpecificSettings.R
#####
# Load previously saved simulation settings 
#####
# You can use the GUI to generate the SimulationSettings.Rdata
# simply run a very small test script (min. repetitions, min. years, min. grid size) and save the project
# You then might want to adapt some of the settings manually (e.g. increase repetitions, number of simulated years etc.)

#NEW ####
library(readr)
proj_dir <- "/home/sagemaker-user/Projects/Research/IBCgrassGUI"
sim_dir <- paste0(proj_dir,"/ExampleAnalyses/verification/")
path <- "Model-files/"
modelpath <- paste0(proj_dir,"/Model-files/")

setwd(proj_dir)

Exp.Matrix <- read_csv(paste0(proj_dir,"/ExampleAnalyses/Exp.Matrix.csv")) # %>% mutate(RUN=1:n())

ind <- which(Exp.Matrix$PCommunity=="Arrhenatheretalia.txt" & Exp.Matrix$EP%in% c("Biomass"))

Exp.Matrix <- Exp.Matrix[ind,]
Exp.Matrix <- Exp.Matrix%>%filter(HCx==0.05)
Exp.Matrix <- rbind(Exp.Matrix,
                    structure(list(HCx = 0.95, EP = "Biomass", HT = "R", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame")),
                    structure(list(HCx = 0.95, EP = "Biomass", HT = "Fixed", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame")),
                    structure(list(HCx = 0.05, EP = "Biomass", HT = "Fixed", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame")),
                    structure(list(HCx = 1, EP = "Biomass", HT = "R", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame")),
                    structure(list(HCx = 1, EP = "Biomass", HT = "CompletelyR", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame")),
                    structure(list(HCx = 0.95, EP = "Mortality", HT = "R", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame")),
                    structure(list(HCx = 0.95, EP = "Mortality", HT = "Fixed", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame")),
                    structure(list(HCx = 0.05, EP = "Mortality", HT = "Fixed", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame")),
                    structure(list(HCx = 1, EP = "Mortality", HT = "R", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame")),
                    structure(list(HCx = 1, EP = "Mortality", HT = "CompletelyR", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame")),
                    structure(list(HCx = 0.95, EP = "Establishment", HT = "R", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame")),
                    structure(list(HCx = 0.95, EP = "Establishment", HT = "Fixed", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame")),
                    structure(list(HCx = 0.05, EP = "Establishment", HT = "Fixed", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame")),
                    structure(list(HCx = 1, EP = "Establishment", HT = "R", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame")),
                    structure(list(HCx = 1, EP = "Establishment", HT = "CompletelyR", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame")),
                    ### SeedlingBiomass	Establishment	SeedSterility	SeedNumber
                    structure(list(HCx = 0.95, EP = "SeedlingBiomass", HT = "Fixed", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame")),
                    structure(list(HCx = 0.05, EP = "SeedlingBiomass", HT = "Fixed", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame")),
                    structure(list(HCx = 0.95, EP = "SeedSterility", HT = "Fixed", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame")),
                    structure(list(HCx = 0.05, EP = "SeedSterility", HT = "Fixed", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame")),
                    structure(list(HCx = 0.95, EP = "SeedNumber", HT = "Fixed", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame")),
                    structure(list(HCx = 0.05, EP = "SeedNumber", HT = "Fixed", PCommunity = "Arrhenatheretalia.txt"), row.names = c(NA,-1L), 
                              class = c("tbl_df", "tbl", "data.frame"))
                    )
Exp.Matrix1 <- Exp.Matrix

Exp.Matrix <- Exp.Matrix %>% mutate(RUN=1:n())
#write_csv(Exp.Matrix,"ExampleAnalyses/verification/Exp.Matrix.csv")
ind <- 17:22
# @Zhenglei we will need to load environment based on Exp.Matrix1$PCommunity
load(paste0(sim_dir,"/SimEnvironments/A/HerbicideSettings/SimulationSettings.Rdata"))# just as an example for this script

###################################
## Change the simulation time
#week_start <- get("IBCweekstart", envir=SaveEnvironment)-10 # start of the herbicide application (calendar week); note that IBC only simulate growing period
SaveEnvironment$IBCweekstart <- 11
SaveEnvironment$IBCRecovery <- 2 # recovery duration [years]
InitDuration <- get("IBCInit", envir=SaveEnvironment) 
SaveEnvironment$IBCInit <- 5 # intial years 
SaveEnvironment$BiomassEff
SaveEnvironment$BiomassEffFile
SaveEnvironment$SurvivalEff
SaveEnvironment$EstablishmentEff
SaveEnvironment$EstablishmentEffFile
SaveEnvironment$SeedNumberEff
SaveEnvironment$EffectData
SaveEnvironment$IBCrepetition <- 20
#Tmax <- InitDuration + HerbDuration + RecovDuration # years to simulate

start_time <- Sys.time()

for(i in ind){
  
  ######################################################################################
  HCx <- Exp.Matrix$HCx[i]
  EP <- Exp.Matrix$EP[i]
  HT <- Exp.Matrix$HT[i]
  PCommunity <- Exp.Matrix$PCommunity[i]
  setwd(proj_dir)
  runtrt(HCx=HCx,EP=EP,HT=HT,PCommunity=PCommunity,
         modelpath="/home/sagemaker-user/Projects/Research/IBCgrassGUI/Model-files/",
         sim_dir="/home/sagemaker-user/Projects/Research/IBCgrassGUI/ExampleAnalyses/verification/",
         sim_subpath=paste0("Run",i),MCruns=20,SaveEnvironment=SaveEnvironment,
         genHerbFile = TRUE,
         proj_dir = "~/Projects/Research/IBCgrassGUI/",
         useHerbFile = "ExampleAnalyses/Arrhenatheretalia/SimEnvironments/A/HerbicideSettings/HerbFact.txt")
  file.remove(paste0(proj_dir,"/",grep("Arrhenatheretalia*",list.files(proj_dir),value=T)))
}

end_time <- Sys.time()
end_time - start_time
timeinfo <- end_time-start_time
save(timeinfo,file=paste0(sim_dir,"/timeinfo.rda"))

