library(foreach)
library(doParallel)
library(labeling)
library(dplyr)
source("~/Projects/Research/IBCgrassGUI/ExampleAnalyses/runTrt_txt.R")
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
sim_dir <- paste0(proj_dir,"/ExampleAnalyses/Calthion/")
path <- "Model-files/"
modelpath <- paste0(proj_dir,"/Model-files/")



Exp.Matrix <- read_csv(paste0(proj_dir,"/ExampleAnalyses/Calthion/Exp.Matrix.csv")) %>% mutate(RUN=1:n())
which(Exp.Matrix$PCommunity=="Calthion.txt" & Exp.Matrix$EP=="Establishment")
i <- 61 # @zhenglei we need a way to subset according to the run number
sim_subpath <- paste0("Run",i)
# @Zhenglei we will need to load environment based on Exp.Matrix1$PCommunity
load(paste0(sim_dir,"/SimEnvironments/C/HerbicideSettings/SimulationSettings.Rdata"))# just as an example for this script

HCx = Exp.Matrix$HCx[i] # Provides HCx (new parameter)
EP = Exp.Matrix$EP[i] # can replace PFTHerbEffectFile
HT = Exp.Matrix$HT[i] # Provides HT (new parameter)
PCommunity = Exp.Matrix$PCommunity[i]
start_time <- Sys.time()
setwd(proj_dir)
runtrt(HCx=HCx,EP=EP,HT=HT,PCommunity=PCommunity,modelpath="/home/sagemaker-user/Projects/Research/IBCgrassGUI/Model-files/",
       sim_dir="/home/sagemaker-user/Projects/Research/IBCgrassGUI/ExampleAnalyses/Calthion/",
       sim_subpath="Run61",MCruns=30,SaveEnvironment=SaveEnvironment)
end_time <- Sys.time()
end_time - start_time
timeinfo <- end_time-start_time
save(timeinfo,file=paste0(sim_dir,sim_subpath,"/timeinfo.rda"))