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
sim_dir <- paste0(proj_dir,"/ExampleAnalyses/Arrhenatheretalia/")
path <- "Model-files/"
modelpath <- paste0(proj_dir,"/Model-files/")



Exp.Matrix <- read_csv(paste0(proj_dir,"/ExampleAnalyses/Calthion/Exp.Matrix.csv")) %>% mutate(RUN=1:n())
which(Exp.Matrix$PCommunity=="Arrhenatheretalia.txt" & Exp.Matrix$EP=="Biomass")

i <- 1 # @zhenglei we need a way to subset according to the run number
sim_subpath <- paste0("Run",i)
# @Zhenglei we will need to load environment based on Exp.Matrix1$PCommunity
load(paste0(sim_dir,"/SimEnvironments/A/HerbicideSettings/SimulationSettings.Rdata"))# just as an example for this script

HCx = Exp.Matrix$HCx[i] # Provides HCx (new parameter)
EP = Exp.Matrix$EP[i] # can replace PFTHerbEffectFile
HT = Exp.Matrix$HT[i] # Provides HT (new parameter)
PCommunity = Exp.Matrix$PCommunity[i]

#####
# read PFT community file and sensitivity file
#####
ModelVersion <- 3
# NEW
## PFTfileName <- get("IBCcommunity", envir=SaveEnvironment) # txt-file with trait parameters of all species ==> need to change here. 
PFTfileName <- as.character(Exp.Matrix[i,"PCommunity"])
PFTHerbEffectFile <- "Input-files/HerbFact.txt" # you need this file if you run IBC with manual effects!! 
## you also need to change this file if you want manual effects! You also need 

AppRateFile <- "Input-files/AppRate.txt" # you need this file if you run IBC with dose-response data
MCruns <- get("IBCrepetition", envir=SaveEnvironment) # nb of repetitions 

GridSize <- get("IBCgridsize", envir=SaveEnvironment) # area of the grid
SeedInput <- get("IBCSeedInput", envir=SaveEnvironment) # external seed input
belowres <- get("IBCbelres", envir=SaveEnvironment) # belowground resources
abres <- get("IBCabres", envir=SaveEnvironment) # aboveground resources
abampl <- get("IBCabampl", envir=SaveEnvironment) # seasonal amplitude of aboveground resources

graz <- get("IBCgraz", envir=SaveEnvironment) # grazing intensity
tramp <- get("IBCtramp", envir=SaveEnvironment) # trampling intensity
cut <- get("IBCcut", envir=SaveEnvironment) # cutting/mowing events

week_start <- get("IBCweekstart", envir=SaveEnvironment)-10 # start of the herbicide application (calendar week); note that IBC only simulate growing period
HerbDuration <- get("IBCDuration", envir=SaveEnvironment) # herbicide duration [years]
RecovDuration <- get("IBCRecovery", envir=SaveEnvironment) # recovery duration [years]
InitDuration <- get("IBCInit", envir=SaveEnvironment) # intial years 
Tmax <- InitDuration + HerbDuration + RecovDuration # years to simulate
HerbEff <- get("IBCherbeffect", envir=SaveEnvironment) # txt or dose response
if(HerbEff=="txt-file") EffectModel <- 0
if(HerbEff=="dose-response") EffectModel <- 2 ## What does this mean??
Scenarios <- as.numeric(get("IBCScenarios", envir=SaveEnvironment)) # number of rates tested
nb_data <- as.numeric(get("nb_data", envir=SaveEnvironment)) # number of test species in the C case 6

#####
# running control
#####
scenario <- 0 # for control runs
# copy community file into Model-folder

write.table(get("IBCcommunityFile", envir=SaveEnvironment), file=paste(modelpath,PFTfileName, sep=""), sep="\t", quote=F,row.names=F) # make sure the Model-files folder includes the file of your PFT community 
setwd(proj_dir)
copy <- file.copy("Input-files/HerbFact.txt",  modelpath) # make sure the Model-files folder includes the HerbFact.txt file 
copy <- file.copy("Input-files/AppRate.txt",  modelpath) # make sure the Model-files folder includes the AppRate.txt file

# Change directory
setwd('Model-files')
# determine the number of cores that can be used (for parallel running)
no_cores <- max(detectCores()-2,1) # you might want to adapt this and give a specific number of cores (e.g. when using a HPC)
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# Start control simulations in parallel with the given settings
MCruns <- 30
foreach(MC = 1:MCruns)  %dopar%
  system(paste('wine ./IBCgrassGUI.exe', ModelVersion, GridSize, Tmax, InitDuration, PFTfileName, SeedInput, belowres, abres, abampl, tramp, graz, cut,
               week_start, HerbDuration, 0, EffectModel, scenario, MC, sep=" "), intern=T)
stopCluster(cl)

# change directory
setwd('..')

# remove files before running treatment
remove <- file.remove(paste(modelpath, PFTfileName, sep=""))
remove <- file.remove("Model-files/HerbFact.txt")
remove <- file.remove("Model-files/AppRate.txt")

#####
# copy control simulations
#####
# create folder
if(!dir.exists(paste0(sim_dir,sim_subpath,"/0"))) dir.create(paste0(sim_dir,sim_subpath,"/0"), recursive=TRUE) # adapt path if needed

# PFT files
file_list <- list.files(path = "Model-files/", pattern="Pt__*")
for (file in file_list){
  path <- paste(sim_dir,sim_subpath,"/", unlist(strsplit(file,"_"))[6], sep="")
  copy <- file.copy(paste("Model-files/" ,file , sep="") ,  path)
  if (copy==T) file.remove(paste("Model-files/" ,file , sep="") )
}

# GRD files
file_list <- list.files(path = "Model-files/", pattern="Grd__*", recursive=TRUE)
for (file in file_list){
  path <- paste(sim_dir,sim_subpath,"/", unlist(strsplit(file,"_"))[6], sep="")
  copy <- file.copy(paste("Model-files/" ,file , sep="") ,  path)
  if (copy==T) file.remove(paste("Model-files/" ,file , sep="") )
}
