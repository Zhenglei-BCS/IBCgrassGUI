# @Ismael
# CHANGES TO THE CODE
# 1. There are two new files loaded:
#     1.1. "Exp.Matrix1.csv" that provides several state variables that define each run: HCx, EP, HT, PCommunity.
#     "Exp.Matrix1" is just one line of "Exp.Matrix1.csv", that correspond to 1 experiment.
#     @Zhenglei: We need a way to input the experimental number of each run in the clowd system
#     1.2. "Cotyledon.csv" provides cotyledon information by sp # Not yet implemented

# 3. New effect classes are created and (>25%, >25%) to assign effect levels to sp below the HCx threshold of interest
# 4. New HT class handeling needed (HT: herbicide tolerance: random("R"), monocot("TMono"), dicot("TDicot"))
#    The idea is that when sensitivity is assign at random, the effect on tolerant is always fix at 0 
#    and count for the total n of species in class "<25%"


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
dir.create(paste0(sim_dir,sim_subpath,"/0"), recursive=TRUE) # adapt path if needed

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

#####
# running treatment based on txt file (not in this example)
#####

if(HerbEff=="txt-file"){
  # cluster for parallel processing
  no_cores <- max(detectCores()-2,1) # you might want to adapt this and give a specific number of cores (e.g. when using a HPC)
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  # generate PFTfile with herbicide sensitivities 
  
  PFTfile <- get("IBCcommunityFile", envir=SaveEnvironment)
  
  # run treatment simulations
  foreach(MC = 1:MCruns, .export=c("PFTfile", "PFTsensitivity", "PFTfileName", "EffectModel",
                                   "ModelVersion", "belowres", "abres", "abampl", "Tmax", "InitDuration", "GridSize", "SeedInput",
                                   "HerbDuration", "tramp", "graz", "cut", "week_start"))  %dopar%
    {
      scenario <- 1 # for treatment simulations
      
      # Generate PFTfile for each MC run
      # NEW - random assignment per MC run of PFT sensitivity ####
      # @Zhenglei should we set.seed???
      # consider HT
      if (HT=="R"){
        # NEW, PFT sensitivity info ####
        # PFTsensitivity <- get("PFTSensitivityFile", envir=SaveEnvironment) # here the sensitivity of all PFTs are stored
        PFTsensitivity <- PFTfile[,c(1:2)] # NEW This is nx2 df species(PFTs), Sensitivity(factor)
        
        # Initial assignment of sensitivite based on HCx (e.g., HCx=005 =>5% =>"High", 95% =>"Low")
        PFTsensitivity$Sensitivity = NA
        PFTsensitivity$Sensitivity[1:round(HCx*nrow(PFTsensitivity),0)] = ">25%"
        PFTsensitivity$Sensitivity[(round(HCx*nrow(PFTsensitivity),0)+1):nrow(PFTsensitivity)] = "<25%"
        PFTsensitivity = PFTsensitivity[2:3] #remove ID column to get a total of two columns for later merge correctly with PFTfile
        # Combined with HerbFact.txt generates all EP-HCx combinations
        # random assignment of sp and sensitivities per MC run
        # <25% = random(0-25%)) (Or sampling an SSD between 0-25%)
        # >25% = random (25% - 100%) (Or sampling an SSD between 25%-100%)
        
        # Assign seisitivity at random
        rows <- sample(nrow(PFTsensitivity))
        PFTsensitivity = PFTsensitivity[rows, ]
        PFTfile<-merge(PFTfile, PFTsensitivity, by="Species")
      }
      if (HT=="TMono"){} # need to split by mono/dicot, Mono get PFTfile$Sensitivity== 0, dicot get the rest (0 count to sum($Sensitivity=="<25%")) 
      if (HT=="TDicot"){} # need to split by mono/dicot, Dicot get PFTfile$Sensitivity== 0, monocot get the rest (0 count to sum($Sensitivity=="<25%"))
      
      
      #make sure, all values are set to 0 (no affect)
      PFTfile[,25] <- 0
      
      # NEW - effect categories #### @Zhenglei this could be replace by sampling from the SSD probability distribution?
      # set <25% values
      PFTfile[PFTfile$Sensitivity=="<25%",25] <- c(runif(nrow(PFTfile[PFTfile$Sensitivity=="<25%",]),  min = 0, max = 0.25))
      # set >25% values
      PFTfile[PFTfile$Sensitivity==">25%",25] <- c(runif(nrow(PFTfile[PFTfile$Sensitivity==">25%",]),  min = 0.25, max = 1))
      
      
      # set random values
      PFTfile[PFTfile$Sensitivity=="random",25] <- c(runif(nrow(PFTfile[PFTfile$Sensitivity=="random",]),  min = 0, max = 1))
      
      # set full values
      PFTfile[PFTfile$Sensitivity=="full",25] <- c(rep(1, nrow(PFTfile[PFTfile$Sensitivity=="full",])))
      
      # set high values
      PFTfile[PFTfile$Sensitivity=="high",25] <- c(runif(nrow(PFTfile[PFTfile$Sensitivity=="high",]),  min = 0.66, max = 1))
      
      # set medium values
      PFTfile[PFTfile$Sensitivity=="medium",25] <- c(runif(nrow(PFTfile[PFTfile$Sensitivity=="medium",]),  min = 0.35, max = 0.65))
      
      # set low values
      PFTfile[PFTfile$Sensitivity=="low",25] <- c(runif(nrow(PFTfile[PFTfile$Sensitivity=="low",]),  min = 0.1, max = 0.35))
      
      # copy others (not affected PFTs get 0)
      PFTfile[PFTfile$Sensitivity=="not affected",25] <- c(rep(0, nrow(PFTfile[PFTfile$Sensitivity=="not affected",])))
      
      # remove temp. column and prepare final PFT file for this repetition
      PFTfile<-PFTfile[,-ncol(PFTfile)]
      PFTfile <- cbind(PFTfile[,c(2,1)],PFTfile[,-c(1:2)])
      
      #save PFT file
      write.table(PFTfile[,-ncol(PFTfile)], paste(unlist(strsplit(PFTfileName,".txt")), MC, ".txt", sep=""), row.names=F, quote=F, sep="\t")
      
      # copy necessary files to Model-files folder
      path <- "Model-files/"
      copy <- file.copy(paste(unlist(strsplit(PFTfileName,".txt")), MC, ".txt", sep=""),  path)
      copy <- file.copy("HerbFact.txt",  path)
      copy <- file.copy("Input-files/AppRate.txt",  path)
      
      # change directory
      setwd('Model-files')
      
      # cmd call for running IBC
      mycall<-paste('./IBCgrassGUI', ModelVersion, GridSize, Tmax, InitDuration, paste("./",unlist(strsplit(PFTfileName,".txt")), MC, ".txt", sep=""), 
                    SeedInput, belowres, abres, abampl, tramp, graz, cut, 
                    week_start, HerbDuration, 1, EffectModel, scenario, MC, sep=" ")
      
      # start treatment run
      system(mycall, intern=TRUE)
      
      # change working directory
      setwd('..')
      
      # remove old files
      remove <- file.remove(paste("Model-files/", unlist(strsplit(PFTfileName,".txt")), MC, ".txt", sep=""))
      remove <- file.remove("Model-files/HerbFact.txt")
      remove <- file.remove("Model-files/AppRate.txt")
    }
  
  # stop cluster for parallel processing
  stopCluster(cl)
  
  #####
  # copy treatment
  #####
  # copy files to directory
  ## dir.create("currentSimulation/1", recursive=TRUE)
  dir.create(paste(sim_dir,sim_subpath,"/1", sep=""))
  file_list <- list.files(path = "Model-files/", pattern="Pt__*")
  for (file in file_list){
    path <- paste(sim_dir,sim_subpath,"/1", sep="")
    copy <- file.copy(paste("Model-files/" ,file , sep="") ,  path)
    #  todo make sure that all files were copied!
    if (copy==T) file.remove(paste("Model-files/" ,file , sep="") )              
  } 
  # GRD files
  file_list <- list.files(path = "Model-files/", pattern="Grd__*")
  for (file in file_list){
    path <- paste(sim_dir,sim_subpath,"/1", sep="")
    copy <- file.copy(paste("Model-files/" ,file , sep="") ,  path)
    #  todo make sure that all files were copied!
    if (copy==T) file.remove(paste("Model-files/" ,file , sep="") )              
  }
  
  # Copy treamtent settings
  dir.create(paste(sim_dir,sim_subpath,"/HerbicideSettings", sep=""), recursive=TRUE)
  file_list <- list.files(pattern=paste(unlist(strsplit(PFTfileName,".txt")),"*",sep=""))
  file_list <- file_list[file_list!=PFTfileName] #TODO: stimmt das mit der Datei???
  file_list <- c(file_list, "SimulationSettings.Rdata")
  copy <- file.copy(file_list ,  paste(sim_dir,sim_subpath,"/HerbicideSettings", sep=""))
  #  todo make sure that all files were copied!
  if (all(copy==T)) file.remove(file_list)
  copy <- file.copy("HerbFact.txt" ,  paste(sim_dir,sim_subpath,"/HerbicideSettings", sep=""))
  if (all(copy==T)) file.remove("HerbFact.txt")
  copy <- file.copy("PFTsensitivity.txt" ,  paste(sim_dir,sim_subpath,"/HerbicideSettings", sep=""))
  if (all(copy==T)) file.remove("PFTsensitivity.txt")
}
