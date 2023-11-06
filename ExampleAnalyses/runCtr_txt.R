#' Title
#'
#' @param HCx 
#' @param EP 
#' @param HT 
#' @param PCommunity 
#' @param modelpath 
#' @param sim_dir 
#' @param sim_subpath 
#' @param MCruns 
#' @param SaveEnvironment 
#' @param genHerbFile 
#' @param proj_dir 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
runctr <- function(HCx=0.05,EP="Establishment",HT="R",PCommunity,
                   modelpath="/home/sagemaker-user/Projects/Research/IBCgrassGUI/Model-files/",
                   sim_dir="/home/sagemaker-user/Projects/Research/IBCgrassGUI/ExampleAnalyses/Calthion/",
                   sim_subpath="Test1",MCruns=2,SaveEnvironment,genHerbFile=FALSE,
                   proj_dir="/home/sagemaker-user/Projects/Research/IBCgrassGUI/",
                   useHerbFile="Input-files/HerbFact.txt",...){
  #####
  # read PFT community file and sensitivity file
  #####
  print(sim_subpath)
  ModelVersion <- 3
  # NEW
  ## PFTfileName <- get("IBCcommunity", envir=SaveEnvironment) # txt-file with trait parameters of all species ==> need to change here. 
  PFTfileName <- as.character(PCommunity)
  PFTHerbEffectFile <- paste0(proj_dir,"Input-files/HerbFact.txt") # you need this file if you run IBC with manual effects!! 
  ## you also need to change this file if you want manual effects! You also need 
  
  AppRateFile <- paste0(proj_dir,"Input-files/AppRate.txt") # you need this file if you run IBC with dose-response data
  #MCruns <- get("IBCrepetition", envir=SaveEnvironment) # nb of repetitions 
  
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
  
  
  
  # generate PFTfile with herbicide sensitivities 
  
  PFTfile <- get("IBCcommunityFile", envir=SaveEnvironment)
  PFTsensitivity <- get("PFTSensitivityFile", envir=SaveEnvironment)
  if(genHerbFile) genHerbFactFile(EP=EP,file=paste0(modelpath,"HerbFact.txt")) else{
    ## use the copied Herbfile already copied before the MCruns.
    file.copy(paste0(proj_dir,useHerbFile),  modelpath)
  }
  
  
  copy <- file.copy(paste0(proj_dir,"Input-files/AppRate.txt"),  modelpath) ### There is no need for this correct!! this line could be commented out. 
  ## file.copy(paste0(proj_dir,"Input-files/HerbFact.txt"),  modelpath)
  
  dir.create(paste(sim_dir,sim_subpath,"/HerbicideSettings", sep=""), recursive=TRUE)
  # # cluster for parallel processing ==> go outside of the function
  ## no_cores <- max(detectCores()-2,1) # you might want to adapt this and give a specific number of cores (e.g. when using a HPC)
  no_cores <- 30
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  # change directory
  setwd(modelpath)
  scenario <- 0 
  foreach(MC = 1:MCruns)  %dopar%
    system(paste('wine ./IBCgrassGUI.exe', ModelVersion, GridSize, Tmax, InitDuration, PFTfileName, SeedInput, belowres, abres, abampl, tramp, graz, cut,
                 week_start, HerbDuration, 0, EffectModel, scenario, MC, sep=" "), intern=T)
  stopCluster(cl)
  
 
  #####
 
  #####
  ## dir.create("currentSimulation/0", recursive=TRUE)
  
  dir.create(paste(sim_dir,sim_subpath,"/0", sep=""))
  file_list <- list.files(path = "Model-files/", pattern="Pt__*")
  for (file in file_list){
    path <- paste(sim_dir,sim_subpath,"/0", sep="")
    copy <- file.copy(paste("Model-files/" ,file , sep="") ,  path)
    #  todo make sure that all files were copied!
    if (copy==T) file.remove(paste("Model-files/" ,file , sep="") )              
  } 
  # GRD files
  file_list <- list.files(path = "Model-files/", pattern="Grd__*")
  for (file in file_list){
    path <- paste(sim_dir,sim_subpath,"/0", sep="")
    copy <- file.copy(paste("Model-files/" ,file , sep="") ,  path)
    #  todo make sure that all files were copied!
    if (copy==T) file.remove(paste("Model-files/" ,file , sep="") )              
  }
  
  file.remove(paste0(modelpath,"HerbFact.txt"))
  remove <- file.remove(paste0(modelpath,"AppRate.txt"))
  # copy <- file.copy(paste0(modelpath,"PFTsensitivity.txt") ,  paste(sim_dir,sim_subpath,"/HerbicideSettings", sep=""))
  # if (all(copy==T)) file.remove("PFTsensitivity.txt")
}