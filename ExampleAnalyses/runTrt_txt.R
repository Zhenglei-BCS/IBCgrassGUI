genHerbFactFile <- function(EP="Establishment",file=NULL){
  allEP <- structure(list(Biomass = 0L, Mortality = 0L, SeedlingBiomass = 0L, 
                          Establishment = 0L, SeedSterility = 0L, SeedNumber = 0L), class = "data.frame", row.names = c(NA, 
                                                                                                                        -1L))
  allEP[EP] <- 1
  if(!is.null(file)) write.table(allEP,file,row.names = F,quote = F)
  return(EP)
}

genSensitivity <- function(n=100,HCx=0.05,HT="R"){
  if(HT=="R"){
    id <- 1:n
    id1 <- 1:round(HCx*n,0)
    id2 <- (round(HCx*n,0)+1):n
    sensi <- rep(NA,n)
    sensi[id1] <- ">25%"
    sensi[id2] <- "<25%"
    sensi1 <- sensi[sample(n)]
  }
  
  return(sensi1)
}

runtrt <- function(HCx=0.05,EP="Establishment",HT="R",PCommunity,
                   modelpath="/home/sagemaker-user/Projects/Research/IBCgrassGUI/Model-files/",
                   sim_dir="/home/sagemaker-user/Projects/Research/IBCgrassGUI/ExampleAnalyses/Calthion/",
                   sim_subpath="Run61",MCruns=2,SaveEnvironment,...){
  #####
  # read PFT community file and sensitivity file
  #####
  ModelVersion <- 3
  # NEW
  ## PFTfileName <- get("IBCcommunity", envir=SaveEnvironment) # txt-file with trait parameters of all species ==> need to change here. 
  PFTfileName <- as.character(PCommunity)
  PFTHerbEffectFile <- "Input-files/HerbFact.txt" # you need this file if you run IBC with manual effects!! 
  ## you also need to change this file if you want manual effects! You also need 
  
  AppRateFile <- "Input-files/AppRate.txt" # you need this file if you run IBC with dose-response data
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
  
  
  # # cluster for parallel processing ==> go outside of the function
  ## no_cores <- max(detectCores()-2,1) # you might want to adapt this and give a specific number of cores (e.g. when using a HPC)
  no_cores <- 30
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  # generate PFTfile with herbicide sensitivities 
  
  PFTfile <- get("IBCcommunityFile", envir=SaveEnvironment)
  PFTsensitivity <- get("PFTSensitivityFile", envir=SaveEnvironment)
  genHerbFactFile(EP=EP,file=paste0(modelpath,"HerbFact.txt"))
  dir.create(paste(sim_dir,sim_subpath,"/HerbicideSettings", sep=""), recursive=TRUE)
  # run treatment simulations
  foreach(MC = 1:MCruns, .export=c("PFTfile", "PFTsensitivity", "PFTfileName", "EffectModel",
                                   "ModelVersion", "belowres", "abres", "abampl", "Tmax", "InitDuration", "GridSize", "SeedInput",
                                   "HerbDuration", "tramp", "graz", "cut", "week_start",
                                   "genSensitivity","genHerbFactFile"))  %dopar%
    {
      scenario <- 1 # for treatment simulations
      
      # Generate PFTfile for each MC run
      # NEW - random assignment per MC run of PFT sensitivity ####
      # @Zhenglei should we set.seed???
      # consider HT
      if (HT=="R"){
        # NEW, PFT sensitivity info ####
        
        # Initial assignment of sensitivite based on HCx (e.g., HCx=005 =>5% =>"High", 95% =>"Low")
        #PFTsensitivity$Sensitivity = genSensitivity(n=nrow(PFTsensitivity),HCx=HCx,HT="R")
        PFTfile$Sensitivity = genSensitivity(n=nrow(PFTfile),HCx=HCx,HT="R")
        
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
      
      
      # # set random values
      # PFTfile[PFTfile$Sensitivity=="random",25] <- c(runif(nrow(PFTfile[PFTfile$Sensitivity=="random",]),  min = 0, max = 1))
      # 
      # # set full values
      # PFTfile[PFTfile$Sensitivity=="full",25] <- c(rep(1, nrow(PFTfile[PFTfile$Sensitivity=="full",])))
      # 
      # # set high values
      # PFTfile[PFTfile$Sensitivity=="high",25] <- c(runif(nrow(PFTfile[PFTfile$Sensitivity=="high",]),  min = 0.66, max = 1))
      # 
      # # set medium values
      # PFTfile[PFTfile$Sensitivity=="medium",25] <- c(runif(nrow(PFTfile[PFTfile$Sensitivity=="medium",]),  min = 0.35, max = 0.65))
      # 
      # # set low values
      # PFTfile[PFTfile$Sensitivity=="low",25] <- c(runif(nrow(PFTfile[PFTfile$Sensitivity=="low",]),  min = 0.1, max = 0.35))
      # 
      # # copy others (not affected PFTs get 0)
      # PFTfile[PFTfile$Sensitivity=="not affected",25] <- c(rep(0, nrow(PFTfile[PFTfile$Sensitivity=="not affected",])))
      # 
      # remove temp. column and prepare final PFT file for this repetition
      PFTfile<-PFTfile[,-ncol(PFTfile)]
      PFTfile <- cbind(PFTfile[,c(2,1)],PFTfile[,-c(1:2)]) ## ??? I don't think it is doing anything other than change the order of the columns. 
      
      #save PFT file
      write.table(PFTfile[,-ncol(PFTfile)], paste(unlist(strsplit(PFTfileName,".txt")), MC, ".txt", sep=""), row.names=F, quote=F, sep="\t")
      
      # copy necessary files to Model-files folder
      path <- modelpath
      copy <- file.copy(paste(unlist(strsplit(PFTfileName,".txt")), MC, ".txt", sep=""),  path)
      
      copy <- file.copy("Input-files/AppRate.txt",  path) ### There is no need for this correct?????
      
      # change directory
      setwd(modelpath)
      
      # cmd call for running IBC
      mycall<-paste('wine ./IBCgrassGUI.exe', ModelVersion, GridSize, Tmax, InitDuration, paste("./",unlist(strsplit(PFTfileName,".txt")), MC, ".txt", sep=""), 
                    SeedInput, belowres, abres, abampl, tramp, graz, cut, 
                    week_start, HerbDuration, 1, EffectModel, scenario, MC, sep=" ")
      
      # start treatment run
      system(mycall, intern=TRUE)
      
      # change working directory
      setwd('..')
      
      # remove old files
      #remove <- file.remove(paste(modelpath, unlist(strsplit(PFTfileName,".txt")), MC, ".txt", sep=""))
      file.remove(paste(modelpath, unlist(strsplit(PFTfileName,".txt")), MC, ".txt", sep=""))
      #remove <- file.remove("Model-files/HerbFact.txt")
      #remove <- file.remove("Model-files/AppRate.txt")
    }
  
  # stop cluster for parallel processing
  stopCluster(cl)
  
  #####
  # copy treatment
  #####
  ## dir.create("currentSimulation/1", recursive=TRUE)
  setwd(modelpath)
  setwd("..")
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
  #dir.create(paste(sim_dir,sim_subpath,"/HerbicideSettings", sep=""), recursive=TRUE)
  file_list <- list.files(pattern=paste(unlist(strsplit(PFTfileName,".txt")),"*",sep=""))
  file_list <- file_list[file_list!=PFTfileName] #TODO: stimmt das mit der Datei???
  #file_list <- c(file_list, "SimulationSettings.Rdata")
  copy <- file.copy(file_list ,  paste(sim_dir,sim_subpath,"/HerbicideSettings", sep=""))
  #  todo make sure that all files were copied!
  if (all(copy==T)) file.remove(file_list)
  copy <- file.copy(paste0(modelpath,"HerbFact.txt") ,  paste(sim_dir,sim_subpath,"/HerbicideSettings", sep=""))
  if (all(copy==T)) file.remove(paste0(modelpath,"HerbFact.txt"))
  remove <- file.remove(paste0(modelpath,"AppRate.txt"))
  # copy <- file.copy(paste0(modelpath,"PFTsensitivity.txt") ,  paste(sim_dir,sim_subpath,"/HerbicideSettings", sep=""))
  # if (all(copy==T)) file.remove("PFTsensitivity.txt")
}