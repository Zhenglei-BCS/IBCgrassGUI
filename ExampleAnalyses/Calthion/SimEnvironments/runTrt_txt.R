

runtrt <- function(ECx=0.05,EP="Establishment",resFolder="SimulationC",PFTfile,MCruns=30,...){
 
  
  # # cluster for parallel processing ==> go outside of the function
  # no_cores <- max(detectCores()-2,1) # you might want to adapt this and give a specific number of cores (e.g. when using a HPC)
  # cl <- makeCluster(no_cores)
  # registerDoParallel(cl)
  
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
  dir.create("currentSimulation/1", recursive=TRUE)
  file_list <- list.files(path = "Model-files/", pattern="Pt__*")
  for (file in file_list){
    path <- "currentSimulation/1"
    copy <- file.copy(paste("Model-files/" ,file , sep="") ,  path)
    #  todo make sure that all files were copied!
    if (copy==T) file.remove(paste("Model-files/" ,file , sep="") )              
  } 
  # GRD files
  file_list <- list.files(path = "Model-files/", pattern="Grd__*")
  for (file in file_list){
    path <- "currentSimulation/1"
    copy <- file.copy(paste("Model-files/" ,file , sep="") ,  path)
    #  todo make sure that all files were copied!
    if (copy==T) file.remove(paste("Model-files/" ,file , sep="") )              
  }
  
  # Copy treamtent settings
  dir.create(paste("currentSimulation/HerbicideSettings", sep=""), recursive=TRUE)
  file_list <- list.files(pattern=paste(unlist(strsplit(PFTfileName,".txt")),"*",sep=""))
  file_list <- file_list[file_list!=PFTfileName] #TODO: stimmt das mit der Datei???
  file_list <- c(file_list, "SimulationSettings.Rdata")
  copy <- file.copy(file_list ,  paste("currentSimulation/HerbicideSettings", sep=""))
  #  todo make sure that all files were copied!
  if (all(copy==T)) file.remove(file_list)
  copy <- file.copy("HerbFact.txt" ,  paste("currentSimulation/HerbicideSettings", sep=""))
  if (all(copy==T)) file.remove("HerbFact.txt")
  copy <- file.copy("PFTsensitivity.txt" ,  paste("currentSimulation/HerbicideSettings", sep=""))
  if (all(copy==T)) file.remove("PFTsensitivity.txt") 
}