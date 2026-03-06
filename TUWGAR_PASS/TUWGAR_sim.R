TUWGAR.sim <- function(param, GAR, L_cat, CI, topology, Z_temp_pre, COR){
  
  #M = number of HUs
  #N = number of simlation timesteps
  
  #param: (numeric Mx16) regional parameters for all HUs;
  #GAR: GAR: (sf) - $ZHYD = ID HUs; $NextDown = HUs immediatly downstream;
  #L_cat: (dataframe) - [,1] HU ID; [,2] level of each HU;
  #CI: climate input (list of NxM dataframes, N=timesteps, M=number of HUs), [[1]] daily precipitation, [[2]] mean daily temperature, [[3]]daily poterntial evapotranspiration;
  #topology: list of M dataframes - each dataframe: [,1] pixels touching the HU, [,2] weight, [,3] area in m^2 of the intersection;
  #Z_temp_pre: list of M dataframes - each dataframe: Z_temp_pre[[M]]$Level = level of the HU for the catchment; Z_temp_pre[[M]]$Z_temp = ID of the HU; Z_temp_pre[[M]]$down = ID of downstream HU;
  #COR: (numeric) Number of cores to use in parallelization;
  
  
  dummy <- vector()
  for(i in 1:length(GAR$geometry)){
    dummy[i] <- which(L_cat[,1] == GAR$ZHYD[i])
  }
  L_cat[,2] <- L_cat[,2][dummy]
  GAR <- cbind(GAR,L_cat[,2])
  colnames(GAR)[dim(GAR)[2]-1] <- "Level"
  
  ######## INITIALIZE SIMUL MATRIX #########
  G <- matrix(NA, ncol=length(GAR$geometry), nrow = dim(CI[[1]])[1])
  colnames(G) <- GAR$ZHYD
  timestamp()
  
  cat("----------- Start Regional Simulation: Generative Module -----------", " \n", sep="")
  
  
  #### GENERATION ON ALL HUs ####
  
  current_cat <- sort(GAR$ZHYD, decreasing = FALSE)
  
  cluster <- makeCluster(COR)
  registerDoParallel(cluster)    
  series <- list()
  
  series <- foreach(i = 1:length(current_cat),.packages = c("TUWmodel")) %dopar% {
    
    ii <- current_cat[i]
    topo <- topology[[which(as.numeric(names(topology)) == ii)]]
    TUW <- TUWmodel(prec = as.matrix(CI[[1]][,which(colnames(CI[[1]]) %in% topo$pxl)]),
                    airt = as.matrix(CI[[2]][,which(colnames(CI[[2]]) %in% topo$pxl)]),
                    ep = as.matrix(CI[[3]][,which(colnames(CI[[3]]) %in% topo$pxl)]),
                    area = topo$weights, param = param[which(rownames(param) == ii),1:15])
    
    FC <- topo$zones*(10^-3)/86400
    
    flag_dummy <- FALSE           
    if(dim(TUW$q0)[1] == 1){      
      TUW <- t(TUW$qzones)*FC
      flag_dummy <- TRUE
    }
    if(flag_dummy == FALSE){    
      conv <- c(1:dim(TUW$qzones)[2])
      TUW$qzones[,conv] <- TUW$qzones[,conv] %*% diag(round(FC[conv],4))
      TUW <- apply(TUW$qzones,1,sum)
    }
    
    series[i] <- TUW
    
  }
  stopCluster(cluster)
  
  for(s in 1:length(series)){
    G[,s] <- series[[s]]
  }
  
  cat("----------- End Regional Simulation: Generative Module -----------", " \n", sep="")
  
  #### ROUTING ON ALL HUs ####
  
  cat("----------- Start Regional Simulation: Routing Module -----------", " \n", sep="")
  
  for(l_cat in 1:length(unique(GAR$Level))){
    cat("Level ",l_cat, " out of ",length(unique(GAR$Level)), " \n",sep="")
    current_cat <- which(GAR$Level == l_cat)
    for(i in current_cat){
      Z_temp <- Z_temp_pre[[i]]
      pp <- max(Z_temp$Level)
      qualicat <- which(Z_temp$Level == pp)
      
      
      
      if(Z_temp$Level[qualicat] != 1){
        dummy <- Z_temp$Z_temp[which(Z_temp$down == qualicat)]
        if(is.null(dim(G[,which(colnames(G)%in%dummy)])) == TRUE){
          add <- G[,which(colnames(G)%in%dummy)]
          add[is.na(add)] <- 0
          G[,which(colnames(G) == Z_temp$Z_temp[qualicat])] <- G[,which(colnames(G) == Z_temp$Z_temp[qualicat])] + add
        }
        if(is.null(dim(G[,which(colnames(G)%in%dummy)])) == FALSE){
          add <- apply(G[,which(colnames(G)%in%dummy)], 1, sum, na.rm=TRUE)
          add[is.na(add)] <- 0
          G[,which(colnames(G) == Z_temp$Z_temp[qualicat])] <- G[,which(colnames(G) == Z_temp$Z_temp[qualicat])] + add
        }
      }
      
      ##NASH
      N <- 1
      t <- seq(1:dim(G)[1])
      UH <- (1/((param[which(rownames(param) == Z_temp$Z_temp[qualicat]),16]^N)*factorial(N-1)))*(t^(N-1))*exp(-(t/param[which(rownames(param) == Z_temp$Z_temp[qualicat]),16]))
      
      UH <- UH[1:1000]
      UH <- UH/sum(UH)
      
      G[,which(colnames(G) == Z_temp$Z_temp[qualicat])][is.na(G[,which(colnames(G) == Z_temp$Z_temp[qualicat])])] <- 0
      Q_out <- convolve(G[,which(colnames(G) == Z_temp$Z_temp[qualicat])], rev(UH), type='open')
      Q_out <- Q_out[1:length(G[,which(colnames(G) == Z_temp$Z_temp[qualicat])])]
      
      G[,which(colnames(G) == Z_temp$Z_temp[qualicat])] <- Q_out
    }
  }
  
  
  
  G <- zoo(G, time(CI[[1]][,1]))
  timestamp()
  
  cat("----------- End Regional Simulation: Routing Module -----------", " \n", sep="")
  timestamp()
  return(G)
}