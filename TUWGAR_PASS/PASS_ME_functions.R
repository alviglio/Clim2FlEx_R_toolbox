ME.TUWmodel <- function(param, GAR, CI, DISC, topology, Z_temp_pre, me.fn, COR){
  
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
  
  ME <- vector()
  for(i in 1:length(DISC)){
    simu <- G[,which(colnames(G) == names(DISC)[i])]
    simu <- window(simu, start=as.Date("2010-01-01"), end = as.Date("2019-12-31"))
    obs_c <- DISC[[i]]
    dummy <- window(CI[[1]][,1], start=as.Date("2010-01-01"), end=as.Date("2019-12-31"))
    obs_c <- window(obs_c, start=as.Date("2010-01-01"), end=as.Date("2019-12-31"))
    obs_c <- apply(merge.zoo(obs_c,dummy),1,sum) - dummy
    
    if(me.fn == "KGE"){
      ME[i] <- KGE(simu,window(obs_c,start=as.Date("2010-01-01"), end = as.Date("2019-12-31")))
    }
    if(me.fn == "wNSE"){
      ME[i] <- wNSE(simu,window(obs_c,start=as.Date("2010-01-01"), end = as.Date("2019-12-31")))
    }
  }
  
  cat("----------- End Regional Simulation: Routing Module -----------", " \n", sep="")
  timestamp()
  return(ME)
}


ME.TUWmodel2 <- function(param, CI, DISC, topology, Z_temp_pre, me.fn, COR){
  
  timestamp()
  cat("----------- Start Local Simulation -----------", " \n", sep="")
  
  #n_cores <- detectCores()
  #cluster <- makeCluster(n_cores -2)
  cluster <- makeCluster(COR)
  registerDoParallel(cluster) 
  
  ME <- list()
  ME <- foreach(i = 1:dim(param)[1], .packages = c("TUWmodel", "zoo","hydroGOF"), .combine = "c") %dopar% {
    
    NOME <- names(DISC)[i]
    obs_c <- DISC[[i]]
    dummy <- window(CI[[1]][,1], start=as.Date("2000-01-01"), end=as.Date("2020-12-31"))
    obs_c <- window(obs_c, start=as.Date("2000-01-01"), end=as.Date("2020-12-31"))
    obs_c <- apply(merge.zoo(obs_c,dummy),1,sum) - dummy
    
    Z_temp <- Z_temp_pre[[as.numeric(NOME)]]
    
    par_un <- as.numeric(param[which(rownames(param) == NOME),])
    
    inter <- matrix(NA, ncol = dim(Z_temp)[1], nrow = length(CI[[1]][,1]))
    UH <- vector()
    
    t <- seq(1:5000)
    UH <- (1/((par_un[16])*factorial(0)))*(t^(0))*exp(-(t/par_un[16]))
    UH <- UH[1:1000]
    UH <- UH/sum(UH)
    
    for(pp in 1:max(Z_temp$Level)){
      qualicat <- which(Z_temp$Level == pp)
      for(jj in 1:length(qualicat)){
        
        topo <- topology[[which(as.numeric(names(topology)) == Z_temp$Z_temp[qualicat[jj]])]]
        TUW <- TUWmodel(prec = as.matrix(CI[[1]][,which(colnames(CI[[1]]) %in% topo$pxl)]),
                        airt = as.matrix(CI[[2]][,which(colnames(CI[[2]]) %in% topo$pxl)]),
                        ep = as.matrix(CI[[3]][,which(colnames(CI[[3]]) %in% topo$pxl)]),
                        area = topo$weights, param = par_un[1:15])
        
        FC <- topo$zones*(10^-3)/86400
        
        flag_dummy <- FALSE           
        if(dim(TUW$qzones)[1] == 1){      
          TUW <- t(TUW$qzones)*FC
          flag_dummy <- TRUE
        }
        
        if(flag_dummy == FALSE){     
          conv <- c(1:dim(TUW$qzones)[2])
          TUW$qzones[,conv] <- TUW$qzones[,conv] %*% diag(round(FC[conv],4))
          TUW <- apply(TUW$qzones,1,sum)
        }
        
        if(Z_temp$Level[qualicat[jj]] != 1){
          dummy <- which(Z_temp$down == qualicat[jj])
          if(is.null(dim(inter[,dummy])) == TRUE){
            add <- inter[,dummy]
            add[is.na(add)] <- 0
            TUW <- TUW + add
          }
          if(is.null(dim(inter[,dummy])) == FALSE){
            add <- apply(inter[,dummy], 1, sum, na.rm=TRUE)
            TUW <- TUW + add
          }
        }
        
        TUW[is.na(TUW)] <- 0
        Q_out <- convolve(TUW, rev(UH), type='open')
        Q_out <- Q_out[1:length(TUW)]
        
        inter[,qualicat[jj]] <- Q_out
        
        TUW[is.na(TUW)] <- 0
        Q_out <- convolve(TUW, rev(UH), type='open')
        Q_out <- Q_out[1:length(TUW)]
        
        inter[,qualicat[jj]] <- Q_out
      }
    }
    
    simu <- zoo(Q_out, time(CI[[1]]))
    simu <- window(simu, start=as.Date("2010-01-01"), end=as.Date("2019-12-31"))
    
    
    if(me.fn == "KGE"){
      ME[[i]] <- KGE(simu, window(obs_c, start=as.Date("2010-01-01"), end=as.Date("2019-12-31")))
    }
    if(me.fn == "wNSE"){
      ME[[i]] <- wNSE(simu, window(obs_c, start=as.Date("2010-01-01"), end=as.Date("2019-12-31")))
    }
  }
  cat("----------- End Local Simulation -----------", " \n", sep="")
  stopCluster(cluster)
  timestamp()
  return(ME)
}