HU_topology <- function(GAR, grid, CI, L_cat){
  
  #GAR: (sf) - $ZHYD = ID HUs; $NextDown = HUs immediatly downstream;
  #grid: (sf) - [,1] Pixel ID; [,2] Pixel geometry;
  #CI: climate input (list of NxM dataframes, N=timesteps, M=number of HUs), [[1]] daily precipitation, [[2]] mean daily temperature, [[3]]daily poterntial evapotranspiration;
  #L_cat: (dataframe) - [,1] HU ID; [,2] level of each HU;
  
  cat("Start calculation HUs topology - It may require some time")
  
  dummy <- vector()
  for(i in 1:length(GAR$geometry)){
    dummy[i] <- which(L_cat[,1] == GAR$ZHYD[i])
  }
  L_cat[,2] <- L_cat[,2][dummy]
  GAR <- cbind(GAR,L_cat[,2])
  colnames(GAR)[dim(GAR)[2]-1] <- "Level"
  
  
  topology <- vector('list', length(GAR$geometry))
  names(topology) <- GAR$ZHYD
  
  count <- 1
  steps <- c(10,20,30,40,50,60,70,80,90)
  
  for(i in 1:length(GAR$geometry)){
    
    if(is.na(steps[count]) == FALSE){
      if(round((i/length(GAR$geometry))*100) != 0){
        if((round((i/length(GAR$geometry))*100,0) %% steps[count]) == 0){
          cat("~", steps[count], " Done \n", sep="")
          count <- count + 1
        }
      }
    }
    
    area_HU <- as.numeric(st_area(GAR$geometry[i]))
    
    dummy <- st_intersection(grid,GAR$geometry[i])
    pxl <- as.numeric(rownames(dummy))
    weights <- as.numeric(st_area(dummy$geometry))/area_HU
    zones <- as.numeric(st_area(dummy$geometry))
    weights <- weights[pxl %in% colnames(CI[[1]])]
    zones <- zones[pxl %in% colnames(CI[[1]])]
    pxl <- pxl[pxl %in% colnames(CI[[1]])]
    weights <- weights/sum(weights)
    
    dummy <- as.data.frame(cbind(pxl, round(weights,2), zones))
    colnames(dummy) <- c("pxl", "weights", "zones")
    topology[[i]] <- dummy
    
  }
  
  return(topology)
}