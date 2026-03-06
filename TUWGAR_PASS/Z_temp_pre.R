Z_temp_pre <- function(GAR, L_cat){
  
  #GAR: (sf) - $ZHYD = ID HUs; $NextDown = HUs immediatly downstream;
  #L_cat: (dataframe) - [,1] HU ID; [,2] level of each HU;
  
  cat("Start input format conversion ---> It may require A LOT of time \n")
  
  dummy <- vector()
  for(i in 1:length(GAR$geometry)){
    dummy[i] <- which(L_cat[,1] == GAR$ZHYD[i])
  }
  L_cat[,2] <- L_cat[,2][dummy]
  GAR <- cbind(GAR,L_cat[,2])
  colnames(GAR)[dim(GAR)[2]-1] <- "Level"
  
  count <- 1
  count2 <- 1
  steps <- c(10,20,30,40,50,60,70,80,90)
  
  
  Z_temp_pre <- vector('list',length(GAR$ZHYD))
  names(Z_temp_pre) <- GAR$ZHYD
  for(l_cat in 1:length(unique(GAR$Level))){
    current_cat <- which(GAR$Level == l_cat)
    
    for(i in current_cat){
      
      if(is.na(steps[count]) == FALSE){
        if(round((count2/length(GAR$geometry))*100) != 0){
          if((round((count2/length(GAR$geometry))*100,0) %% steps[count]) == 0){
            cat("~", steps[count], " Done \n", sep="")
            count <- count + 1
          }
        }
      }
      
      Z_temp <- GAR$ZHYD[i]
      flag <- FALSE
      temp <- Z_temp
      while(flag == FALSE){
        temp <- GAR$ZHYD[which(GAR$NextDown %in% temp)]
        if(length(temp) != 0){
          Z_temp <- c(Z_temp, temp)
        }
        if(length(temp) == 0){
          flag <- TRUE
        }
      }
      
      flag <- vector()
      for(ff in 1:length(Z_temp)){
        flag[ff] <- which(GAR$ZHYD == Z_temp[ff])
      }
      
      down <- vector()
      for(ff in 1:length(flag)){
        down[ff] <- GAR$NextDown[flag[ff]]
      }
      Z_temp <- as.data.frame(cbind(Z_temp, (1:length(down))))
      Z_temp[,2] <- as.numeric(Z_temp[,2])
      
      for(ff in 1:dim(Z_temp)[1]){
        if((down[ff] %in% Z_temp[,1]) == FALSE){
          down[ff] <- NA
        }
        if((down[ff] %in% Z_temp[,1]) == TRUE){
          down[ff] <- Z_temp[which(Z_temp[,1] == down[ff]),2]
        }
      }
      down <- as.numeric(down)
      
      ##AREA OF THE SUB-BASINS FOUND
      area_cats <- vector()
      for(yy in 1:length(flag)){
        area_cats[yy] <- as.numeric(st_area(GAR$geometry[flag[yy]]))
      }
      
      Level <- rep(0, length(down))
      L <- unique(c(1:length(down))[!c(1:length(down)) %in% down])
      Level[L] <- 1
      flag <- 2
      k <- FALSE
      
      while(k==FALSE){
        L <- unique(down[L])
        Level[L] <- flag
        flag <- flag + 1
        if(sum(!is.na(L))==0){
          k <- TRUE
        }
      }
      
      Z_temp <- cbind(Z_temp, Level, down)
      Z_temp_pre[[i]] <- Z_temp
      count2 <- count2 + 1
    }
  }
  return(Z_temp_pre)
  
}