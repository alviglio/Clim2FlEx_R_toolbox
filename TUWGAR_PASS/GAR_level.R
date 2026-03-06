GAR_level <- function(GAR){
  
  #GAR: (sf) - $ZHYD = ID HUs; $NextDown = HUs immediatly downstream;
  
  level <- vector()
  for(i in 1:dim(GAR)[1]){
    temp <- GAR$ZHYD[i]
    flag <- FALSE
    count <- 1
    while(flag == FALSE){
      temp <- GAR$ZHYD[which(GAR$NextDown %in% temp)]
      if(length(temp) == 0){
        flag <- TRUE
      }
      if(length(temp) != 0){
        count <- count + 1
      }
    }
    level[i] <- count
  }
  L_cat <- as.data.frame(cbind(GAR$ZHYD,level))
  L_cat[,2] <- as.numeric(L_cat[,2])
  colnames(L_cat) <- c("ZHYD","Level")
  return(L_cat)
}