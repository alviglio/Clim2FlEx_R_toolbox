PASS <- function (Y,            # list (N.cat) of data.frames (XXX x N.par) with locally lumped calibrated model parameters OR previously obtained PASS output
                  X.cat,        # matrix or data.frame (N.cat x N.dsc) of catchment descriptors
                  X.grd,        # matrix or data.frame (N.grd x N.dsc) of model unit/pixel descriptors
                  wcat,
                  L_cat,        # dataframe (cat and Level)
                  model.eff.fn, # function to calculate the regional efficiency
                  model.eff.fn2,# function to calculate the lumped efficiency
                  lower,        # two vectors specifying scalar real lower and upper bounds on each parameter to be optimized, so that the i-th element
                  upper,        # of 'lower' and 'upper' applies to the i-th parameter.
                  options=PASS.options(), # options for PASS
                  GAR,          # (sf) - $ZHYD = ID HUs; $NextDown = HUs immediatly downstream;
                  CI,           # climate input (list of NxM dataframes, N=timesteps, M=number of HUs), [[1]] daily precipitation, [[2]] mean daily temperature, [[3]]daily poterntial evapotranspiration;
                  DISC,         # list(M catchments) - discharge series (zoo) in m^3/s
                  topology,     # list of M dataframes - each dataframe: [,1] pixels touching the HU, [,2] weight, [,3] area in m^2 of the intersection;
                  Z_temp_pre,   # list of M dataframes - each dataframe: Z_temp_pre[[M]]$Level = level of the HU for the catchment; Z_temp_pre[[M]]$Z_temp = ID of the HU; Z_temp_pre[[M]]$down = ID of downstream HU;
                  COR,          # (numeric) Number of cores to use in parallelization;
                  me.fn,        # (character) model efficiency function. Can be "KGE" for Kling-Gupta or "wNSE" for eighted Nash-Sutcliffe   
                  ...){         # further arguments to be passed 
  #require(rpart)
  # select parameter sets
  N.cat <- nrow(X.cat)  # number of CALIBRATED sub-catchments
  N.grd <- nrow(X.grd)  # number of TOTAL sub-catchments
  if (class(Y) != 'PASS') {
    
    train.parameters <- Y
    options <- do.call(PASS.options, as.list(options))
    overall.eff.threshold <- rep(-999, options$nGroups)
    PASS.group <- list()
    PASS.save <- list()
  } else {
    previousPASSout <- Y
    overall.eff.threshold <- sapply(previousPASSout$groups, function (x) x$overall.eff)
    options <- do.call(PASS.options, as.list(options))
    options$nGroups <- length(previousPASSout$groups)
    train.parameters <- previousPASSout$train.parameters.updated
    PASS.group <- list()
    PASS.save <- previousPASSout
  }
  
  # parameter normalisation (to calculate distances in the regional consistency algorithm)
  train.parameters.norm <- train.parameters
  for (i.cat in 1:N.cat) {
    dummy_max <- t(upper %*% t(rep(1, nrow(train.parameters[[i.cat]]))))
    dummy_min <- t(lower %*% t(rep(1, nrow(train.parameters[[i.cat]]))))
    train.parameters.norm[[i.cat]][,-1] <- (train.parameters[[i.cat]][,-1] - dummy_min)/(dummy_max - dummy_min)  # transformation  (the first column contains efficiencies)
  }
  dummy_max.cat <- t(upper %*% t(rep(1, N.cat)))
  dummy_min.cat <- t(lower %*% t(rep(1, N.cat)))
  dummy_max.grd <- t(upper %*% t(rep(1, N.grd)))
  dummy_min.grd <- t(lower %*% t(rep(1, N.grd)))
  
  cat('%%%% ------------- START ------------- %%%%\n')
  time0 <- Sys.time()
  for (i.PASS in 1:options$maxLoops) {
    #if (i.PASS %% round(options$maxLoops/10) == 0) cat('     Loop', i.PASS, 'out of', options$maxLoops, 'loops\n')
    cat(i.PASS)
    i.Group <- sample(options$nGroups, 1)
    if (options$sampling == 'random') {
      selected.parameters <- t(sapply(train.parameters, function (x) x[sample(nrow(x), 1),]))   # matrix with N rows (n. of catchments) and M columns (n. of parameters)
      selected.parameters <- as.matrix(selected.parameters)
    } else if (options$sampling == 'optim') {
      selected.parameters <- previousPASSout$groups[[i.Group]]$selected.parameters
      select.cat <- sample(N.cat, round(N.cat*options$optim.subset.cat))
      selected.parameters[select.cat,] <- t(sapply(train.parameters[select.cat], function (x) x[sample(nrow(x), 1),]))
    }
    # parameter normalisation (to calculate distances in the regional consistency algorithm)
    selected.parameters[,-1] <- (selected.parameters[,-1] - dummy_min.cat)/(dummy_max.cat - dummy_min.cat)  # transformation
    # Regional consistency algorithm:
    for (i.REG in 1:options$REGloops) {
      # regionalization
      regionalization <- .DT.app(par.in=selected.parameters[, -1], catch_CDs=X.cat, grds_CDs=X.grd)  # output = cat.par.pred and grd.par.pred
      for (i.cat in 1:N.cat) {
        par.reg.cat <- regionalization$cat.par.pred[i.cat,]
        distances <- (train.parameters.norm[[i.cat]][,-1] - t(par.reg.cat %*% t(rep(1, nrow(train.parameters.norm[[i.cat]])))))^2
        selected.parameters[i.cat,] <- train.parameters.norm[[i.cat]][which.min(apply(distances, 1, sum, na.rm=TRUE)),]
      }
    }
    colnames(regionalization$cat.par.pred) <- colnames(selected.parameters)[-1]
    rownames(regionalization$cat.par.pred) <- rownames(X.cat)
    colnames(regionalization$grd.par.pred) <- colnames(selected.parameters)[-1]
    rownames(regionalization$grd.par.pred) <- rownames(X.grd)
    regionalization$cat.par.pred[regionalization$cat.par.pred > 1] <- 1
    regionalization$cat.par.pred[regionalization$cat.par.pred < 0] <- 0
    regionalization$grd.par.pred[regionalization$grd.par.pred > 1] <- 1
    regionalization$grd.par.pred[regionalization$grd.par.pred < 0] <- 0
    regionalization$cat.par.pred <- regionalization$cat.par.pred*(dummy_max.cat - dummy_min.cat) + dummy_min.cat
    regionalization$grd.par.pred <- regionalization$grd.par.pred*(dummy_max.grd - dummy_min.grd) + dummy_min.grd
    eff.lump <- rep(NA, N.cat)
    eff.dist <- rep(NA, N.cat)
    
    
    eff.lump <- model.eff.fn2(param = regionalization$cat.par.pred,CI,DISC,topology,Z_temp_pre, me.fn = me.fn, COR) 
    
    for (i.cat in 1:length(DISC)){
      qc <- which(names(train.parameters) == names(DISC)[i.cat])
      if (eff.lump[i.cat] > options$proportion.max.eff.update*max(train.parameters[[qc]][, 1])) {    
        train.parameters[[qc]] <- rbind(train.parameters[[qc]], c(eff.lump[i.cat], regionalization$cat.par.pred[qc,]))
      }   
    }
    
    param.grd <- regionalization$grd.par.pred
    param.grd[is.na(param.grd)] <- -999
    
    eff.dist <- model.eff.fn(param = param.grd, GAR, CI, DISC, topology, Z_temp_pre, me.fn = me.fn, COR)
    
    overall.eff.dist <- median(eff.dist)
    
    if (options$sampling == 'random') {
      i.Group <- which.min(overall.eff.threshold) 
    }   # in optim stay in the Group you are optimizing
    if (overall.eff.dist > overall.eff.threshold[i.Group]) {
      overall.eff.threshold[i.Group] <- overall.eff.dist
      PASS.group[['overall.eff']] <- overall.eff.dist
      PASS.group[['selected.parameters']] <- selected.parameters
      PASS.group[['regionalized.parameters']] <- regionalization
      PASS.group[['cat.eff.lump']] <- eff.lump
      PASS.group[['cat.eff.dist']] <- eff.dist
      PASS.save[['groups']][[paste('Group', i.Group, sep='')]] <- PASS.group
    }
  }  # end of for (i.PASS in 1:maxLoops)
  PASS.save[['train.parameters.updated']] <- train.parameters
  PASS.save[['PASS.options']] <- options
  cat('%%%% -------------  END  ------------- %%%%\n')
  time1 <- Sys.time()
  print(time1 - time0)
  class(PASS.save) <- 'PASS'
  return(PASS.save)
}










