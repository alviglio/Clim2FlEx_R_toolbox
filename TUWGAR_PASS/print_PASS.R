print.PASS <- function (x, ...) {
  cat('Output of the PArameter Set Shuffling algorithm:\n\n')
  cat(' number of catchments:', nrow(x$groups[[1]]$regionalized.parameters$cat.par.pred), '\n')
  cat(' number of model units (e.g. pixels):', nrow(x$groups[[1]]$regionalized.parameters$grd.par.pred), '\n')
  cat('sampling strategy:', x$PASS.options$sampling, '\n')
  cat(' number of loops:', x$PASS.options$maxLoops, '\n')
  cat(' number of groups:', x$PASS.options$nGroups, '\n')
  cat(' number of loops for regional consistency:', x$PASS.options$REGloops, '\n')
  if (x$PASS.options$sampling == 'optim') {
    cat(' proportion of randomized parameters (when optim):', x$PASS.options$optim.subset.cat, '\n')
  }
  cat('lumped regionalization efficiencies:', '\n')
  print(round(apply(sapply(x$groups, function(x) x$cat.eff.lump), 2, summary), 3))
  cat('distributed regionalization efficiencies:', '\n')
  print(round(apply(sapply(x$groups, function(x) x$cat.eff.dist), 2, summary), 3))
  cat(' updated number of train parameter sets:', sum(sapply(x$train.parameters.updated, nrow)), '\n')
}