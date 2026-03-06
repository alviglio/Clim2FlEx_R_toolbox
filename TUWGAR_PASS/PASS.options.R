PASS.options <- function(maxLoops=100, nGroups=10, REGloops=5, 
                         generalized.mean.power=-1, proportion.max.eff.update=0.95, 
                         sampling='random', optim.subset.cat=0.7) {
  list(maxLoops=maxLoops, nGroups=nGroups, REGloops=REGloops, 
       generalized.mean.power=generalized.mean.power, proportion.max.eff.update=proportion.max.eff.update,
       sampling=sampling, optim.subset.cat=optim.subset.cat)
}