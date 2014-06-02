#' fn.metaSIM
#' 
#' @title A metacommunity simulation for ecologists
#' 
#' @description This function is a lottery-based, zero-sum, spatially explicit 
#' simulation that can include neutral and/or niche-based dynamics.  Results are 
#' written to a .csv file in a SIM_OUTPUT directory.  Output includes parameter 
#' settings, diversity partitioning outcomes, and variation partitioning outcomes.
#' 
#' @usage 
#' fn.metaSIM(scenario.ID = NA, alpha.fisher = 2, nu = 1e-04, 
#' speciation.limit = NA, n.timestep = 10, landscape.edge.size = 5, 
#' ave.JL = 200, JL.min.proportion = 1, target.m = 0.25, IL.scale = NA, 
#' JL.scale = NA, Ef.scale = NA, Ef.specificity = 0, 
#' Ef.specificity.scale = NA, SWM.slope = 0, trait.dispersal.median = 1, 
#' trait.dispersal.range = 0, trait.Ef.sd = 0.1, q.order = c(0, 2), 
#' save.sim = TRUE, output.dir.path = "SIM_OUTPUT", 
#' keep.timesteps = c(1, 10))
#' 
#' @param scenario.ID A character string used to identify a simulation.  Used in the 
#' simulation data file name (\code{SIM_{scenario.ID}_yyyymmdd_hhmmss_{6digit_random_number}.rda}) 
#' and appended to the end of the metadata file name (\code{sim.metadata_{scenario.ID}.csv}).
#' @param alpha.fisher Numeric, Fisher's alpha used to calculate the regional rank 
#' abundance curve that is used to initiate the metacommunity simulation.  Must be > 0. 
#' @param nu Numeric, Hubbell's "speciation rate", but can be interpreted as the probability 
#' of the appearance of a novel species.  If set to 0, no novel species will appear during 
#' the simulation and regional extinctions will be permanent.
#' @param speciation.limit Numeric, a limit to the number of novel species that can appear 
#' in the simulation.  It is necessary to set a limit for simulations with a large initial 
#' richness.  
#' @param n.timestep Numeric, number of generations in the simulation.  Must be a positive 
#' integer.  
#' @param landscape.edge.size Numeric, sets the edge size for the simulation "landscape".  
#' For example, there will be 16 "sites" in a simulation in which 
#' \code{landscape.edge.size = 4}.  Must be a positive integer > 1.
#' @param ave.JL Numeric, sets the average assemblage size.  Must be a positive integer.  
#' @param JL.min.proportion Numeric, sets the minimum assemblage size that can occur in 
#' the simulation as a proportion of the average \code{JL}.  Must be a value that is > 0 
#' and <= 1.  A value of 1 will result in all assemblages having a size of \code{ave.JL}.
#' @param target.m Numeric.  \code{m} determines the relative contribution of local and 
#' regional pools to the local recruitment pool (Hubbell 2001, Etienne 2005).  Must be a 
#' value >=0 and < 1.  \code{m = IL / (IL + JL - 1)} where \code{IL} is the size of the 
#' immigrant pool and \code{JL} is the size of the extant local assemblage.  Note that 
#' if \code{m = 1}, then \code{IL} is \code{Inf}.
#' @param IL.scale Numeric.  Scale of spatial heterogeneity (among sites) in immigrant 
#' pool size (\code{IL}).  See details below.  
#' @param JL.scale Numeric.  Spatial heterogeneity (among sites) in assemblage size 
#' (\code{JL}).  Same as \code{IL.scale}.  See details below.
#' @param Ef.scale Numeric.  Spatial scale of heterogeneity (among sites) in the 
#' environmental filter (\code{Ef}).  See details below.  
#' @param Ef.specificity Numeric.  Determines how broad or narrow a trait set will be 
#' selected by an environmental filter at a site.  A value of 0 will result in highly 
#' selective environmental filters.  Values >> 0 will result in a simulation resembling 
#' a neutral model.  An assigned value < 0 results in a special case that creates a 
#' convergent filter in half of the sites in the metacommunity.  
#' @param Ef.specificity.scale Numeric, must be -1, 1, or NA.  Describes the spatial 
#' heterogeneity in \code{Ef.specificity}.  A value of -1 is high spatial asynchrony, a 
#' value of 1 is high spatial synchrony, and a value of NA will assign the same 
#' \code{Ef.specificity} to all sites.  
#' @param SWM.slope Numeric, must be >= 0.  The slope of the dispersal kernel (this sets 
#' the value for \code{w} in Eq. 6 in Gravel et al. 2006).  A value of 0 will result in 
#' no dispersal limitation, and thus, a spatially implicit simulation.  Larger values 
#' increase dispersal limitation.  
#' @param trait.dispersal.median Numeric in the range [0,1].  Sets the relative probabilities
#'  that each taxonomic group will emigrate (enter the regional species pool and contribute 
#'  to IL).   
#' @param trait.dispersal.range Numeric in the range [0,1].  A value of 0 will result 
#' in equivalent emigration among all taxonomic groups in the simulation.  Values greater 
#' than 0 create interspecific variation in emigration.  
#' @param trait.Ef.sd Numeric, >= 0.  Niche breadth for all species in the simulation 
#' (i.e., the value of sigma in Fig 1. in Gravel et al. 2006).  Large values result in 
#' a simulation resembling neutral community models (Hubbell 2001).  
#' @param q.order Numeric, can be scalar or a vector of integers >= 0.  This parameter 
#' sets the order q used to calculate diversity metrics.  q = 0 calculates diversity 
#' based on presence/absence (biased toward rare species), q = 1 is unbiased (i.e., Shannon 
#' diversity), and q = 2 is biased to reflect differences in the composition of dominant 
#' species (i.e., Simpson diversity) (Jost 2007, Chao et al. 2012).  See 
#' \link[vegetarian]{d} in the \pkg{vegetarian} package.
#' @param save.sim Boolean (\code{TRUE/FALSE}), should the simulation be saved?  Metadata 
#' will be saved, but you may not want to save the entire simulation if it is very large, 
#' or if you are running many simulations.
#' @param output.dir.path Character string naming the output directory.  Default is 
#' \code{SIM_OUTPUT}.  The simulation will create the directory if it doesn't already exist.  
#' @param keep.timesteps Must be an integer, can be a vector of integers.  Determines which 
#' time steps to include in the metadata file.  
#' 
#' @references 
#' Chao, A., C. H. Chiu, and T. C. Hsieh. 2012. Proposing a resolution to debates on diversity 
#' partitioning. Ecology 93:2037--2051.
#' 
#' Etienne, R. S. 2005. A new sampling formula for neutral biodiversity. Ecology Letters 8:253--260.
#' 
#' Gravel, D., C. D. Canham, M. Beaudet, and C. Messier. 2006. Reconciling niche and neutrality: 
#' the continuum hypothesis. Ecology Letters 9:399--409.
#' 
#' Hubbell, S. P. 2001. A unified theory of biodiversity and biogeography. Princeton University 
#' Press.
#' 
#' Jost, L. 2007. Partitioning diversity into independent alpha and beta components. 
#' Ecology 88:2427--2439.
#' 
#' @export
#' 
fn.metaSIM<-function(
  scenario.ID = NA, 
  alpha.fisher = 2, 
  nu = 1e-04, 
  speciation.limit = NA, 
  n.timestep = 10, 
  landscape.edge.size = 5, 
  ave.JL = 200, 
  JL.min.proportion = 1,         
  target.m = 0.25, 
  IL.scale = NA, 
  JL.scale = NA, 
  Ef.scale = NA, 
  Ef.specificity = 0,
  Ef.specificity.scale = NA,
  SWM.slope = 0, 
  trait.dispersal.median = 1, 
  trait.dispersal.range = 0, 
  trait.Ef.sd = 0.1, q.order = c(0,2), 
  save.sim = TRUE, 
  output.dir.path = "SIM_OUTPUT", 
  keep.timesteps = c(1, 10)
  ){
  n.sites <- landscape.edge.size^2
  JM <- n.sites * ave.JL
  JL.min <- JL.min.proportion * ave.JL
  JL.min <- ifelse(JL.min < 1, 1, JL.min)
  IL.intensity <- round((target.m * (ave.JL - 1))/(1 - target.m), 
                        0)
  landscape <- fn.make.landscape(JM = JM, n.sites = n.sites, 
                                 JL.scale = JL.scale, 
                                  Ef.scale = Ef.scale, 
                                  Ef.specificity = Ef.specificity,
                                  Ef.specificity.scale=Ef.specificity.scale,
                                 IL.scale = IL.scale, IL.intensity = IL.intensity, JL.min)
  dat.gamma.t0 <- fn.set.regional.species.pool(n.timestep = n.timestep, 
                                               nu = nu, speciation.limit = speciation.limit, JM = JM, 
                                               alpha.fisher = alpha.fisher, trait.dispersal.median = trait.dispersal.median, 
                                               trait.dispersal.range = trait.dispersal.range)
  taxa.list <- as.character(dat.gamma.t0$taxa.list)
  d.temp <- expand.grid(Ef = landscape$dat$Ef, stringsAsFactors = FALSE, 
                        trait.Ef = dat.gamma.t0$trait.Ef)
  d.temp$Ef.specificity <- landscape$dat$Ef.specificity
  lambda.Ef.siteBYspp <- matrix(data = mapply(FUN = fn.lambda, 
                                              trait.optimum = d.temp$trait.Ef, Ef = d.temp$Ef, Ef.specificity = d.temp$Ef.specificity, 
                                              MoreArgs = list(niche.breadth = trait.Ef.sd)), nrow = nrow(landscape$dat), 
                                ncol = length(taxa.list), byrow = FALSE)
  lambda.Ef.siteBYspp <- lambda.Ef.siteBYspp/rowSums(lambda.Ef.siteBYspp)
  n.sites <- nrow(landscape$dat)
  R.probs.t0 <- lambda.Ef.siteBYspp * (rep(1, n.sites) %o% 
                                         dat.gamma.t0$regional.RA)
  R.probs.t0 <- R.probs.t0/rowSums(R.probs.t0)
  R.probs.list <- as.list(data.frame(t(R.probs.t0)))
  J.t0 <- data.frame(row.names = NULL, t(mapply(FUN = fn.lottery.recruit, 
                                                vect.recruitment.weights = R.probs.list, scalar.JL = as.list(landscape$dat$JL), 
                                                MoreArgs = list(vect.taxa.list = taxa.list))))
  J <- list()
  J[[1]] <- J.t0
  J.t.minus.1 <- J.t0
  for (t.index in 2:n.timestep) {
    J.t <- fn.recruit.Jt(landscape.xy = landscape$dat[, c("x", 
                                                          "y")], nu = nu, SWM.slope = SWM.slope, J.t.minus.1 = J.t.minus.1, 
                         taxa.list = taxa.list, traits.Ef = dat.gamma.t0$trait.Ef, 
                         trait.Ef.sd = trait.Ef.sd, traits.dispersal = dat.gamma.t0$trait.dispersal, 
                         landscape.m = landscape$dat$m, landscape.Ef = landscape$dat$Ef, 
                         landscape.Ef.specificity = landscape$dat$Ef.specificity, 
                         landscape.JL = landscape$dat$JL)
    J[[t.index]] <- J.t
    J.t.minus.1 <- J.t
    print(paste("Timestep:", t.index))
  }
  sim.result.name <- paste("SIM_", scenario.ID, "_", format(Sys.time(), 
                                                            "%Y%m%d_%H%M%S"), "_", trunc(runif(1, 1e+05, 999999)), 
                           sep = "")
  sim.result <- list(scenario.ID = scenario.ID, sim.result.name = sim.result.name, 
                     x = landscape$dat$x, y = landscape$dat$y, JL = landscape$dat$JL, 
                     Ef = landscape$dat$Ef, Ef.specificity = Ef.specificity, 
                     m = landscape$dat$m, alpha.fisher = alpha.fisher, nu.sim = nu, 
                     trait.Ef.sd = trait.Ef.sd, trait.dispersal = dat.gamma.t0$trait.dispersal, 
                     trait.Ef = dat.gamma.t0$trait.Ef, landscape = landscape, 
                     dat.gamma.t0 = dat.gamma.t0, SWM.slope = SWM.slope, J = J, 
                     n.timestep = n.timestep, taxa.list = taxa.list)
  fn.sim.metadata.archive4(sim.result = sim.result, q.order = q.order, 
                           save.sim = save.sim, var.dir = output.dir.path, 
                           keep.timesteps = keep.timesteps)
  return(sim.result)
}