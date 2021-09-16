#' @title A metacommunity simulation for ecologists
#' 
#' @aliases fn.metaSIM
#' 
#' @usage 
#' metasim(landscape, ...)
#' 
#' @description metasim() initiates a metacommunity simulation based on a landscape
#' created by the fn.make.landscape() function.
#' 
#' @param landscape Landscape object created by function fn.make.landscape()
#' @param scenario.ID A name for the simulation scenario. All simulations with the same scenario name will have metadata collated in a single .csv file. Default is NA
#' @param sim.ID A name for this particular simulation. Default will provide a random ID.
#' @param alpha.fisher User can use Fisher's alpha to seed a simulation's initial regional source pool.
#' @param nu Probability a novel species will appear in a single recruitment event
#' @param speciation.limit Set's a limit on the number of novel taxa that can appear in a simulation
#' @param JM.limit set's an upper limit for the number of individuals in a simulation
#' @param n.timestep Number of generations in a simulation
#' @param W.r Dispersal kernel slope
#' @param save.sim Binary, will save the simulation as a .rda file if TRUE. Default is FALSE.
#' @param output.dir.path Name of directory to save simulation results and metadata. Default is "SIM_OUTPUT". Simulation will create a sub-directory in the working directory if none exists.
#' @param trait.dispersal Vector of dispersal traits. Default is NULL
#' @param trait.dispersal.median Scalar value for dispersal applied to all species if no vector is provided. Default is 1.
#' @param trait.dispersal.range Range of dispersal values given to species if dispersal traits are randomly assigned.
#' @param trait.Ef Vector of species' niche positions
#' @param trait.Ef.sd Vector of species' niche breadths
#' @param gamma.abund Vector of regional abundances, can be used to seed a simulation
#' @param J.t0 Site by species data.table or matrix that can be used to seed a simulation
#' @param taxa.list A character vector of names for species
#' 
#' @export
#' 
#' @details There are two steps to creating a metacommunity simulation in MCSim:
#' 1. Make a "landscape" -- The landscape is the “game board” on which the simulation plays out, and it is created using the fn.make.landscape function.
#' 2. Run the simulation -- Once the landscape is created, you can pass the landscape object to fn.metaSIM along with parameter settings that define the rules for how metacommunity dynamics will play out in the metacommunity simulation. Note that the current version of MCSim is zero sum, which means there will always be JM individuals in the simulation during each generation.
#' For a tutorial, see \url{http://rpubs.com/sokole/159425}
#' 
#' Note that a user can choose to save simulation output to a directory set by output.dir.path
#' by setting save.sim = TRUE
#' 
metasim <- function (
  # -- unchanged vars
  landscape = NA, 
  scenario.ID = NA,
  sim.ID = NA,
  alpha.fisher = 1, 
  nu = 1e-04, 
  speciation.limit = 0, 
  JM.limit = 1e+05, 
  n.timestep = 10, 
  W.r = 0, 
  save.sim = FALSE, 
  output.dir.path = "SIM_OUTPUT",
  
  # -- New or edited vars
  trait.dispersal = NULL, # a vector of dispersal traits
  trait.dispersal.median = 1,# or set a median and range of dispersals
  trait.dispersal.range = 0, 
  
  trait.Ef = NULL,
  trait.Ef.sd = NULL, #NULL creates a nearly neutral simulation
  
  gamma.abund = NULL,
  J.t0 = NULL,
  taxa.list = NULL,
  ...) {
  if (class(landscape) != "list") {
    print("This simulation needs a landscape!")
  } else {
    try({
      JM <- sum(landscape$site.info$JL)
      n.sites <- length(landscape$site.info$JL)
      
      if(!is.null(J.t0)){
        gamma.abund<-colMeans(J.t0)
        
        J.t0.RAs<-J.t0/rowSums(J.t0)
        J.t0 <- round(J.t0.RAs * landscape$site.info$JL, 0)
        
        if(!is.null(taxa.list)){
          # -- if taxa.list is provided, use it to name species
          try({names(J.t0) <- taxa.list})
        }else if(!is.numeric(names(J.t0))){
          # -- if no taxa.list, use names from J.t0
          try({taxa.list <- names(J.t0)})
        }else{
          # -- if names in J.t0 are numierc, add 'spp' as a prefix
          names(J.t0) <- paste0('spp',1:ncol(J.t0))
          taxa.list <- names(J.t0)
        }
      }
      
      dat.gamma.t0 <- fn.set.regional.species.pool(n.timestep = n.timestep, 
                                                   nu = nu, 
                                                   speciation.limit = speciation.limit, 
                                                   JM = ifelse(JM > JM.limit, JM.limit, JM), 
                                                   alpha.fisher = alpha.fisher, 
                                                   trait.dispersal = trait.dispersal, # a vector of dispersal traits
                                                   trait.dispersal.median = trait.dispersal.median,# or set a median and range of dispersals
                                                   trait.dispersal.range = trait.dispersal.range, 
                                                   
                                                   gamma.abund = gamma.abund,
                                                   
                                                   Ef.min = min(landscape$site.info$Ef)-abs(0.01*min(landscape$site.info$Ef)), 
                                                   Ef.max = max(landscape$site.info$Ef)+abs(0.01*max(landscape$site.info$Ef)),
                                                   
                                                   trait.Ef = trait.Ef,
                                                   trait.Ef.sd = trait.Ef.sd,
                                                   taxa.list = taxa.list)
      
      # -- extend J.t0 to include possible new species 
      taxa.list <- as.character(dat.gamma.t0$taxa.list)
      if(!is.null(J.t0)){
        if(length(taxa.list)>ncol(J.t0)){
          J.t0[, c( (ncol(J.t0)+1):length(taxa.list) )]<-0
          names(J.t0)<-taxa.list
        }
      }
      
      # -- combine traits and haibtats in a long data table
      d.temp <- data.frame(
        site.id = c(1:length(landscape$site.info$Ef.specificity)),
        spp.id = rep(dat.gamma.t0$taxa.list, each = length(landscape$site.info$Ef)),
        expand.grid(Ef = landscape$site.info$Ef, 
                    stringsAsFactors = FALSE, 
                    trait.Ef = dat.gamma.t0$trait.Ef),
        niche.breadth = rep(dat.gamma.t0$trait.Ef.sd, each = length(landscape$site.info$Ef)),
        Ef.specificity = landscape$site.info$Ef.specificity #from landscape
      )
      
      # -- estimate lottery weights based on habitat preferences, niche breadths, and site habitat chcaracteristics
      lambda.Ef.siteBYspp <- matrix(data = mapply(FUN = fn.lambda, 
                                                  trait.optimum = d.temp$trait.Ef, 
                                                  Ef = d.temp$Ef, 
                                                  Ef.specificity = d.temp$Ef.specificity,
                                                  niche.breadth = d.temp$niche.breadth),
                                    nrow = n.sites, 
                                    ncol = length(taxa.list), 
                                    byrow = FALSE)
      
      lambda.Ef.siteBYspp <- lambda.Ef.siteBYspp/rowSums(lambda.Ef.siteBYspp)
      R.probs.t0 <- lambda.Ef.siteBYspp * (rep(1, n.sites) %o% 
                                             dat.gamma.t0$regional.RA)
      R.probs.t0 <- R.probs.t0/rowSums(R.probs.t0)
      R.probs.list <- as.list(data.frame(t(R.probs.t0)))
      
      # ---------------------------------
      # -- lottery to initiate simulation at time t0
      # -------------
      if(is.null(J.t0)){
        J.t0 <- data.frame(row.names = NULL, t(mapply(FUN = fn.lottery.recruit, 
                                                      vect.recruitment.weights = R.probs.list, 
                                                      scalar.JL = as.list(landscape$site.info$JL), 
                                                      MoreArgs = list(vect.taxa.list = taxa.list))))
      }
      
      # ---------------------------------
      # -- initialize J list for time 0
      # -------------
      J <- list()
      J[[1]] <- J.t0
      J.t.minus.1 <- J.t0
      
      # -- make J a long matrix
      suppressMessages({
        J.long<-data.frame(
          timestep=1,
          reshape2::melt(as.matrix(J.t0),
                         value.name = 'count'))
      })
      names(J.long)<-c('timestep','site','spp','count')
      
      # ---------------------------------
      # -- loterry recruits for all subsequent time steps
      # -------------
      for (t.index in 2:n.timestep) {
        J.t <- fn.recruit.Jt(mat.geodist = landscape$dist.mat, 
                             nu = nu, 
                             SWM.slope = W.r, 
                             J.t.minus.1 = J.t.minus.1, 
                             taxa.list = taxa.list, 
                             traits.Ef = dat.gamma.t0$trait.Ef, 
                             trait.Ef.sd = dat.gamma.t0$trait.Ef.sd, 
                             traits.dispersal = dat.gamma.t0$trait.dispersal, 
                             m = landscape$site.info$m, 
                             Ef = landscape$site.info$Ef, 
                             Ef.specificity = landscape$site.info$Ef.specificity, 
                             JL = landscape$site.info$JL)
        
        # -- keep spp counts in long format
        suppressMessages({
          J.long.temp<-data.frame(
            timestep=t.index,
            reshape2::melt(as.matrix(J.t),
                           value.name = 'count')
          )
        })
        names(J.long.temp)<-c('timestep','site','spp','count')
        
        J.long<-rbind(J.long,
                      J.long.temp)
        
        J[[t.index]] <- J.t
        J.t.minus.1 <- J.t
        print(paste("Timestep:", t.index))
      }
      
      # ---------------------------------
      # -- name sim result
      # -------------
      if(is.na(sim.ID)){
        sim.result.name <- paste("SIM_", scenario.ID, "_", 
                                 format(Sys.time(), 
                                        "%Y%m%d_%H%M%S"), 
                                 "_", 
                                 trunc(runif(1, 1e+05, 999999)), 
                                 sep = "")
      }else{
        sim.result.name <- paste(
          "SIM", 
          sim.ID,
          as.character(format(Sys.time(), 
                              "%Y%m%d_%H%M%S")),
          sep='_')
      }
      
      # ---------------------------------
      # -- collate info into sim.result
      # -------------
      sim.result <- list(scenario.ID = scenario.ID, 
                         sim.result.name = sim.result.name, 
                         landscape = landscape, 
                         dat.gamma.t0 = dat.gamma.t0, 
                         W.r = W.r, 
                         J.long = J.long)
      
      sim.result.filename <- paste(output.dir.path, "/", 
                                   sim.result.name, ".rda", sep = "")
      
      # -- summarise site metadata
      siteinfo.metadata<-dplyr::select(landscape$site.info, -site.ID)
      gamma.metadata<-dplyr::select(dat.gamma.t0, -taxa.list)
      sim.result.metadata.wide <- data.frame(
        n.sites = n.sites, 
        n.timestep = n.timestep, 
        alpha.fisher = alpha.fisher, 
        nu.sim = nu, 
        JM = JM, 
        W.r = W.r,
        stringsAsFactors = FALSE)

      metadata.long<-data.frame(
        sim.ID = sim.result.name, 
        rbind(
          data.frame(
            param.name=names(sim.result.metadata.wide),
            param.stat='fixed_value',
            param.stat.val=t(sim.result.metadata.wide),
            row.names=NULL
          ),
          data.frame(
            param.name=names(siteinfo.metadata),
            param.stat='mean',
            param.stat.val=t(dplyr::summarise_all(siteinfo.metadata, dplyr::funs(mean))),
            row.names=NULL
          ),
          data.frame(
            param.name=names(siteinfo.metadata),
            param.stat='sd',
            param.stat.val=t(dplyr::summarise_all(siteinfo.metadata, dplyr::funs(sd))),
            row.names=NULL
          ),
          data.frame(
            param.name=names(gamma.metadata),
            param.stat='mean',
            param.stat.val=t(dplyr::summarise_all(gamma.metadata, dplyr::funs(mean))),
            row.names=NULL
          ),
          data.frame(
            param.name=names(gamma.metadata),
            param.stat='sd',
            param.stat.val=t(dplyr::summarise_all(gamma.metadata, dplyr::funs(sd))),
            row.names=NULL
          ))
      )
      
      # -- check for output directory, create if necessary
      if (!file.exists(output.dir.path)){
        dir.create(output.dir.path)
      } 
      
      # -- save sim if requested
      if (save.sim) 
        save(sim.result, file = sim.result.filename)
      filname.sim.metadata <- paste(output.dir.path, "/sim.metadata_", 
                                    sim.result$scenario.ID, ".csv", sep = "")
      
      # -- save metadata
      if (file.exists(filname.sim.metadata)) {
        dat.sim.metadata <- utils::read.csv(filname.sim.metadata, 
                                     header = TRUE,
                                     stringsAsFactors = FALSE)
        if (!sim.result.name %in% row.names(dat.sim.metadata)) {
          dat.sim.metadata <- rbind(dat.sim.metadata, 
                                    metadata.long)
          utils::write.csv(dat.sim.metadata, 
                    filname.sim.metadata,
                    row.names = FALSE)
        }
      } else {
        utils::write.csv(metadata.long, 
                  filname.sim.metadata,
                  row.names = FALSE)
      }
      try(detach(landscape$site.info), silent = TRUE)
      try(detach(landscape), silent = TRUE)
      return(sim.result)
    })
  }
}

#' assign alias function
#' @export
fn.metaSIM <- metasim