# ---------------------------------------------------------------------------------------
# -- v0.4.0.9000 dev functions
# ---------------------------------------------------------------------------------------
#' fn.metaSIM
#' 
#' @title A metacommunity simulation for ecologists
#' 
#' @usage 
#' fn.metaSIM(landscape)
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
fn.metaSIM <- function (
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
            param.stat.val=t(dplyr::summarise_each(siteinfo.metadata, dplyr::funs(mean))),
            row.names=NULL
          ),
          data.frame(
            param.name=names(siteinfo.metadata),
            param.stat='sd',
            param.stat.val=t(dplyr::summarise_each(siteinfo.metadata, dplyr::funs(sd))),
            row.names=NULL
          ),
          data.frame(
            param.name=names(gamma.metadata),
            param.stat='mean',
            param.stat.val=t(dplyr::summarise_each(gamma.metadata, dplyr::funs(mean))),
            row.names=NULL
          ),
          data.frame(
            param.name=names(gamma.metadata),
            param.stat='sd',
            param.stat.val=t(dplyr::summarise_each(gamma.metadata, dplyr::funs(sd))),
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
        dat.sim.metadata <- read.csv(filname.sim.metadata, 
                                     header = TRUE,
                                     stringsAsFactors = FALSE)
        if (!sim.result.name %in% row.names(dat.sim.metadata)) {
          dat.sim.metadata <- rbind(dat.sim.metadata, 
                                    metadata.long)
          write.csv(dat.sim.metadata, 
                    filname.sim.metadata,
                    row.names = FALSE)
        }
      } else {
        write.csv(metadata.long, 
                  filname.sim.metadata,
                  row.names = FALSE)
      }
      try(detach(landscape$site.info), silent = TRUE)
      try(detach(landscape), silent = TRUE)
      return(sim.result)
    })
  }
}
# ---------------------------------------------------------------------------------------
#' fn.make.landscape
#' 
#' @title make a simulation landscape
#' 
#' @param site.coords A data.frame of site coordinates. Can be 1, 2, or more dimensions
#' @param dist.mat Alternative to site.coords. Can be a distance matrix or a network map from the igraph package
#' @param site.info A data frame with site information
#' @param JL Scalar or vector number of individuals at each site, overrides JM
#' @param JM Total number of individuals to include in a MCSim simulation.
#' @param m Immigration rate paramter, from Hubbells neutral model. Overrides I.rate.m2.
#' @param I.rate.m2 Alternative to m, immigration rate in number of individuals / m2 / timestep. Default is 1.
#' @param area.m2 Area of each site in m2. Default is 1.
#' @param Ef.specificity Vector of specificity values for environmental filters at each site. If 0 (default), site habitat value is modeled as a single point along an environmental gradient. If > 0, a site's habitat is modeled as a normal curve around a point on an environmental gradient.
#' @param Ef Vector of habitat scores for each site.
#' @param guess.site.coords Binary. If TRUE, Uses PCoA to extract site coordinates if given a distance matrix or network map. Useful to make a map to display sites. Not necessary if igraph input is used because igraph has a function to plot a network map. Default is FALSE.
#' @param list.of.stuff A list that can be used to store other landscape attributes in the landscape object. Useful for storing igraph properties when igraph is used. 
#' 
#' @export
#' 
fn.make.landscape<-function(
  # -------------------------------
  # -------------------------------
  # -------------------------------
  # -- data frame inputs
  # -- need one or the other, or they have to match. Priority given to dist if they don't
  site.coords = c(1:10),
  dist.mat = data.frame(),
  
  # -- will fill in info if none given
  site.info = data.frame(),
  
  # -- default metacommunity parameters if none given
  JL = NULL,
  JM = NULL,
  m = NULL,
  I.rate.m2 = NULL,
  area.m2 = 1,
  Ef.specificity = 0, # 0 is point specificity along env. gradient
  Ef = 0,
  guess.site.coords = FALSE,
  list.of.stuff = NA,
  ...
){
  # -------------------------------
  # -------------------------------
  # -------------------------------
  
  # -- ESCAPE VAR
  get.the.f.out<-FALSE
  
  try({site.coords<-as.data.frame(site.coords)},
      silent = TRUE)
  
  # -------------------------------
  # -- Make geo landscape data.frame
  # -------------------------------
  dist.mat<-as.data.frame(as.matrix(dist.mat)) #make into a data.frame
  
  # --  check info
  if(nrow(dist.mat)>0){
    dat.geo.dist.out<-dist.mat
    if(nrow(site.coords)!=nrow(dist.mat) & guess.site.coords){
      dat.geo.out<-data.frame(cmdscale(dat.geo.dist.out))
      print('I assigned dat.geo for you')
    }else if(nrow(site.coords)!=nrow(dist.mat)){
      dat.geo.out<-data.frame(site.label=c(1:nrow(dat.geo.dist.out)))
      print('Rock on')
    }else if(nrow(site.coords)==nrow(dist.mat)){
      dat.geo.out<-site.coords
    }
  }else if(nrow(site.coords)>0){
    dat.geo.out<-site.coords
    dat.geo.dist.out<-data.frame(as.matrix(dist(site.coords)))
    print('gangsta')
  }else if(nrow(site.info)>0){
    get.the.f.out<-TRUE
    print('no geo info!')
  }else{
    get.the.f.out<-TRUE
    print('no landscape info!')
  }
  
  n.obs<-nrow(dat.geo.dist.out)
  
  if(!get.the.f.out){
    # -- turn area into a vector if it is not a vector, otherwise, it remains the same, get's fed up if it's the wrong length
    area.m2<-data.frame(
      dummy=c(1:n.obs),
      area.m2=area.m2)$area.m2
    
    # -------------------------------
    # -- calculate assemblage sizes at sites, JL 
    # -------------------------------
    if( is.null(JL) | !(length(JL)%in%c(1, n.obs)) ){
      JL.wts <- area.m2 / sum(area.m2)
      JL.wts <- JL.wts/sum(JL.wts)
      JL <- round(JL.wts * JM,0)
    }else{
      JL <- data.frame(
        dummy=c(1:n.obs),
        JL=JL)$JL
      JM <- sum(JL)
    }
    
    # -------------------------------
    # -- calculate immigration at sites, IL 
    # -------------------------------
    
    if(!is.null(m) & length(m)%in%c(1,n.obs)){
      m.site <- data.frame(
        dummy=c(1:n.obs),
        m=m)$m
      I.site <- m.site * (JL - 1) / (1 - m.site)
      I.site <- ifelse(is.infinite(I.site), 10*JL, I.site)
    }else if(!is.null(I.rate.m2)){
      I.site <- I.rate.m2 * area.m2
      I.site <- round(I.site,0)
      m.site <- I.site/(I.site + JL - 1)
    }else{
      print('aaaaggghhh, need immigration rate')
    }
    
    
    # -------------------------------
    # -- dat with info
    # -------------------------------
    dat.info.default<-data.frame(
      site.ID = c(1:nrow(dat.geo.out)),
      area.m2 = area.m2,
      JL = JL,
      Ef = Ef,
      Ef.specificity = Ef.specificity,
      IL = I.site,
      m = m.site
    )
    
    if(nrow(site.info)>0){
      dat.info.out<-site.info  
      # -- check to see if specific vars need to be filled in
      for(i.var in c('site.ID','area.m2','JL','Ef','Ef.specificity','IL','m')){
        if(!i.var%in%names(dat.info.out)) dat.info.out[,i.var] <- dat.info.default[,i.var]
      }
    }else{
      dat.info.out<-dat.info.default
    }
    
    return(
      list(site.info=dat.info.out,
           site.coords=dat.geo.out,
           dist.mat=dat.geo.dist.out,
           list.of.stuff=list.of.stuff)
    )
  }else{
    print('no landscape for you!')
  }
}


# ---------------------------------------------------------------------------------------
#' fn.recruit.Jt
#' 
#' @title Recruitment in a metacommunity with context
#' 
#' @description Internal function called by \link{fn.metaSIM} to calculate local 
#' recruitment pools for all sites in a metacommunity, and call \link{fn.lottery.recruit}
#'to create a new assemblages (J) at time (t) after accounting for metacommunity 
#'composition at time (t-1), dispersal dynamics, landscape topology, and environmental 
#'filtering.
#' 
#' @usage 
#' fn.recruit.Jt(landscape.site.coords = landscape$dat[, c("x", "y")], nu = nu,
#'               SWM.slope = SWM.slope, J.t.minus.1 = J.t.minus.1, 
#'               taxa.list = taxa.list, traits.Ef = dat.gamma.t0$trait.Ef, 
#'               trait.Ef.sd = trait.Ef.sd, traits.dispersal = dat.gamma.t0$trait.dispersal, 
#'               landscape.m = landscape$dat$m, landscape.Ef = landscape$dat$Ef, 
#'               landscape.Ef.specificity = landscape$dat$Ef.specificity, 
#'               landscape.JL = landscape$dat$JL)
#'               
#' @param landscape.site.coords A data.frame with xy-coordinates for each site.  Each row is a 
#' site in the metacommunity landscape.
#' @param nu Numeric, Hubbell's "speciation rate", but can be interpreted as the 
#' probability of the appearance of a novel species.  If set to 0, no novel taxa 
#' will appear during the simulation.
#' @param SWM.slope Slope of dispersal kernel, see \link{fn.metaSIM}.  This value 
#' is used for \code{w} in eqn. 6 from Gravel et al. (2006) to model dispersal 
#' limitation.
#' @param J.t.minus.1 Species composition of sites in a metacommunity at previous 
#' time step.
#' @param taxa.list A vector of character strings of species' names.
#' @param traits.Ef A vector of species' trait scores, numeric.
#' @param trait.Ef.sd A scalar value representing species' niche widths (numeric, 
#' currently one value is used for all species).
#' @param traits.dispersal A vector of species' dispersal trait scores.
#' @param landscape.m A vector of values for \code{m} for each site in the landscape.
#' @param landscape.Ef A vector of environmental filter \code{Ef} values for each site 
#' in the landscape.
#' @param landscape.Ef.specificity Environmental specificity at each site (numeric). 
#' @param landscape.JL A vector of assemblage sizes (positive integers) for all sites 
#' in the landscape.
#' 
#' @seealso \link{fn.lottery.recruit}, \link{fn.make.landscape}, \link{fn.metaSIM}
#' 
#' @export
#' 
fn.recruit.Jt <- function(
  # calls calculates R.t (expected pool) to calculate extant pool
  mat.geodist = as.matrix(dist(1:3)),
  nu = .001,
  SWM.slope = 0,
  JL = 10,
  J.t.minus.1 = matrix(JL/length(taxa.list), 
                       nrow=nrow(mat.geodist), 
                       ncol=length(taxa.list)),
  taxa.list = letters[1:5],
  traits.Ef = runif(length(taxa.list),0,1),
  trait.Ef.sd = 0.15,
  traits.dispersal = array(1, length(taxa.list)),
  m = 1,
  Ef = array(.5, nrow(mat.geodist)),
  Ef.specificity = 0){
  
  # -- Calculate R.t, estimates of locally available pools at each site based on 
  # local relative abundances at time t-1 (or t0)
  # regional relative abundances at time t-1 (or t0), which is calculated from neighbors weighted by the SWM
  # the balance of local and regional pools is determined by m (where m = 0 is only local)
  
  # deal with sites with 0s
  J.t.minus.1.RAs<-as.matrix(J.t.minus.1/rowSums(J.t.minus.1))
  J.t.minus.1.RAs[!is.finite(J.t.minus.1.RAs)]<-0
  J.t.minus.1.RAs<-as.data.frame(J.t.minus.1.RAs)
  
  # ----------------------------
  # -- calculate recruitment probabilities from I for each site
  # create SWM to weight all neighboring sites contributions to "Immigrant" pool
  # use eq 6 from Gravel et al. 2006, but set diags to 0
  #   mat.geodist<-as.matrix(dist(landscape.site.coords))
  
  max.val<-max(mat.geodist[mat.geodist!=Inf]) #max finite value
  mat.geodist.scaled<-mat.geodist/max.val
  SWM<-exp(-1*SWM.slope*mat.geodist.scaled^2)
  diag(SWM)<-0
  SWM<-SWM/rowSums(SWM)
  
  
  # -- calculate I for each site based on composition of neighboring sites
  I.RAs<-SWM%*%as.matrix(J.t.minus.1.RAs) #distance decay without species bias
  
  # -- Alter RAs for I based on dispersal traits -- create bias toward species with better dispersal
  I.RAs.dispersal.biased<-I.RAs*(t(traits.dispersal%o%rep(1,nrow(I.RAs))))
  I.RAs.dispersal.biased<-I.RAs.dispersal.biased/rowSums(I.RAs.dispersal.biased)
  I.RAs.dispersal.biased[is.nan(I.RAs.dispersal.biased)]<-0
  # -- calculate recruitment pool by combining local RAs and regional RAs, weighted by m
  R.t<-m*I.RAs.dispersal.biased+(1-m)*J.t.minus.1.RAs
  
  # -- all species that have 0 regional RA at time t-1 are assigned a non-zero probability of recruitment,
  # which is nu / (count of unobserved species in the list).  A "speciation event" in this simulation is 
  # recruitment of a previously unobserved species. 
  unobserved.spp.count <- sum(colSums(R.t)==0)
  speciation.recruitment.prob <- ifelse(
    unobserved.spp.count == 0,
    0,
    nu / unobserved.spp.count)
  
  # -- alter all the probabilities accordingly, all rows must add to 1
  mat.nu<-R.t
  mat.nu[,colSums(mat.nu)>0]<-mat.nu[,colSums(mat.nu)>0]*(1-nu)
  mat.nu[,colSums(mat.nu)==0]<-speciation.recruitment.prob
  R.t.probs<-mat.nu
  
  # -- Local environmental filtering
  d.temp<-expand.grid(Ef=Ef,stringsAsFactors=FALSE,
                      trait.Ef=traits.Ef)
  d.temp$Ef.specificity<-Ef.specificity
  
  lambda.Ef.siteBYspp<-matrix(
    data=mapply(FUN=fn.lambda,
                trait.optimum=d.temp$trait.Ef,
                Ef=d.temp$Ef,
                Ef.specificity=d.temp$Ef.specificity,
                MoreArgs=list(niche.breadth=trait.Ef.sd)
    ),
    nrow=length(Ef), #number sites
    ncol=length(traits.Ef), #number spp
    byrow=FALSE
  )
  lambda.Ef.siteBYspp<-lambda.Ef.siteBYspp/rowSums(lambda.Ef.siteBYspp) #rescale so row sums are 1
  
  # reweight R.t.probs based on local env. filters
  R.t.envfiltered<-lambda.Ef.siteBYspp*R.t.probs
  R.t.envfiltered<-R.t.envfiltered/rowSums(R.t.envfiltered)
  
  # make lists of recruitment weights for each site for lottery call below
  R.t.list<-as.list(data.frame(t(R.t.envfiltered)))
  
  # -- Recruit extant community for time t from R.t
  J.t1<-data.frame(
    row.names=NULL,
    t(mapply(
      FUN=fn.lottery.recruit,
      vect.recruitment.weights=R.t.list,
      scalar.JL=as.list(JL),
      MoreArgs=list(vect.taxa.list=taxa.list)
    )))
  return(J.t1)
}

# --------------------------------------------------------------------------------------------
#' fn.set.regional.species.pool
#' 
#' @title Initialize species pool for fn.metaSIM
#' 
#' @description Initialize species pool for fn.metaSIM
#' 
#' @usage 
#' fn.set.regional.species.pool(alpha.fisher = 1)
#' fn.set.regional.species.pool(gamma.abund = c(.9, .4, .2, .3))

#' @param alpha.fisher Fisher's alpha used to initiate the simulation.
#' @param n.timestep Number of generations (time steps) in the simulation, can affect number of novel species. Default is 0.
#' @param nu Scalar value representing the probability that a novel species 
#' is recruited.  Hubbell's \dQuote{speciation rate}.
#' @param speciation.limit A limit to the number of novel species that can occur 
#' in the simulation, default is 0.  
#' @param JM Total number of individuals in the metacommunity.
#' @param gamma.abund Vector of relative abundances or counts used to set regional species pool RAs
#' @param trait.dispersal Vector of trait scores for dispersal, should be in range 0 to 1, or NULL 
#' @param trait.dispersal.median Median value for species' dispersal trait values.
#' @param trait.dispersal.range Range in variation allowed for species' trait 
#' dispersal values.
#' @param Ef.min Min limit for environmental gradient
#' @param Ef.max Max limit for environmental gradient
#' @param trait.Ef Vector or value of niche positions for species in regional pool
#' @param trait.Ef.sd Vector or value of niche breadths
#' @param taxa.list.prefix A character string to use as a prefix in species' names.
#' @param taxa.list.prefix Text appended as prefix to species names to prevent errors for column names when site by species matrices are created. Default is "spp".
#' 
#' @export
#' 
fn.set.regional.species.pool <- function ( alpha.fisher = 1,
                                           n.timestep = 0, 
                                           nu = 0.001, 
                                           speciation.limit = 0, 
                                           JM = 1000, 
                                           
                                           gamma.abund = NULL,
                                           
                                           trait.dispersal = NULL, # a vector of dispersal traits
                                           trait.dispersal.median = 1,# or set a median and range of dispersals
                                           trait.dispersal.range = 0, 
                                           
                                           Ef.min = 0, 
                                           Ef.max = 1,
                                           
                                           trait.Ef = NULL,
                                           trait.Ef.sd = NULL,
                                           
                                           taxa.list = NULL,
                                           taxa.list.prefix = "spp"){
  
  # -- set regional RAs if given, otherwise sample dist based on alpha.fisher
  if(is.array(gamma.abund)|is.vector(gamma.abund)){
    alpha.fisher <- NA
  } else {
    # -- default alpha.fisher to 1 when no gamma.abund is provided
    if(is.na(alpha.fisher)) {alpha.fisher <- 1}
    gamma.abund <- fisher.ecosystem(N = JM, nmax = JM, alpha = alpha.fisher)
  }
  
  # -- add columsn for potential new species
  if (is.na(speciation.limit)) {
    n.new.spp <- length(gamma.abund)
  } else {
    n.new.spp <- nu * n.timestep * JM
    if (n.new.spp > speciation.limit) 
      n.new.spp <- speciation.limit
  }
  
  # -- Rel Abundances
  regional.RA <- c(gamma.abund/sum(gamma.abund), rep(0, n.new.spp))
  total.richness <- length(regional.RA)
  
  # -- dispersal traits -- 1 is max dispersal, 0 is min
  while(length(trait.dispersal)<total.richness){
    if(!(is.array(trait.dispersal)|is.vector(trait.dispersal))){
      trait.dispersal <- runif(total.richness, 
                               trait.dispersal.median - trait.dispersal.range, 
                               trait.dispersal.median + trait.dispersal.range)
      trait.dispersal[trait.dispersal > 1] <- 1
      trait.dispersal[trait.dispersal < 0] <- 0
    } else {
      trait.dispersal<-c(trait.dispersal,
                         runif(total.richness-length(trait.dispersal), 
                               trait.dispersal.median - trait.dispersal.range, 
                               trait.dispersal.median + trait.dispersal.range))
      trait.dispersal[trait.dispersal > 1] <- 1
      trait.dispersal[trait.dispersal < 0] <- 0
    }
  }
  
  # -- assign habitat preference traits
  while(length(trait.Ef)<total.richness){
    if(!(is.array(trait.Ef)|is.vector(trait.Ef))){
      trait.Ef <- runif(total.richness, Ef.min, Ef.max)
    } else {
      trait.Ef <- c(trait.Ef,
                  runif(total.richness-length(trait.Ef), 
                        Ef.min, Ef.max))
    }
  }
  
  # -- assign niche breadths
  while(length(trait.Ef.sd)<total.richness){
    if(!(is.array(trait.Ef.sd)|is.vector(trait.Ef.sd))){ #what to do if trait.Ef.sd array is length 0, make NCM
      
      # -- default is neutral community model
      trait.range <- max(trait.Ef)-min(trait.Ef)
      if(!trait.range>0) trait.range <- 1 #make non-zero trait range
      trait.Ef.sd <- rep(1000*trait.range, total.richness) #make neutral community
      
    } else {
      trait.Ef.sd<-c(trait.Ef.sd,
                     sample(x = c(trait.Ef.sd, trait.Ef.sd),
                            size = total.richness-length(trait.Ef.sd), 
                            replace = TRUE))
    }
  }
  
  # -- taxa.list
  while(length(taxa.list)<total.richness){
    if(!(is.array(taxa.list)|is.vector(taxa.list))){
      taxa.list <- paste(taxa.list.prefix, c(1:total.richness), 
                         sep = "")
    } else {
      taxa.list <- c(taxa.list,
                     paste(taxa.list.prefix, 
                           c((length(taxa.list)+1):total.richness), 
                           sep = ""))
    }
  }
    
  if(length(taxa.list)>total.richness){taxa.list<-taxa.list[1:total.richness]}
  
  # -- return dataframe
  return(data.frame(row.names = taxa.list, 
                    taxa.list = taxa.list, 
                    trait.dispersal = trait.dispersal[1:length(taxa.list)], 
                    trait.Ef = trait.Ef[1:length(taxa.list)], 
                    trait.Ef.sd = trait.Ef.sd[1:length(taxa.list)],
                    regional.RA = c(regional.RA), 
                    stringsAsFactors = FALSE))
}
