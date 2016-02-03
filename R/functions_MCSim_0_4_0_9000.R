# ---------------------------------------------------------------------------------------
# -- v0.4.0.9000 dev functions
# ---------------------------------------------------------------------------------------
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
#' fn.metaSIM3(landscape)
#' 
#' @param landscape landscape list from \link{fn.make.landscape3}
#' @param scenario.ID see \link{fn.metaSIM}
#' @param alpha.fisher see \link{fn.metaSIM}
#' @param nu see \link{fn.metaSIM}
#' @param speciation.limit see \link{fn.metaSIM}
#' @param n.timestep see \link{fn.metaSIM}
#' @param SWM.slope see \link{fn.metaSIM}
#' @param trait.dispersal.median see \link{fn.metaSIM}
#' @param trait.dispersal.range see \link{fn.metaSIM}
#' @param trait.Ef.sd see \link{fn.metaSIM}
#' @param save.sim see \link{fn.metaSIM}
#' @param output.dir.path see \link{fn.metaSIM}
#' @param \dots other parameters to be passed to internal functions
#' 
#' @seealso \link{fn.metaSIM}
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
fn.metaSIM<-function(
  landscape = NA,
  scenario.ID = NA, 
  alpha.fisher = 2,
  nu = 1e-04,
  speciation.limit = NA,
  JM.limit = 1e5, # may need to put a limit on how many individuals can be in a simulation
  n.timestep = 10,
  SWM.slope = 0,
  trait.dispersal.median = 1 ,
  trait.dispersal.range = 0,
  trait.Ef.sd = 0.3,
  save.sim = TRUE, 
  output.dir.path = "SIM_OUTPUT", 
  ...
){
  if(class(landscape)!='list'){
    print('This simulation needs a landscape!')
  }else{
    try({  
      # --------------------------------------------------------
      # -- calculate JM from landscape
      # --------------------------------------------------------
      attach(landscape$site.info)
      JM <- sum(JL)
      n.sites <- length(JL)
      
      #' -- create regional pool
      dat.gamma.t0 <- fn.set.regional.species.pool(n.timestep = n.timestep, 
                                                   nu = nu, 
                                                   speciation.limit = speciation.limit, 
                                                   JM = ifelse(JM>JM.limit,JM.limit,JM), 
                                                   alpha.fisher = alpha.fisher, 
                                                   trait.dispersal.median = trait.dispersal.median, 
                                                   trait.dispersal.range = trait.dispersal.range)
      taxa.list <- as.character(dat.gamma.t0$taxa.list)
      d.temp <- expand.grid(Ef = Ef, stringsAsFactors = FALSE, 
                            trait.Ef = dat.gamma.t0$trait.Ef)
      d.temp$Ef.specificity <- Ef.specificity
      lambda.Ef.siteBYspp <- matrix(data = mapply(FUN = fn.lambda, 
                                                  trait.optimum = d.temp$trait.Ef, Ef = d.temp$Ef, Ef.specificity = d.temp$Ef.specificity, 
                                                  MoreArgs = list(niche.breadth = trait.Ef.sd)), nrow = n.sites, 
                                    ncol = length(taxa.list), byrow = FALSE)
      lambda.Ef.siteBYspp <- lambda.Ef.siteBYspp/rowSums(lambda.Ef.siteBYspp)
      
      R.probs.t0 <- lambda.Ef.siteBYspp * (rep(1, n.sites) %o% 
                                             dat.gamma.t0$regional.RA)
      R.probs.t0 <- R.probs.t0/rowSums(R.probs.t0)
      R.probs.list <- as.list(data.frame(t(R.probs.t0)))
      
      J.t0 <- data.frame(row.names = NULL, t(mapply(FUN = fn.lottery.recruit, 
                                                    vect.recruitment.weights = R.probs.list, scalar.JL = as.list(JL), 
                                                    MoreArgs = list(vect.taxa.list = taxa.list))))
      J <- list()
      J[[1]] <- J.t0
      J.t.minus.1 <- J.t0
      
      # ----------------------------------------------------------------------
      # loop for generational turnover in metacommunity
      # ----------------------------------------------------------------------
      for (t.index in 2:n.timestep) {
        J.t <- fn.recruit.Jt(mat.geodist=landscape$dist.mat,
                                   nu=nu,
                                   SWM.slope=SWM.slope,
                                   J.t.minus.1=J.t.minus.1,
                                   taxa.list=taxa.list,
                                   traits.Ef=dat.gamma.t0$trait.Ef,
                                   trait.Ef.sd=trait.Ef.sd,
                                   traits.dispersal=dat.gamma.t0$trait.dispersal,
                                   m=m,
                                   Ef=Ef,
                                   Ef.specificity=Ef.specificity,
                                   JL=JL)
        J[[t.index]] <- J.t
        J.t.minus.1 <- J.t
        print(paste("Timestep:", t.index))
      }
      
      # ----------------------------------------------------------------------
      sim.result.name <- paste("SIM_", scenario.ID, "_", format(Sys.time(), 
                                                                "%Y%m%d_%H%M%S"), "_", trunc(runif(1, 1e+05, 999999)), 
                               sep = "")
      
      sim.result <- list(scenario.ID = scenario.ID, 
                         sim.result.name = sim.result.name, 
                         alpha.fisher = alpha.fisher, 
                         nu.sim = nu, 
                         trait.Ef.sd = trait.Ef.sd, 
                         trait.dispersal = dat.gamma.t0$trait.dispersal, 
                         trait.Ef = dat.gamma.t0$trait.Ef, 
                         landscape = landscape, 
                         dat.gamma.t0 = dat.gamma.t0, 
                         SWM.slope = SWM.slope, 
                         J = J, 
                         n.timestep = n.timestep, 
                         taxa.list = taxa.list)
      
      #'   fn.sim.metadata.archive4(sim.result = sim.result, 
      #'                            save.sim = save.sim, var.dir = output.dir.path, 
      #'                            keep.timesteps = keep.timesteps,
      #'                            q.order=NA,
      #'                            ...)
      sim.result.filename<-paste(output.dir.path,"/",sim.result.name,".rda",sep="")
      sim.result.metadata <- data.frame(row.names = sim.result.name, 
                                        scenario.ID = scenario.ID, 
                                        sim.ID = sim.result.name, 
                                        sim.result.filename = sim.result.filename, 
                                        n.sites = n.sites, 
                                        n.timestep = n.timestep, 
                                        alpha.fisher = alpha.fisher, 
                                        nu.sim = nu, 
                                        JM = JM, 
                                        JL.mean = mean(JL), 
                                        JL.sd = sd(JL), 
                                        m.mean = mean(m), 
                                        m.sd = sd(m), 
                                        IL.mean = mean(IL), 
                                        IL.sd = sd(IL), 
                                        SWM.slope = SWM.slope, 
                                        Ef.mean = mean(Ef), 
                                        Ef.sd = sd(Ef), 
                                        Ef.specificity.mean = mean(Ef.specificity), 
                                        Ef.specificity.sd = sd(Ef.specificity), 
                                        Tr.disp.mean = mean(dat.gamma.t0$trait.dispersal), 
                                        Tr.disp.sd = sd(dat.gamma.t0$trait.dispersal), 
                                        Niche.breadth = sim.result$trait.Ef.sd, 
                                        stringsAsFactors = FALSE)
      
      #' -- check for director
      if(!output.dir.path%in%list.files())  dir.create(output.dir.path)
      
      #' -- save sim
      if(save.sim) save(sim.result,file=sim.result.filename)
      
      #' -- check to see if data for other reps from this scenario exist
      filname.sim.metadata <- paste(output.dir.path, "/sim.metadata_", 
                                    sim.result$scenario.ID, ".csv", sep = "")
      
      #' -- create new file if none exists, write results to file, otherwise append to existing file
      if (file.exists(filname.sim.metadata)) {
        dat.sim.metadata <- read.csv(filname.sim.metadata, row.names = 1, 
                                     header = TRUE)
        if (!sim.result.name %in% row.names(dat.sim.metadata)) {
          dat.sim.metadata <- rbind(dat.sim.metadata, sim.result.metadata)
          write.csv(dat.sim.metadata, filname.sim.metadata)
        }
      }else {
        write.csv(sim.result.metadata, filname.sim.metadata)
      }
      
      try(detach(landscape$site.info), silent = TRUE)
      try(detach(landscape), silent = TRUE)
      
      #' -- return results
      return(sim.result)
    })
  }
} 
#end function

# ---------------------------------------------------------------------------------------
#' fn.make.landscape
#' 
#' @title make a simulation landscape
#' 
fn.make.landscape<-function(
  # -------------------------------
  # -------------------------------
  # -------------------------------
  # -- data frame inputs
  # -- need one or the other, or they have to match. Priority given to dist if they don't
  site.coords = data.frame(),
  dist.mat = data.frame(),
  
  # -- will fill in info if none given
  site.info = data.frame(),
  
  # -- default metacommunity parameters if none given
  JM = 1000,
  I.rate.m2 = 1,
  area.m2 = 1,
  Ef.specificity = 0, # 0 is point specificity along env. gradient
  Ef = .5,
  guess.site.coords = FALSE,
  list.of.stuff = NA
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
    # -- JL influenced by management
    # -------------------------------
    JL.wts <- area.m2 / sum(area.m2)
    JL.wts <- JL.wts/sum(JL.wts)
    JL <- round(JL.wts * JM,0)
    
    # -------------------------------
    # -- calculate immigration at sites, IL 
    # -------------------------------
    I.site <- I.rate.m2 * area.m2
    I.site <- round(I.site,0)
    m.site <- I.site/(I.site + JL - 1)
    
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
  SWM<-SWM/rowSums(SWM)
  diag(SWM)<-0
  
  
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
