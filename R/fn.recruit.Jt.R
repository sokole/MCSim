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
#' fn.recruit.Jt(landscape.xy = landscape$dat[, c("x", "y")], nu = nu,
#'               SWM.slope = SWM.slope, J.t.minus.1 = J.t.minus.1, 
#'               taxa.list = taxa.list, traits.Ef = dat.gamma.t0$trait.Ef, 
#'               trait.Ef.sd = trait.Ef.sd, traits.dispersal = dat.gamma.t0$trait.dispersal, 
#'               landscape.m = landscape$dat$m, landscape.Ef = landscape$dat$Ef, 
#'               landscape.Ef.specificity = landscape$dat$Ef.specificity, 
#'               landscape.JL = landscape$dat$JL)
#'               
#' @param landscape.xy A data.frame with xy-coordinates for each site.  Each row is a 
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
  landscape.xy=landscape$dat[,c("x","y")],
  nu=nu,
  SWM.slope=SWM.slope,
  J.t.minus.1=J.t.minus.1,
  taxa.list=taxa.list,
  traits.Ef=dat.gamma.t0$trait.Ef,
  trait.Ef.sd=trait.Ef.sd,
  traits.dispersal=dat.gamma.t0$trait.dispersal,
  landscape.m=landscape$dat$m,
  landscape.Ef=landscape$dat$Ef,
  landscape.Ef.specificity=landscape$dat$Ef.specificity,
  landscape.JL=landscape$dat$JL){
  
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
  mat.geodist<-as.matrix(dist(landscape.xy))
  mat.geodist.scaled<-mat.geodist/max(mat.geodist)
  SWM<-exp(-1*SWM.slope*mat.geodist.scaled^2)
  diag(SWM)<-0
  SWM<-SWM/rowSums(SWM)
  
  # -- calculate I for each site based on composition of neighboring sites
  I.RAs<-SWM%*%as.matrix(J.t.minus.1.RAs) #distance decay without species bias
  
  # -- Alter RAs for I based on dispersal traits -- create bias toward species with better dispersal
  I.RAs.dispersal.biased<-I.RAs*(t(traits.dispersal%o%rep(1,nrow(I.RAs))))
  I.RAs.dispersal.biased<-I.RAs.dispersal.biased/rowSums(I.RAs.dispersal.biased)
  
  # -- calculate recruitment pool by combining local RAs and regional RAs, weighted by m
  R.t<-landscape.m*I.RAs.dispersal.biased+(1-landscape.m)*J.t.minus.1.RAs
  
  # -- all species that have 0 regional RA at time t-1 are assigned a non-zero probability of recruitment,
  # which is nu / (coung of unobserved species in the list).  A "speciation event" in this simulation is 
  # recruitment of a previously unobserved species.  
  speciation.recruitment.prob<-nu/(sum(colSums(R.t)==0))
  
  # -- alter all the probabilities accordingly, all rows must add to 1
  mat.nu<-R.t
  mat.nu[,colSums(mat.nu)>0]<-mat.nu[,colSums(mat.nu)>0]*(1-nu)
  mat.nu[,colSums(mat.nu)==0]<-speciation.recruitment.prob
  R.t.probs<-mat.nu
  
  # -- Local environmental filtering
  d.temp<-expand.grid(Ef=landscape.Ef,stringsAsFactors=FALSE,
                      trait.Ef=traits.Ef)
  d.temp$Ef.specificity<-landscape.Ef.specificity
  
  lambda.Ef.siteBYspp<-matrix(
    data=mapply(FUN=fn.lambda,
                trait.optimum=d.temp$trait.Ef,
                Ef=d.temp$Ef,
                Ef.specificity=d.temp$Ef.specificity,
                MoreArgs=list(niche.breadth=trait.Ef.sd)
    ),
    nrow=length(landscape.Ef), #number sites
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
      scalar.JL=as.list(landscape.JL),
      MoreArgs=list(vect.taxa.list=taxa.list)
    )))
  return(J.t1)
}
