#' fn.set.regional.species.pool
#' 
#' @title Initialize species pool for fn.metaSIM
#' 
#' @description Initialize species pool for fn.metaSIM
#' 
#' @usage 
#' fn.set.regional.species.pool(n.timestep = n.timestep, nu = nu,
#' speciation.limit = NA, JM, alpha.fisher, trait.dispersal.median,
#' trait.dispersal.range, taxa.list.prefix = "spp")

#' @param n.timestep Number of generations (time steps) in the simulation.
#' @param nu Scalar value representing the probability that a novel species 
#' is recruited.  Hubbell's \dQuote{speciation rate}.
#' @param speciation.limit A limit to the number of novel species that can occur 
#' in the simulation.  
#' @param JM Total number of individuals in the metacommunity.
#' @param alpha.fisher Fisher's alpha used to initiate the simulation.
#' @param trait.dispersal.median Median value for species' dispersal trait values.
#' @param trait.dispersal.range Range in variation allowed for species' trait 
#' dispersal values.
#' @param taxa.list.prefix A character string to use as a prefix in species' names.
#' 
#' @export
#' 
fn.set.regional.species.pool <- function(
  # Function to initialize the regional pool.  Set richness and traits.
  n.timestep=n.timestep,
  nu=nu,
  speciation.limit=NA,
  JM,
  alpha.fisher,
  trait.dispersal.median,
  trait.dispersal.range,
  taxa.list.prefix="spp"){
  
  gamma.abund<-fisher.ecosystem(N=JM,nmax=JM,alpha=alpha.fisher) #use Fisher's logseries to create a regional rank abundance curve
  
  # -- estimate number of new species via speciation
  if(is.na(speciation.limit)){
    n.new.spp<-length(gamma.abund)
  }else{
    n.new.spp<-nu*n.timestep*JM
    if(n.new.spp>speciation.limit) n.new.spp<-speciation.limit
  }
  regional.RA<-c(gamma.abund/sum(gamma.abund),rep(0,n.new.spp))
  total.richness<-length(regional.RA)
  
  # ---------
  # -- TRAITS
  # -- dispersal, Regional bias based on dispersal ability
  # dispersal vals range [0,1], where 0 is no dispersal and 1 
  # is best dispersal (dispersal = m)
  # default is m = 1 for all species
  trait.dispersal<-runif(total.richness,
                         trait.dispersal.median-trait.dispersal.range,
                         trait.dispersal.median+trait.dispersal.range)
  trait.dispersal[trait.dispersal>1]<-1         #set range to [0,1]
  trait.dispersal[trait.dispersal<0]<-0
  
  # -- local bias based on environmental filtering (Ef)
  # current implimentation is...
  # optimal niche value in [0,1] (using gausian fucntion for filtering), see Gravel et al. 2006, eq 3
  trait.Ef<-runif(total.richness,0,1)
  
  # -- make taxa list
  taxa.list<-paste(taxa.list.prefix,c(1:total.richness),sep="")
  
  return(
    data.frame(
      row.names=taxa.list,
      taxa.list=taxa.list,
      trait.dispersal=trait.dispersal,
      trait.Ef=trait.Ef,
      regional.RA=c(regional.RA),
      stringsAsFactors=FALSE
    )
  )
}
