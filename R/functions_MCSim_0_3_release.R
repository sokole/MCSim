# --------------------------------------------------------------------------------------------
#' fn.lambda
#' 
#' @title Calculate recruitment bias based on environmental filtering
#' 
#' @description Used by \link{fn.metaSIM}.  Uses niche breadth and an 
#' environmental gradient to determine how to bias species selection 
#' weights based on environmental filtering dynamics.  See Gravel et al. 2006.
#' 
#' @usage fn.lambda(trait.optimum, niche.breadth, Ef, Ef.specificity)
#' 
#' @param trait.optimum Optimum value for a given species along an environmental 
#' gradient.  Gradients in \link{fn.metaSIM} are restricted to the range [0,1], so 
#' this value must be in the range [0,1].
#' @param niche.breadth Niche breadth around a species' trait optimum.  The value of sigma 
#' in Fig 1 in Gravel et al. (2006).
#' @param Ef Value of the environmental filter at the site for which lambda is being 
#' calculated.
#' @param Ef.specificity The selection specificity of the environmental filter.
#' 
#' @references 
#' Gravel, D., C. D. Canham, M. Beaudet, and C. Messier. 2006. Reconciling niche and 
#' neutrality: the continuum hypothesis. Ecology Letters 9:399--409.
#' 
#' @export
#' 
fn.lambda <- function(
  trait.optimum = .5,
  niche.breadth = 1,
  Ef = .5,
  Ef.specificity = 0){
  if(Ef.specificity!=0){
    integrate(f=function(e,t,niche)exp(-1 * ((e -  t)^2) / (2 * niche^2) ),
              lower=Ef-(Ef.specificity/2),
              upper=Ef+(Ef.specificity/2),
              t=trait.optimum,
              niche=niche.breadth)$value
  }else{
    exp(-1 * ((Ef -  trait.optimum)^2) / (2 * niche.breadth^2) )
  }
}

# --------------------------------------------------------------------------------------------
#' fn.lottery.recruit
#' 
#' @title Lottery recruitment
#' 
#' @description Lottery recruitment function for a single assemblage used 
#' by \link{fn.metaSIM}.
#' 
#' @usage fn.lottery.recruit(vect.recruitment.weights, vect.taxa.list, scalar.JL)
#' 
#' @param vect.recruitment.weights
#' @param vect.taxa.list A vector of characters (or character strings) that is a 
#' list of species names.
#' @param scalar.JL Size of the assemblage (i.e., number of individuals to recruit).
#' 
#' @export
#' 
fn.lottery.recruit <- function(
  vect.recruitment.weights = rep(1,length(vect.taxa.list)),
  vect.taxa.list = "spp1",   # -- recruitment pool with weights
  scalar.JL = 100L
){
  # -- function to recruit from a source pool, given site carrying capacities and a recruitment pool
  dat.abund<-table(as.character(
    sample(x=vect.taxa.list,                            
           size=scalar.JL,           
           prob=vect.recruitment.weights,
           replace=TRUE)))[vect.taxa.list]
  names(dat.abund)<-vect.taxa.list
  dat.abund[is.na(dat.abund)]<-0
  return(dat.abund) 
}

# --------------------------------------------------------------------------------------------
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
  n.timestep = 0,
  nu=.001,
  speciation.limit=NA,
  JM = 1000,
  alpha.fisher = 2,
  trait.dispersal.median = 1,
  trait.dispersal.range = 0,
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
