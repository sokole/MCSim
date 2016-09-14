# --------------------------------------------------------------------------------------------
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
#' @details Used by \link{fn.metaSIM}.  
#' 
#' @references 
#' Gravel, D., C. D. Canham, M. Beaudet, and C. Messier. 2006. Reconciling niche and 
#' neutrality: the continuum hypothesis. Ecology Letters 9:399--409.
#' 
#' @export
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
#' 
#' @title Lottery recruitment
#' 
#' @description Lottery recruitment function for a single assemblage used 
#' by \link{fn.metaSIM}.
#' 
#' @usage fn.lottery.recruit(vect.recruitment.weights, vect.taxa.list, scalar.JL)
#' 
#' @param vect.recruitment.weights A vector of recruitment weights.
#' @param vect.taxa.list A vector of characters (or character strings) that is a 
#' list of species names.
#' @param scalar.JL Size of the assemblage (i.e., number of individuals to recruit).
#' 
#' @details Used by \link{fn.metaSIM}.
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