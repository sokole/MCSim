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