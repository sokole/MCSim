#' fn.true.dist
#' 
#' @title Custom function to calculate true distances (pairwise beta-diversities)
#' 
#' @description Custom function to calculate true distances based on Jost 2006, 
#' Jost 2007, and Chao et al. 2012.
#' 
#' @usage 
#' fn.true.dist(d = dat.comm, q.order = 1)
#' 
#' @param d Site by species data matrix or data frame.
#' @param q.order Order \code{q}
#' 
#' @references Chao, A., C.-H. Chiu, and T. C. Hsieh. 2012. Proposing a resolution 
#' to debates on diversity partitioning. Ecology 93:2037--2051.
#' 
#' @export
#' 
fn.true.dist <- function(
  # used in variation partitioning when 
  d=dat.comm, #site X species community matrix
  q.order=1 ){ #scale for order q
  # -- note calls fn.CqN.overlap, calcluates overlap as described by Jost et al. 2011, Ch6 in Biological Diversity 
  require(plyr)
  # -- expanded matrix
  d2<-expand.grid(site1=c(1:nrow(d)),site2=c(1:nrow(d)))
  dist.CqN<-maply(
    .data=d2,
    .fun=function(site1,site2,d,q.order) fn.CqN.overlap(d[c(site1,site2),],q.order),
    d=d,
    q.order=q.order)
  dimnames(dist.CqN)<-NULL
  return(1-dist.CqN)
}
