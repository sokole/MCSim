#' fn.diversity.partition
#' 
#' @title Diversity partitioning
#' 
#' @description Function called by \link{fn.metaSIM} to calculate alpha, 
#' multiplicative beta, and gamma diversity for a site by species matrix 
#' based on Jost (2007) and Chao et al. (2012).  This function relies on 
#' the \link[vegetarian]{d} function in the \pkg{vegetarian} package.
#' 
#' @usage fn.diversity.partition(dat.comm, q.order = 1)
#' 
#' @param dat.comm
#' @param q.order
#' 
#' @references 
#' Chao, A., C. H. Chiu, and T. C. Hsieh. 2012. Proposing a resolution to debates 
#' on diversity partitioning. Ecology 93:2037--2051.
#' 
#' Jost, L. 2007. Partitioning diversity into independent alpha and beta components. 
#' Ecology 88:2427--2439.
#' 
#' see the \href{http://cran.r-project.org/web/packages/vegetarian/index.html}{vegetarian} 
#' package
#' 
#' @seealso \link[vegetarian]{d} in the \pkg{vegetarian} package
#' 
#' @export
fn.diversity.partition <- function(     # Diversity partitioning
  dat.comm,
  q.order=1){
  
  require(vegetarian)
  
  div.results<-data.frame(
    d.alpha=d(dat.comm,lev="alpha",
              wts=rowSums(dat.comm),
              q=q.order),
    d.beta=d(dat.comm,lev="beta",
             wts=rowSums(dat.comm),
             q=q.order),
    d.gamma=d(dat.comm,lev="gamma",
              wts=rowSums(dat.comm),
              q=q.order)
  )
  
  names(div.results)<-paste(names(div.results),".q",q.order,sep="")
  
  return(div.results)
}
