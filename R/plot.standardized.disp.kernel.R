#' @title plot standardized dispersal kernel
#' 
#' @description plot standardized dispersal kernel. Requires either (1) a sim.result, (2) a landsca;e and dispersal kernel slope(s) (w), or (3) a distance.matrix and dispersal kernel slope(s) (w).
#' 
#' @usage 
#' plot.standardized.disp.kernel(sim.result)
#' 
#' @param w dispersal kernel slope, default is NULL
#' @param sim.result simulation output, default is NULL
#' @param landscape a landscape object, default is NULL
#' @param distance.matrix a distance matrix, default is NULL
#' 
#' 
#' 
#' @export plot.standardized.disp.kernel
#' 
plot.standardized.disp.kernel <- function(
  sim.result = NULL,
  w = NULL,
  landscape = NULL,
  distance.matrix = NULL,
  ...){
  
  # requireNamespace("magrittr")
  
  w.ordered <- c(w, sim.result$W.r) %>% unique() %>% sort(decreasing = TRUE) %>% na.omit()
  stopifnot({length(w.ordered) > 0})
  
  w.max <- max(w.ordered)
  
  curve(exp(-1*w.max*x^2),0,1,
        ylim=c(0,1.05),
        ylab='W(r)',
        xlab='r')
  
  if(length(w.ordered) > 2){
    for(i in 2:length(w.ordered)){
      
      curve(exp(-1*w.ordered[i]*x^2),0,1,
            ylim = c(0,1.15),
            lty = i,
            lwd = 1.25,
            add = TRUE)
    }
    
    legend(x = 'topright',
           legend = w.ordered,
           lty = c(1:length(w.ordered)),
           lwd = 1.25,
           bg = ('white'),
           title = 'w = ',
           cex = .85)
  }

  if(!is.null(sim.result)){
    dist.vect <- sim.result$landscape$dist.mat %>% unlist()
    rug(dist.vect/max(dist.vect),
        side=1)
  }else if(!is.null(landscape)){
    dist.vect <- landscape$dist.mat %>% unlist()
    rug(dist.vect/max(dist.vect),
        side=1)
  }else if(!is.null(distance.matrix)){
    dist.vect <- distance.matrix %>% unlist()
    rug(dist.vect/max(dist.vect),
        side=1)
  }
    
}#END FUNCTION

# median(dist.vect/max(dist.vect))