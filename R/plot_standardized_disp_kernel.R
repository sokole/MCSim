#' @title plot standardized dispersal kernel
#' 
#' @description plot standardized dispersal kernel
#' 
#' @param w dispersal kernel slope, default is NULL
#' @param sim_result simulation output, default is NULL
#' @param landscape a landscape object, default is NULL
#' @param distance_matrix a distance matrix, default is NULL
#' 
#' 
#' 
#' @export
plot_standardized_disp_kernel <- function(
  sim_result = NULL,
  w = NULL,
  landscape = NULL,
  distance_matrix = NULL,
  ...){
  
  # requireNamespace("magrittr")
  
  w_ordered <- c(w, sim_result$W.r) %>% unique() %>% sort(decreasing = TRUE) %>% na.omit()
  stopifnot({length(w_ordered) > 0})
  
  w_max <- max(w_ordered)
  
  curve(exp(-1*w_max*x^2),0,1,
        ylim=c(0,1.05),
        ylab='W(r)',
        xlab='r')
  
  if(length(w_ordered) > 2){
    for(i in 2:length(w_ordered)){
      
      curve(exp(-1*w_ordered[i]*x^2),0,1,
            ylim = c(0,1.15),
            lty = i,
            lwd = 1.25,
            add = TRUE)
    }
    
    legend(x = 'topright',
           legend = w_ordered,
           lty = c(1:length(w_ordered)),
           lwd = 1.25,
           bg = ('white'),
           title = 'w = ',
           cex = .85)
  }

  if(!is.null(sim_result)){
    dist.vect <- sim_result$landscape$dist.mat %>% unlist()
    rug(dist.vect/max(dist.vect),
        side=1)
  }else if(!is.null(landscape)){
    dist.vect <- landscape$dist.mat %>% unlist()
    rug(dist.vect/max(dist.vect),
        side=1)
  }else if(!is.null(distance_matrix)){
    dist.vect <- distance_matrix %>% unlist()
    rug(dist.vect/max(dist.vect),
        side=1)
  }
    
}#END FUNCTION

# median(dist.vect/max(dist.vect))