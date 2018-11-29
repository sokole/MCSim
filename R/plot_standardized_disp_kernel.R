plot_standardized_disp_kernel <- function(
  Wr = NULL,
  sim_result = NULL,
  landscape = NULL,
  distance_matrix = NULL,
  ...){
  
  w <- max(c(sim_result$W.r, Wr))
  
  if(!is.null(sim_result)) Wr <- c(w, Wr) %>% unique()
  
  
  stopifnot({length(Wr) > 0})
  
  curve(exp(-1*w*x^2),0,1,
        ylim=c(0,1.05),
        ylab='W(r)',
        xlab='r')
  
  if(length(Wr) > 1){
    for(i in 2:length(Wr)){
      
      curve(exp(-1*Wr[i]*x^2),0,1,
            ylim = c(0,1.15),
            lty = i,
            lwd = 1.25,
            add = TRUE)
    }
    
    legend(x = 'topright',
           legend = Wr[order(Wr)],
           lty = c(1:length(Wr))[order(Wr)],
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