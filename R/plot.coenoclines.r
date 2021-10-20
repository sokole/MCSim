#' @title plot standardized dispersal kernel
#' 
#' @description plot standardized dispersal kernel
#' 
#' @usage 
#' plot.coenoclines(sim.result)
#' 
#' @param sim.result output from a simulation, default is NULL
#' @param landscape a landscape object, default is NULL
#' @param Ef numeric vector of environmental filter values for sites in a simulation result, default is NULL
#' @param trait.Ef numeric vector of niche positions of species in metacommunity, default is NULL
#' @param trait.Ef.sd numeric vector of niche breadths, default is NULL
#' @param plot.trait.Ef.sd.0s should species with niche-breadth of 0 be plotted? Default is FALSE
#' @param y.max.val (optional) max value for y-axis
#' 
#' 
#' @export
#' 
plot.coenoclines <- function(
    sim.result = NULL,
    landscape = NULL,
    Ef = NULL,
    trait.Ef = NULL,
    trait.Ef.sd = NULL,
    plot.trait.Ef.sd.0s = FALSE,
    y.max.val = NULL,
    ...){
    
  # check for necessary info
  if(is.null(trait.Ef)){
    
    if(!is.null(sim.result)){
      if("dat.gamma.t0" %in% names(sim.result)){
        trait.Ef <- sim.result$dat.gamma.t0$trait.Ef
      }else if("dat.gamma.t0.list" %in% names(sim.result)){
        trait.Ef <- sim.result$dat.gamma.t0.list[[1]]$dat.gamma.t0$trait.Ef
      }
    }else{
      stop('no value provided for trait.Ef')
    }
  }
  
  # check for necessary info
  if(is.null(trait.Ef.sd)){
    if(!is.null(sim.result)){
      
      if("dat.gamma.t0" %in% names(sim.result)){
        trait.Ef.sd <- sim.result$dat.gamma.t0$trait.Ef.sd
      }else if("dat.gamma.t0.list" %in% names(sim.result)){
        trait.Ef.sd <- sim.result$dat.gamma.t0.list[[1]]$dat.gamma.t0$trait.Ef.sd
      }
    }else{
      stop('no value provided for trait.Ef.sd')
    }
  }
  
  # find Ef data in landscape or sim.result
  if(is.null(Ef)){
    if(!is.null(landscape)){
      Ef <- landscape$site.info$Ef
    }else if(!is.null(sim.result)){
      if("landscape" %in% names(sim.result)){
        site.info <- sim.result$landscape$site.info %>%
          dplyr::mutate(Ef.rank = rank(.data$Ef))
      }else if("landscape.list" %in% names(sim.result)){
        site.info <- sim.result$landscape.list[[1]]$site.info %>%
          dplyr::mutate(Ef.rank = rank(.data$Ef))
        message("WARNING: this sim result includes a changing landscape, this plotting function is has not been optimized for changing landscapes and the plot will be based on the initial landscape configuration")
      }
      Ef <- site.info$Ef
    }
  }
  
  

  n.spp <- length(trait.Ef)

  # -- function for plotting bell curves
  fn.norm.curve <- function(sigma=1, mu = 0,...) {
    curve(
      (1/sigma * sqrt(2 * pi)) * exp((-1  *(x - mu)^2) / (2 * sigma^2)), ...) #formula for bell curve
  }

  # # get y.max.val, minimum is 1
  if(is.null(y.max.val)){
    y.max.val <- 1
    for (i.spp in 1:n.spp){
      if(trait.Ef.sd[i.spp] > 0){
        max.val.tmp <- (fn.norm.curve(
          mu = trait.Ef[i.spp], 
          sigma = trait.Ef.sd[i.spp],
          from = min(Ef),
          to = max(Ef),
          add = FALSE))$y %>% max()
        y.max.val <- max(c(y.max.val, max.val.tmp), na.rm = TRUE)
      }
    }
  }

  
  # -- Initialize plot of coenoclines

  plot(1,1,
       xlim = c(min(Ef), max(Ef)),
       ylim = c(0, y.max.val),
       type = 'n',
       xlab = 'Environmental gradient',
       ylab = 'Prob. dens.',
       main = 'Niche positions and widths',
       ...)
  
  mypal <- grDevices::rainbow(n.spp)
  
  

  # -- loop to plot each species' habitat preference
  for (i.spp in 1:n.spp){
    if(trait.Ef.sd[i.spp] > 0){
      fn.norm.curve(
        mu = trait.Ef[i.spp], 
        sigma = trait.Ef.sd[i.spp], 
        add = TRUE,
        col = mypal[i.spp])
    }else if(plot.trait.Ef.sd.0s & trait.Ef.sd[i.spp] == 0){
      # plot verticle lines when trait.Ef.sd == 0, if desired
      abline(v = trait.Ef[i.spp],
             col = mypal[i.spp])
    }
  }

  
  # -- plot sites along the x-axis
  rug(Ef)
}
