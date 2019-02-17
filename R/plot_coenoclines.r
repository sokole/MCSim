#' @title plot standardized dispersal kernel
#' 
#' @description plot standardized dispersal kernel
#' 
#' @param sim_result output from a simulation, default is NULL
#' @param landscape a landscape object, default is NULL
#' @param Ef numeric vector of environmental filter values for sites in a simulation result, default is NULL
#' @param trait_Ef numeric vector of niche positions of species in metacommunity, default is NULL
#' @param trait_Ef_sd numeric vector of niche breadths, default is NULL
#' 
#' 
#' 
#' @export
plot_coenoclines <- function(
    sim_result = NULL,
    landscape = NULL,
    Ef = NULL,
    trait_Ef = NULL,
    trait_Ef_sd = NULL,
    ...){
    
  requireNamespace("magrittr")
  
  # check for necessary info
  if(is.null(trait_Ef)){
    if(!is.null(sim_result)){
      trait_Ef <- sim_result$dat.gamma.t0$trait.Ef
    }else{
      stop('no value provided for trait_Ef')
    }
  }
  
  # check for necessary info
  if(is.null(trait_Ef_sd)){
    if(!is.null(sim_result)){
      trait_Ef_sd <- sim_result$dat.gamma.t0$trait.Ef.sd
    }else{
      stop('no value provided for trait_Ef_sd')
    }
  }
  
  # find Ef data in landscape or sim_result
  if(is.null(Ef)){
    if(!is.null(landscape)){
      Ef <- landscape$site.info$Ef
    }else if(!is.null(sim.result)){
      Ef <- sim.result$landscape$site.info$Ef
    }
  }
  
  n_spp <- length(trait_Ef)

  # -- function for plotting bell curves
  fn_norm_curve <- function(sigma=1, mu = 0,...) {
    curve(
      (1/sigma * sqrt(2 * pi)) * exp((-1  *(x - mu)^2) / (2 * sigma^2)), ...) #formula for bell curve
  }

  # # get y_max_val
  y_max_val <- 1

  for (i.spp in 1:n_spp){
    max_val_tmp <- (fn_norm_curve(
      mu = trait_Ef[i.spp], 
      sigma = trait_Ef_sd[i.spp], 
      add = FALSE))$y %>% max()
    y_max_val <- max(c(y_max_val, max_val_tmp), na.rm = TRUE)
  }

  # -- Initialize plot of coenoclines
  plot(1,1,
       xlim = c(min(Ef), max(Ef)),
       ylim = c(0, y_max_val),
       type = 'n',
       xlab = 'Environmental gradient',
       ylab = 'Prob. dens.',
       main = 'Niche positions and widths',
       ...)

  mypal <- rainbow(n_spp)

  # -- loop to plot each species' habitat preference
  for (i.spp in 1:nrow(sim_result$dat.gamma.t0)){
    fn_norm_curve(
      mu = trait_Ef[i.spp], 
      sigma = trait_Ef_sd[i.spp], 
      add = TRUE,
      col = mypal[i.spp])
  }

  # -- plot sites along the x-axis
  rug(Ef)
}
