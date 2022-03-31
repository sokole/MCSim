#' @title plot dot plots
#' 
#' @description plot dot plots of selected timesteps from an MCSim simulation
#' 
#' @usage 
#' plot.dot.plots(sim.result)
#' 
#' @param sim.result output from a simulation
#' @param timesteps (array of integers) timesteps to plot. Default is "NA" and will auto select timesteps
#' @param spp.label.threshold.RA (numeric) the relative abundance (RA) cutoff for taxa to label on the plot. Devault is to label taxa with a max RA >= 0.40
#' 
#' @export plot.dot.plots
#' 
#' @examples 
#' \dontrun{
#' plot.dot.plots(my_sim_result) 
#' }
#' 
plot.dot.plots <- function(
  sim.result,
  timesteps = NA_integer_,
  spp.label.threshold.RA = 0.40
  ){
  
  # clean up timesteps argument
  if(all(is.na(timesteps))){
    timesteps <- c(1,2,10, 100, max(sim.result$J.long$timestep)) %>% unique() %>% stats::na.omit()
  }else{
    timesteps <- timesteps %>% unique() %>% stats::na.omit()
  }
  
  # make sure timesteps exist, only use those in the dataset
  timesteps <- timesteps %>% dplyr::intersect(unique(sim.result$J.long$timestep)) 
  
  
  # -- extract timesteps to plot
  J.long <- dplyr::filter(sim.result$J.long, .data$timestep %in% timesteps)
  
  # -- group spp counts by time and site to calculate site count totals
  J.time.site <- dplyr::group_by(J.long, .data$timestep, .data$site)
  JLs <- dplyr::summarize(J.time.site, site.totals=sum(.data$count))
  J.JLs <- dplyr::full_join(J.long, JLs, by=c('timestep','site'))
  
  # -- calculate relative abundances (RAs) from counts and site totals, remove
  # observations with RAs of 0
  J.RAs.long <- J.JLs
  J.RAs.long$RA <- J.RAs.long$count/J.RAs.long$site.totals
  J.RAs.long <- dplyr::filter(J.RAs.long, .data$RA > 0)
  
  # -- add environmental data and species names for plotting
  J.RAs.long <- dplyr::mutate(
    J.RAs.long, 
    spp.no = as.numeric(as.factor(.data$spp)))
  
  # extract site info from the landscape object
  # rank sites locations along the env gradient
  
  # browser()
  
  if("landscape" %in% names(sim.result)){
    site.info <- sim.result$landscape$site.info %>%
      dplyr::mutate(Ef.rank = rank(.data$Ef,ties.method = "first"))
  }else if("landscape.list" %in% names(sim.result)){
    site.info <- sim.result$landscape.list[[1]]$site.info %>%
      dplyr::mutate(Ef.rank = rank(.data$Ef,ties.method = "first"))
    message("WARNING: this sim result includes a changing landscape, this plotting function is has not been optimized for changing landscapes and the plot will be based on the initial landscape configuration")
  }
  
  

  # extract starting regional species pool from sim object
  # rank spp positions along env gradient
  
  if("dat.gamma.t0" %in% names(sim.result)){
    spp.info <- sim.result$dat.gamma.t0 %>%
      dplyr::mutate(trait.rank = rank(trait.Ef,ties.method = "first"))
  }else if("dat.gamma.t0.list" %in% names(sim.result)){
    spp.info <- sim.result$dat.gamma.t0.list[[1]]$dat.gamma.t0 %>%
      dplyr::mutate(trait.rank = rank(.data$trait.Ef,ties.method = "first"))
  }

  
  # J.RAs.long$region <- d.site.1D[J.RAs.long$site, 'region']
  J.RAs.long$Ef.rank <- site.info[J.RAs.long$site, 'Ef.rank']
  J.RAs.long$trait.rank <- spp.info[J.RAs.long$spp,'trait.rank']
  
  # -- make a species characteristic data frame for labeling the spp axis
  d.spp.RAs <- as.data.frame(
    J.RAs.long %>% 
      dplyr::group_by(.data$spp) %>% 
      dplyr::summarize(max.RA = max(.data$RA),
                       spp.no = .data$spp.no[1],
                       trait.rank = .data$trait.rank[1]))
  
  # -- only label species with RAs over 0.40
  spp.labels <- dplyr::filter(d.spp.RAs, max.RA >= spp.label.threshold.RA)
  
  # -- make plot
  p <- ggplot2::ggplot(
    J.RAs.long, 
    ggplot2::aes(
      as.factor(.data$trait.rank), 
      as.factor(.data$Ef.rank), 
      size = .data$RA))
  
  browser()
  
  p <- p + 
    ggplot2::geom_point() + 
    ggplot2::facet_grid(. ~ .data$timestep) +
    ggplot2::scale_size('Relative\nAbundance', range = c(.5,4)) +
    ggplot2::labs(
      title = paste(
        stats::na.omit(c(
        "Metacommunity Simulation",
        sim.result$scenario.ID,
        "\nNo. generations ( ---> )")), collapse = " "),
         y = 'Site ID',
         x = 'Species ID') +
    ggplot2::scale_x_discrete(
      breaks = spp.labels$trait.rank,
      labels = spp.labels$spp) +
    ggplot2::theme_bw() +
    ggplot2::coord_flip()
  
  return(p)
}