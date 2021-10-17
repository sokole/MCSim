#' @title plot dot plots
#' 
#' @description plot dot plots of selected timesteps from an MCSim simulation
#' 
#' @param sim_result output from a simulation
#' @param timesteps (array of integers) timesteps to plot. Default is "NA" and will auto select timesteps
#' @param spp_label_threshold_RA (numeric) the relative abundance (RA) cutoff for taxa to label on the plot. Devault is to label taxa with a max RA >= 0.40
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' plot_dot_plots(my_sim_result) 
#' }
plot_dot_plots <- function(
  sim_result,
  timesteps = NA_integer_,
  spp_label_threshold_RA = 0.40
  ){
  
  # clean up timesteps argument
  if(all(is.na(timesteps))){
    timesteps <- c(1,2,10, 100, max(sim_result$J.long$timestep)) %>% unique() %>% stats::na.omit()
  }else{
    timesteps <- timesteps %>% unique() %>% stats::na.omit()
  }
  
  # make sure timesteps exist, only use those in the dataset
  timesteps <- timesteps %>% dplyr::intersect(unique(sim_result$J.long$timestep)) 
  
  
  # -- extract timesteps to plot
  J.long <- dplyr::filter(sim_result$J.long, .data$timestep %in% timesteps)
  
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
  if("landscape" %in% names(sim_result)){
    site_info <- sim_result$landscape$site.info %>%
      dplyr::mutate(Ef_rank = rank(.data$Ef))
  }else if("landscape.list" %in% names(sim_result)){
    site_info <- sim_result$landscape_list[[1]]$site.info %>%
      dplyr::mutate(Ef_rank = rank(.data$Ef))
    message("Your sim result includes a changing landscape, this plotting function is has not been optimized for changing landscapes and the plot will be based on the initial landscape configuration")
  }
  
  
  # extract starting regional species pool from sim object
  # rank spp positions along env gradient
  
  if("dat.gamma.t0" %in% names(sim_result)){
    spp_info <- sim_result$dat.gamma.t0 %>%
      mutate(trait_rank = rank(trait.Ef))
  }else if("dat.gamma.t0.list" %in% names(sim_result)){
    spp_info <- sim_result$dat.gamma.t0.list[[1]]
      mutate(trait_rank = rank(trait.Ef))
  }

  
  # J.RAs.long$region <- d.site.1D[J.RAs.long$site, 'region']
  J.RAs.long$Ef_rank <- site_info[J.RAs.long$site, 'Ef_rank']
  J.RAs.long$trait_rank <- spp_info[J.RAs.long$spp,'trait_rank']
  
  # -- make a species characteristic data frame for labeling the spp axis
  d.spp.RAs <- as.data.frame(
    J.RAs.long %>% 
      dplyr::group_by(.data$spp) %>% 
      dplyr::summarize(max.RA = max(.data$RA),
                       spp.no = .data$spp.no[1],
                       trait_rank = .data$trait_rank[1]))
  
  # -- only label species with RAs over 0.40
  spp.labels <- dplyr::filter(d.spp.RAs, max.RA >= spp_label_threshold_RA)
  
  # -- make plot
  p <- ggplot2::ggplot(
    J.RAs.long, 
    ggplot2::aes(
      as.factor(.data$trait_rank), 
      as.factor(.data$Ef_rank), 
      size = .data$RA))
  
  p <- p + 
    ggplot2::geom_point() + 
    ggplot2::facet_grid(. ~ .data$timestep) +
    ggplot2::scale_size('Relative\nAbundance', range = c(.5,4)) +
    ggplot2::labs(
      title = 'Simulated Metacommunity Simulation\nNo. generations ( ---> )',
         y = 'Site ID',
         x = 'Species ID') +
    ggplot2::scale_x_discrete(
      breaks = spp.labels$trait_rank,
      labels = spp.labels$spp) +
    ggplot2::theme_bw() +
    ggplot2::coord_flip()
  
  return(p)
}