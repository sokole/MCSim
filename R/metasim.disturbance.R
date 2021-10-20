#############
# send function landscape list and interval duration list, along with other metaSim parameters that will be fixed across all simulations

#' @title A metacommunity simulation for ecologists
#' 
#' @aliases fn.metaSIM.disturbance
#' 
#' @usage 
#' metasim.disturbance(scenario.name, landscape.list, time.interval.durations)
#' 
#' @description Wrapper function for MCSim::metasim, sends landscape list and interval duration list, along with other metasim parameters that will be fixed across all simulations
#' 
#' @param scenario.name text string naming the scenario
#' @param landscape.list list of landscapes to run in series
#' @param time.interval.durations vector of integers representing the number of timesteps to use for each landscape in the landscape_list, MUST be same length as landscape_list
#' @param J.t0.initial Initial species abundance matrix for timestep "t0"
#' 
#' @param scenario_name DEPRECATED: use \code{scenario.name}
#' @param landscape_list DEPRECATED: use \code{landscape.list}
#' @param time_interval_durations DEPRECATED: use \code{time.interval.durations}
#' @param J.t0_initial Deprecated: use \code{J.t0.initial}

#' @param ... parameters sent to metasim()
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' set.seed(1234) #set random seed
#' 
#' # make a landscapes
#' my.landscape.1 <- make.landscape(
#'   site.coords = data.frame(
#'     x = c(1, 2, 3, 4, 5),
#'     y = c(1, 3, 1, 5, 2)),
#'   Ef = c(-.8, -.6, 0, .25, .9),
#'   m = 0.5,
#'   JM = 10000)
#' 
#' my.landscape.2 <- make.landscape(
#'   site.coords = data.frame(
#'     x = c(1, 2, 3, 4, 5),
#'     y = c(1, 3, 1, 5, 2)),
#'   Ef = c(-.8, -.6, 0, 0, 0),
#'   m = 0.5,
#'   JM = 10000)
#' 
#' # niche positions, niche breadths, and relative abundances for three species
#' niche.positions <-  c(-.5, 0, .5)
#' niche.breadths <- c(.2, .2, 5)
#' regional.rel.abund <- c(.8, .1, .1)
#' 
#' sim.result <- metasim.disturbance(
#'   scenario.name = "niche_shift",
#'   landscape.list = list(my.landscape.1,my.landscape.2),
#'   time.interval.durations = c(10,10),
#'   trait.Ef = niche.positions,
#'   trait.Ef.sd = niche.breadths,
#'   gamma.abund = regional.rel.abund,
#'   W.r = 0,
#'   nu = 0.001,
#'   output.dir.path = "my_disturbance_sim_output_directory")
#' }
metasim.disturbance <- function(
  
  scenario.name, # text string
  landscape.list, # list of landscapes to run in series
  time.interval.durations, # vector of integers representing the number of timesteps to use for each landscape in the landscape_list, same length as landscape_list
  J.t0.initial = NULL,
  
  # Deprecated arguments
  scenario_name = NULL, # text string
  landscape_list = NULL, # list of landscapes to run in series
  time_interval_durations = NULL, # vector of integers representing the number of timesteps to use for each landscape in the landscape_list, same length as landscape_list
  J.t0_initial = NULL,
  
  ...){
  
  if(!is.null(scenario_name)){
    scenario.name <- scenario_name 
    warning('Input argument "scenario_name" is deprecated, use "scenario.name" instead', 
            call. = FALSE)    
  }
  if(!is.null(landscape_list)){
    landscape.list <- landscape_list
    warning('Input argument "landscape_list" is deprecated, use "landscape.list" instead', 
            call. = FALSE)    
  }
  if(!is.null(time_interval_durations)){
    time.interval.durations <- time_interval_durations
    warning('Input argument "time_interval_durations" is deprecated, use "time.interval.durations" instead', 
            call. = FALSE)    
  }
  if(!is.null(J.t0_initial)){
    J.t0.initial <- J.t0_initial
    warning('Input argument "J.t0_initial" is deprecated, use "J.t0.initial" instead', 
            call. = FALSE)    
  }
  
  my_scenario_ID <- paste0(scenario.name,'_', Sys.time() %>% format('%Y%m%d_%H%M%S'))
  # A data frame with coordinates for 5 sites
  
  sim_result_list <- list()
  ###################
  for (i_ts in 1:length(landscape.list)){
    
    #################
    if(!is.null(names(landscape.list))){
      sim_id_name <- names(landscape.list)[i_ts]
    }else{
      sim_id_name <- paste0('TS',i_ts)
    }
    
    if(i_ts == 1){
      
      if(is.null(J.t0.initial)){
        sim_result_list[[i_ts]] <- MCSim::metasim(
          landscape = landscape.list[[i_ts]],
          n.timestep = time.interval.durations[i_ts], 
          scenario.ID = my_scenario_ID,
          sim.ID = sim_id_name,
          ...) 
      }else{
        sim_result_list[[i_ts]] <- MCSim::metasim(
          landscape = landscape.list[[i_ts]],
          n.timestep = time.interval.durations[i_ts], 
          scenario.ID = my_scenario_ID,
          sim.ID = sim_id_name,
          J.t0 = J.t0.initial,
          ...) 
      }
      
    }else if(i_ts > 1){
      J_last_timestep_long <- sim_result_list[[i_ts - 1]]$J.long %>% 
        dplyr::filter(timestep == max(.data$timestep)) 
      
      last_time_step_number <- sim_result_list[[i_ts - 1]]$J.long$timestep %>% max()
      
      J_last_timestep_wide <- J_last_timestep_long %>%
        dplyr::select(.data$site, .data$spp, .data$count) %>%
        dplyr::group_by(site) %>%
        dplyr::ungroup() %>%
        # dplyr::spread(spp, count, fill = 0) %>%
        tidyr::pivot_wider(id_cols = .data$site, names_from = .data$spp,values_from = .data$count) %>%
        tibble::column_to_rownames('site') %>%
        as.data.frame()
      
      sim_result_list[[i_ts]] <- MCSim::metasim(
        landscape = landscape.list[[i_ts]],
        J.t0 = J_last_timestep_wide,
        n.timestep = time.interval.durations[i_ts], 
        scenario.ID = my_scenario_ID,
        sim.ID = sim_id_name,
        ...) 
      
      # modify timestep to be global for combined time series
      sim_result_list[[i_ts]]$J.long <- sim_result_list[[i_ts]]$J.long %>% 
        dplyr::mutate(timestep = .data$timestep + last_time_step_number)
    }
    
    # add simulation info to link output to metadata
    sim_result_list[[i_ts]]$J.long <- sim_result_list[[i_ts]]$J.long %>% 
      dplyr::mutate(
        scenario.ID = sim_result_list[[i_ts]]$scenario.ID,
        sim.result.name = sim_result_list[[i_ts]]$sim.result.name)
  }
  
  return(
    list(
      scenario.ID = sim_result_list[[1]]$scenario.ID,
      sim.result.name.list = sim_result_list %>% purrr::map(`[`, "sim.result.name"),
      landscape.list = landscape.list,
      dat.gamma.t0.list = sim_result_list %>% purrr::map(`[`, "dat.gamma.t0"),
      W.r.list = sim_result_list %>% purrr::map(`[`, "W.r"),
      J.long = sim_result_list %>% purrr::map(`[`, "J.long") %>%
        purrr::flatten() %>% 
        dplyr::bind_rows()
    ))
} #END FUNCTION


# define alias functions
#' @export
fn.metaSIM.disturbance <- metasim.disturbance