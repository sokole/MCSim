#############
# send function landscape list and interval duration list, along with other metaSim parameters that will be fixed across all simulations

#' @title A metacommunity simulation for ecologists
#' 
#' @aliases fn.metaSIM.disturbance
#' 
#' @usage 
#' metasim.disturbance(scenario_name, landscape_list, time_interval_durations)
#' 
#' @description Wrapper function for MCSim::metasim, sends landscape list and interval duration list, along with other metasim parameters that will be fixed across all simulations
#' 
#' @param scenario_name text string naming the scenario
#' @param landscape_list list of landscapes to run in series
#' @param time_interval_durations vector of integers representing the number of timesteps to use for each landscape in the landscape_list, MUST be same length as landscape_list
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
#'   scenario_name = "niche_shift",
#'   landscape_list = list(my.landscape.1,my.landscape.2),
#'   time_interval_durations = c(10,10),
#'   trait.Ef = niche.positions,
#'   trait.Ef.sd = niche.breadths,
#'   gamma.abund = regional.rel.abund,
#'   W.r = 0,
#'   nu = 0.001,
#'   output.dir.path = "my_disturbance_sim_output_directory")
#' }
metasim.disturbance <- function(
  
  scenario_name, # text string
  landscape_list, # list of landscapes to run in series
  time_interval_durations, # vector of integers representing the number of timesteps to use for each landscape in the landscape_list, same length as landscape_list
  J.t0_initial = NULL,
  ...){

my_scenario_ID <- paste0(scenario_name,'_', Sys.time() %>% format('%Y%m%d_%H%M%S'))
# A data frame with coordinates for 5 sites

sim_result_list <- list()
###################
for (i_ts in 1:length(landscape_list)){
  
  #################
  if(!is.null(names(landscape_list))){
    sim_id_name <- names(landscape_list)[i_ts]
  }else{
    sim_id_name <- paste0('TS',i_ts)
  }
  
  if(i_ts == 1){
    
    if(is.null(J.t0_initial)){
      sim_result_list[[i_ts]] <- MCSim::metasim(
        landscape = landscape_list[[i_ts]],
        n.timestep = time_interval_durations[i_ts], 
        scenario.ID = my_scenario_ID,
        sim.ID = sim_id_name,
        ...) 
    }else{
      sim_result_list[[i_ts]] <- MCSim::metasim(
        landscape = landscape_list[[i_ts]],
        n.timestep = time_interval_durations[i_ts], 
        scenario.ID = my_scenario_ID,
        sim.ID = sim_id_name,
        J.t0 = J.t0_initial,
        ...) 
    }

  }else if(i_ts > 1){
    J_last_timestep_long <- sim_result_list[[i_ts - 1]]$J.long %>% 
      filter(timestep == max(timestep)) 
    
    last_time_step_number <- sim_result_list[[i_ts - 1]]$J.long$timestep %>% max()
    
    J_last_timestep_wide <- J_last_timestep_long %>%
      select(site, spp, count) %>%
      group_by(site) %>%
      ungroup() %>%
      spread(spp, count, fill = 0) %>%
      column_to_rownames('site') %>%
      as.data.frame()
    
    sim_result_list[[i_ts]] <- MCSim::metasim(
      landscape = landscape_list[[i_ts]],
      J.t0 = J_last_timestep_wide,
      n.timestep = time_interval_durations[i_ts], 
      scenario.ID = my_scenario_ID,
      sim.ID = sim_id_name,
      ...) 
    
    # modify timestep to be global for combined time series
    sim_result_list[[i_ts]]$J.long <- sim_result_list[[i_ts]]$J.long %>% mutate(timestep = timestep + last_time_step_number)
  }
  
  # add simulation info to link output to metadata
  sim_result_list[[i_ts]]$J.long <- sim_result_list[[i_ts]]$J.long %>% mutate(
    scenario.ID = sim_result_list[[i_ts]]$scenario.ID,
    sim.result.name = sim_result_list[[i_ts]]$sim.result.name)
}

return(
  list(
    scenario.ID = sim_result_list[[1]]$scenario.ID,
    sim.result.name.list = sim_result_list %>% purrr::map(`[`, "sim.result.name"),
    landscape.list = landscape_list,
    dat.gamma.t0.list = sim_result_list %>% purrr::map(`[`, "dat.gamma.t0"),
    W.r.list = sim_result_list %>% purrr::map(`[`, "W.r"),
    J.long = sim_result_list %>% purrr::map(`[`, "J.long") %>%
      flatten() %>% bind_rows()
  ))
} #END FUNCTION


# define alias functions
#' @export
fn.metaSIM.disturbance <- metasim.disturbance