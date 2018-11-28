# setwd("C:/Users/esokol/Documents/Git/MCSim/R_script_in_development")

library(tidyverse)

source('fn.metaSIM.disturbance.R')

# full time series name, string together pre and post disturbance together using this character string
my_scenario_name <- 'mySim'

# generating pre and post disturbance landscapes
xy_coordinates <- data.frame(
  x = c(1, 1.5, 3, 4, 5),
  y = c(1, 1.5, 1, 6, 2))

plot(xy_coordinates)

#############
# send function landscape list and interval duration list, along with other metaSim parameters that will be fixed across all simulations

my_landscape_list <- 
  list(
    landscape_1 = MCSim::fn.make.landscape( # pre-disturbance landscape
      site.coords = xy_coordinates,
      m = c(0.5, 0.5, 0.1, 0.1, 0.01),
      Ef = c(-1, -.25, .1, 1, 2),
      JL = c(100, 100, 100, 100, 100)),
    landscape_2 = MCSim::fn.make.landscape( # disturbed landscape with reduced carrying capacity
      site.coords = xy_coordinates,
      m = c(0.5, 0.5, 0.1, 0.1, 0.01),
      Ef = c(-1, -.25, .1, 1, 2),
      JL = c(50, 50, 100, 100, 100)), # <-- decreased carrying capacity at 2 sites
    landscape_3 = MCSim::fn.make.landscape(
      site.coords = xy_coordinates,
      m = c(0.5, 0.5, 0.1, 0.1, 0.01),
      Ef = c(-1, -.25, .1, 1, 2),
      JL = c(100, 100, 100, 100, 100))) # back to pre-disturbance conditions

my_time_interval_durations <- c(10, 3, 10 )

# -- fixed across simulations
# niche positions, niche breadths, and relative abundances for three species
my_niche_positions <-  c(-.25, 0, .25)
my_niche_breadths <- c(.2, 1, 1)
my_regional_rel_abund <- c(.8, .1, .1)
trait.dispersal <- c(1, 50, 100)

# -- call the wrapper function
my_sim_result <- fn.metaSIM.disturbance(
  ###### below are the parameters that are new to this wrapper function
  scenario_name = my_scenario_name,
  landscape_list = my_landscape_list,
  time_interval_durations = my_time_interval_durations,
  ####### below is stuff that would be similar to a regular fn.metaSIM() function call
  gamma.abund = my_regional_rel_abund,
  trait.Ef = my_niche_positions,
  trait.Ef.sd = my_niche_breadths,
  W.r = 0,
  nu = 0.001,
  save.sim = FALSE
)

my_sim_result$J.long %>% ggplot(aes(timestep, count, color = spp)) +
  geom_line() +
  facet_wrap(~ as.factor(site))

data_wide <- my_sim_result$J.long %>% group_by(timestep, site) %>% spread(spp, count)