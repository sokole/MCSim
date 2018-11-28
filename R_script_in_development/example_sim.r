# setwd("C:/Users/esokol/Documents/Git/MCSim/R_script_in_development")

library(tidyverse)

source('fn.metaSIM.disturbance.R')

# full time series name, string together pre and post disturbance together using this character string
my_scenario_name <- 'mySim'

# generating pre and post disturbance landscapes
xy_coordinates <- data.frame(
  x = c(1, 2, 3, 4, 5),
  y = c(1, 3, 1, 5, 2))

#############
# send function landscape list and interval duration list, along with other metaSim parameters that will be fixed across all simulations

my_landscape <- MCSim::fn.make.landscape( # pre-disturbance landscape
      site.coords = xy_coordinates,
      m = c(0.5, 0.5, 0.1, 0.1, 0.01),
      Ef = c(-1, -.25, .1, 1, 2),
      JL = c(100, 100, 100, 100, 100))


# -- fixed across simulations
# niche positions, niche breadths, and relative abundances for three species
my_niche_positions <-  rep(.1, 6)
my_niche_breadths <- rep(.1, 6)
my_regional_rel_abund <- rep(.1, 6)

# -- call the wrapper function
my_sim_result <- fn.metaSIM(
  scenario.ID = my_scenario_name,
  landscape = my_landscape,
  gamma.abund = my_regional_rel_abund,
  trait.Ef = my_niche_positions,
  trait.Ef.sd = my_niche_breadths,
  W.r = 0,
  nu = 0.001,
  save.sim = FALSE
)

data_wide <- my_sim_result$J.long %>% group_by(timestep, site) %>% spread(spp, count)
print(data_wide)


my_sim_result$J.long %>% ggplot(aes(timestep, count, color = spp)) +
  geom_line() +
  facet_wrap(~ as.factor(site))
