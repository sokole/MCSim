#' @title  simple MCSim scenarios for Ben
#' @author Eric Sokol (sokole@gmail.com)
#' @description  Distinct scenarios -- 
#'   1. species sorting with no dispersal (everything everywhere but environment selects)
#'   2. species sorting with strong dispersal limitation
#'   3. neutral no dispersal limitation
#'   4. neutral with dispersal limitation

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# required packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

# use devtools package to install MCSim from github
# devtools::install_github('sokole/MCSim')
library(MCSim)

library(tidyverse)

my_scenario_name <- 'neutral_unlimited_dispersal'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# landscape design common to all sims -- two sets of clumped sites
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(123)

# generating pre and post disturbance landscapes
xy_coordinates <- data.frame(
  x = c(1, 2, 3, 8, 9, 10),
  y = c(1, 2, 1, 8, 9, 8))

plot(xy_coordinates)

#############
# send function landscape list and interval duration list, along with other metaSim parameters that will be fixed across all simulations
n_sites <- nrow(xy_coordinates)

my_landscape <- MCSim::fn.make.landscape( # pre-disturbance landscape
  site.coords = xy_coordinates,
  m = rep(.1, n_sites),
  Ef = runif(n_sites, -1, 1),
  JL = rep(100, n_sites))


# -- fixed across simulations
# niche positions, niche breadths, and relative abundances for 6 species
my_regional_rel_abund <- c(9, 2, 1, 1)
my_regional_rel_abund <- my_regional_rel_abund / sum(my_regional_rel_abund)

n_spp <- length(my_regional_rel_abund)
my_niche_positions <-  rep(.1, n_spp)
my_niche_breadths <- rep(.1, n_spp)


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

