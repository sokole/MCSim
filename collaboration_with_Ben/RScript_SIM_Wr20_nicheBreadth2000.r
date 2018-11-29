#' @title  simple MCSim scenarios for Ben
#' @author Eric Sokol (sokole@gmail.com)
#' @description  Distinct scenarios -- 
#'   1. species sorting with no dispersal (everything everywhere but environment selects)
#'   2. species sorting with strong dispersal limitation
#'   3. neutral no dispersal limitation
#'   --> 4. neutral with dispersal limitation

rm(list=ls())
gc()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# required packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

# use devtools package to install MCSim from github
# devtools::install_github('sokole/MCSim')
library(MCSim)

library(tidyverse)


# user provided vars
my_W.r <- 20 #dispersal kernel slope, 0 = no limitation

my_niche_breadth_multiplier <- 2000 #increases neutrality by increasing niche bredth, numbers > 2 should be approx neutral

my_scenario_name <- paste0('Wr',my_W.r,'_nicheBreadth',my_niche_breadth_multiplier)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# arguments common to all sims -- two sets of clumped sites
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(123)

# default is no invasion
# my_nu <- 0.001 #metacommunity invasion probability
# my_speciation.limit <- 5 #set limit to number of potential invaders
my_nu <- 0
my_speciation.limit <- 0 #closed system with no invaders

# number of time steps
my_n.timestep <- 200

# -- fixed across simulations -- metacommunity characteristics
# niche positions, niche breadths, and relative abundances for 6 species
my_J.t0 <- data.frame(
  spp1 = c(48, 48, 48,
           48, 48, 48,
           1, 0, 1,
           0, 1, 0),
  spp2 = c(1, 0, 1,
           48, 48, 48,
           48, 48, 48,
           0, 1, 0),
  spp3 = c(5, 5, 5,
           5, 5, 5,
           5, 5, 5,
           5, 5, 5),
  spp4 = c(48, 48, 48,
           0, 1, 0,
           1, 0, 1,
           48, 48, 48),
  spp5 = c(0, 1, 0,
           1, 0, 1,
           48, 48, 48,
           48, 48, 48))

my_niche_positions <-  c(.2, -.5, 0, -.15, .5)
my_niche_breadths <- c(.15, .15, .3, .15, .15) * my_niche_breadth_multiplier

# -- fixed across simulations -- landscape characteristics
# defining site data for sim landscape
site_data <- data.frame(
  x = c(1, 2, 3, 
        1, 2, 3, 
        11, 12, 13,
        11, 12, 13),
  y = c(1, 2, 1, 
        11, 12, 11,
        1, 2, 1,
        11, 12, 11),
  Ef = c(.1, .2, .3,
         -.1, -.2, -.3,
         -.5, .5, -.9,
         .2, .3, .2
         ),
  m = rep(.1, 12),
  JL = rowSums(my_J.t0))


#############
# send function landscape list and interval duration list, along with other metaSim parameters that will be fixed across all simulations

my_landscape <- MCSim::fn.make.landscape(site_data)

# -- call the wrapper function
my_sim_result <- fn.metaSIM(
  scenario.ID = my_scenario_name,
  landscape = my_landscape,
  J.t0 = my_J.t0,
  trait.Ef = my_niche_positions,
  trait.Ef.sd = my_niche_breadths,
  W.r = my_W.r,
  nu = my_nu,
  speciation.limit = my_speciation.limit,
  n.timestep = my_n.timestep,
  save.sim = TRUE
)

data_wide <- my_sim_result$J.long %>% group_by(timestep, site) %>% spread(spp, count)
print(data_wide)

# plot dispersal kernel
file_name <- paste0('SIM_OUTPUT/Fig_disp_kern_',my_scenario_name,'.pdf')
pdf(file_name)
plot_standardized_disp_kernel(my_W.r, landscape = my_landscape)
dev.off()

# plot sites
file_name <- paste0('SIM_OUTPUT/Fig_landscape_',my_scenario_name,'.pdf')
pdf(file_name)
site_data %>%
  ggplot(aes(x, y,
             size = JL,
             color = Ef)) +
  geom_point() +
  scale_color_gradient2(
    low = "red",
    mid = 'gray',
    high = "darkblue",
    midpoint = 0)
dev.off()


# plot dispersal kernel
file_name <- paste0('SIM_OUTPUT/Fig_comm_TS_',my_sim_result$sim.result.name,'.pdf')
pdf(file_name)
my_sim_result$J.long %>% ggplot(aes(timestep, count, color = spp)) +
  geom_line() +
  facet_wrap(~ as.factor(site))
dev.off()

#write out sim result as csv
file_name <- paste0('SIM_OUTPUT/RESULT_JT_',my_sim_result$sim.result.name,'.csv')
write_csv(
  my_sim_result$J.long,
  path = file_name)
