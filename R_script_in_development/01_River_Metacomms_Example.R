# library(MCSim)
# library(tidyverse)

# read in source data, calc regional relative abundances
dat_comm_wide <- read_csv("years1419.csv")[,-1]

dat_comm_long <- dat_comm_wide %>%
  tidyr::pivot_longer(
    cols = -c(location_id, year))

dat_comm_summary <- dat_comm_long %>%
  group_by(name) %>%
  summarize(mean_count = mean(value)) %>%
  dplyr::filter(mean_count > 0)


dat_comm_gamma_RA <- dat_comm_summary %>%
  mutate(RA = mean_count/sum(mean_count))

spp_list <- dat_comm_gamma_RA$name %>% unique()

obs_JL <- dat_comm_wide[,spp_list] %>% rowSums()

# Make MCSim landscape
n_sites <- 10

my_landscape <- MCSim::make.landscape(
  site.coords = data.frame(
    x = runif(n_sites),
    y = runif(n_sites)),
  m = 10^runif(n_sites,-4,0),
  Ef = rnorm(n_sites),
  JL = sample(obs_JL, size = n_sites, replace = TRUE))

# niche positions, niche breadths, and relative abundances for three species
niche_positions <- rnorm(length(spp_list)) # norm distribution
niche_breadths <- exp(rnorm(length(spp_list))) # lognorm dist
regional_RA <- dat_comm_gamma_RA$RA # mean RAs from data

# run a simulation with 10 generations
my_sim_result <- MCSim::metasim(
  landscape = my_landscape,
  trait.Ef = niche_positions,
  trait.Ef.sd = niche_breadths, 
  gamma.abund = regional_RA,
  W.r = 0,
  nu = 0.001,
  n.timestep = 10, 
  sim.ID = "my_sim",
  output.dir.path = "my_sim_output_directory")


# # need to fix plot_ceonoclines
# plot_coenoclines(my_sim_results)

# plot dispersal kernel
MCSim::plot_standardized_disp_kernel(sim_result = my_sim_result)

# explore other potential dispersal kernels on the landscape
MCSim::plot_standardized_disp_kernel(landscape = my_sim_result$landscape,
                                     w = c(0,1,10,100,1000))


# # plot whole sim -- takes a min or more
# MCSim::plot_coenoclines(sim_result = my_sim_result)

# plot subset
MCSim::plot_coenoclines(
  landscape = my_sim_result$landscape,
  trait_Ef = my_sim_result$dat.gamma.t0$trait.Ef[1:20],
  trait_Ef_sd = my_sim_result$dat.gamma.t0$trait.Ef.sd[1:20])

# dot plot of sim results
MCSim::plot_dot_plots(my_sim_result)

# choolse timesteps and RA threshold for taxa label
MCSim::plot_dot_plots(
  my_sim_result, 
  timesteps = c(1:3,8,10),
  spp_label_threshold_RA = 0.30)


