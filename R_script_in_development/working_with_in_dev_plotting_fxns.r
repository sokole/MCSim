library(MCSim)

# load some sim output
load('R_script_in_development/SIM_METACOMM_nicheScaling-1_m-0.001_nu-0.001_w-1e+05_rep-168_2019-08-29_181004_20190829_181231.rda')

plot_coenoclines(sim_result = sim.result,
                 plot_trait_Ef_sd_0s = TRUE)
#############################
######
my_niche_breadth_coef <- 1
my_timestep <- 1
# y_max_val <- 7
#########
# these should be the same for all simulations
dat.gamma.t0 <- sim.result$dat.gamma.t0
dat.comm.long <- sim.result$J.long %>% filter(timestep == my_timestep, count > 0)
dat.comm.wide <- dat.comm.t0.long %>% spread(spp, count, fill = 0)

dat.gamma_no_zeros <- dat.gamma.t0 %>% filter(taxa.list %in% dat.comm.long$spp)

# -------------------------------------------------------------
# -- plot the coenoclines, species' niche positions along the environmental gradient
# --
# -- Note that the "rug" marks along the x-axis denote the sites' locations
# -- on the environmental gradient
# ------------------------

env.axis <- sim.result$landscape$site.info$Ef

my_mu <- dat.gamma_no_zeros$trait.Ef
my_sigma <- dat.gamma_no_zeros$trait.Ef.sd

my_sigma_adj <- (my_sigma * my_niche_breadth_coef)

graphics.off()
windows(6,4)
plot_coenoclines(
  Ef = env.axis,
  trait_Ef = my_mu,
  trait_Ef_sd = my_sigma_adj,
  plot_trait_Ef_sd_0s = FALSE
)

graphics.off()
windows(6,4)
plot_coenoclines(
  Ef = env.axis,
  trait_Ef = my_mu,
  trait_Ef_sd = my_sigma_adj,
  plot_trait_Ef_sd_0s = FALSE,
  y_max_val = 10
)

graphics.off()
windows(6,4)
plot_coenoclines(
  Ef = env.axis,
  trait_Ef = my_mu,
  trait_Ef_sd = my_sigma_adj,
  plot_trait_Ef_sd_0s = TRUE,
  y_max_val = 10
)

graphics.off()
windows(6,4)
plot_coenoclines(
  Ef = env.axis,
  trait_Ef = my_mu,
  trait_Ef_sd = my_sigma_adj)
###############################################################################
###############################################################################
###############################################################################
###############################################################################

# plotting disp kernels

graphics.off()
windows(6,4)

plot_standardized_disp_kernel(
  w = c(0, 1, 100),
  distance_matrix = sim.result$landscape$dist.mat
)

plot_standardized_disp_kernel(
  w = c(0, 1, 100)
)

plot_standardized_disp_kernel(
  w = 100
)

plot_standardized_disp_kernel(
  sim_result = sim.result
)
