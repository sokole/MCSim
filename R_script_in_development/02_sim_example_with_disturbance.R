set.seed(1234) #set random seed

# make a landscapes
my.landscape.1 <- make.landscape(
  site.coords = data.frame(
    x = c(1, 2, 3, 4, 5),
    y = c(1, 3, 1, 5, 2)),
  Ef = c(-.8, -.6, 0, .25, .9),
  m = 0.5,
  JM = 10000)

my.landscape.2 <- make.landscape(
  site.coords = data.frame(
    x = c(1, 2, 3, 4, 5),
    y = c(1, 3, 1, 5, 2)),
  Ef = c(-.8, -.6, 0, 0, 0),
  m = 0.5,
  JM = 10000)

# niche positions, niche breadths, and relative abundances for three species
niche.positions <-  c(-.5, 0, .5)
niche.breadths <- c(.2, .2, 5)
regional.rel.abund <- c(.8, .1, .1)

sim.result <- metasim.disturbance(
  scenario_name = "niche_shift",
  landscape_list = list(my.landscape.1,my.landscape.2),
  time_interval_durations = c(10,10),
  trait.Ef = niche.positions,
  trait.Ef.sd = niche.breadths,
  gamma.abund = regional.rel.abund,
  W.r = 0,
  nu = 0.001,
  output.dir.path = "my_disturbance_sim_output_directory")


plot_coenoclines(sim.result)
plot_standardized_disp_kernel(sim.result)
plot_dot_plots(sim.result)
