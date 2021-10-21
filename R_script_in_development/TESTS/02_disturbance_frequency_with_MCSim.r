# Example 2
# Modeling the effect of disturbance frequency

set.seed(1234) #set random seed



# niche positions, niche breadths, and relative abundances for three species
niche.positions <-  c(-.5, 0, .5)
niche.breadths <- c(.3, .3, .3)

d.comm.t0 <- data.frame(
  spp_1 = c(10000,10000,10000,0, 0),
  spp_2 = c(0, 0, 10000, 10000, 10000),
  spp_3 = c(0, 0, 20000, 0, 0))


# make a landscapes -- initial landscape
my.landscape.1 <- make.landscape(
  site.coords = data.frame(
    x = c(1, 2, 3, 4, 5),
    y = c(1, 3, 1, 5, 2)),
  Ef = c(-.8, -.5, 0, .25, .5),
  m = 0.5,
  JL = rowSums(d.comm.t0))

# disturbed landscape
my.landscape.2 <- make.landscape(
  site.coords = data.frame(
    x = c(1, 2, 3, 4, 5),
    y = c(1, 3, 1, 5, 2)),
  Ef = rep(.9, 5), # Note change in Ef
  m = 0.5,
  JL = rowSums(d.comm.t0))


# Example 2A -- high frequency disturbance

# How many times should the disturbance occur?
n_cycles <- 10

# make a list of landscapes to cycle disturbance and recovery 
my_landscape_list <- c(
  list(my.landscape.1 = my.landscape.1), # initial landscape
  rep(
    list(my.landscape.2 = my.landscape.2, # disturbed landscape
         my.landscape.1 = my.landscape.1), # back to initial landscape
    n_cycles)) # how many times to switch back and forth?

# a vector for time intervals for each landscape type
recovery_duration <- 15 # n timesteps
disturbance_duration <- 3 # n timesteps -- currently min allowable duration is 3 timesteps
my_time_interval_durations <- c(recovery_duration, rep(c(disturbance_duration, recovery_duration), n_cycles))

dist_frequency = n_cycles / sum(my_time_interval_durations)

# run the simulation
sim.result <- metasim.disturbance(
  scenario.name = "disturbanc_high_freq",
  landscape.list = my_landscape_list,
  time.interval.durations = my_time_interval_durations,
  trait.Ef = niche.positions,
  trait.Ef.sd = niche.breadths,
  J.t0.initial = d.comm.t0,
  W.r = 20,
  nu = 0.00000001,
  output.dir.path = "my_disturbance_sim_output_directory")

# look at dispersal kernel
plot.standardized.disp.kernel(w = 20, landscape = my.landscape.1)

# view map of sites
my.landscape.1$site.coords %>%
  ggplot(aes(x,y,size = my.landscape.1$site.info$JL)) +
  geom_point()

# Landscape 1 coenoclines
plot.coenoclines(landscape = my.landscape.1,
                 trait.Ef = niche.positions,
                 trait.Ef.sd = niche.breadths)

# Landscape 2 coenoclines -- note the shift of site positions on the x-axis
plot.coenoclines(landscape = my.landscape.2,
                 trait.Ef = niche.positions,
                 trait.Ef.sd = niche.breadths)

# look at change in composition over time
plot.dot.plots(sim.result,timesteps = c(1,2,15:22,195))





# 2B -- Lower frequency disturbance

# How many times should the disturbance occur?
n_cycles <- 1

# make a list of landscapes to cycle disturbance and recovery 
my_landscape_list <- c(
  list(my.landscape.1 = my.landscape.1), # initial landscape
  rep(
    list(my.landscape.2 = my.landscape.2, # disturbed landscape
         my.landscape.1 = my.landscape.1), # back to initial landscape
    n_cycles)) # how many times to switch back and forth?

# a vector for time intervals for each landscape type
recovery_duration <- 95 # n timesteps
disturbance_duration <- 3 # n timesteps -- currently min allowable duration is 3 timesteps
my_time_interval_durations <- c(recovery_duration, rep(c(disturbance_duration, recovery_duration), n_cycles))

dist_frequency = n_cycles / sum(my_time_interval_durations)

# run the simulation
sim.result <- metasim.disturbance(
  scenario.name = "disturbanc_low_freq",
  landscape.list = my_landscape_list,
  time.interval.durations = my_time_interval_durations,
  trait.Ef = niche.positions,
  trait.Ef.sd = niche.breadths,
  J.t0.initial = d.comm.t0,
  W.r = 20,
  output.dir.path = "my_disturbance_sim_output_directory")

# view map of sites
my.landscape.1$site.coords %>%
  ggplot(aes(x,y,size = my.landscape.1$site.info$JL)) +
  geom_point()

# Landscape 1 coenoclines
plot.coenoclines(landscape = my.landscape.1,
                 trait.Ef = niche.positions,
                 trait.Ef.sd = niche.breadths)

# Landscape 2 coenoclines -- note the shift of site positions on the x-axis
plot.coenoclines(landscape = my.landscape.2,
                 trait.Ef = niche.positions,
                 trait.Ef.sd = niche.breadths)

# look at change in composition over time
plot.dot.plots(sim.result)





# Example 2C -- No disturbance, reference sim:
sim.result <- metasim(
  scenario.name = "no_disturbance",
  landscape = my.landscape.1,
  n.timestep = 195,
  trait.Ef = niche.positions,
  trait.Ef.sd = niche.breadths,
  J.t0.initial = d.comm.t0,
  W.r = 0,
  nu = 0,
  output.dir.path = "my_disturbance_sim_output_directory")

# view map of sites
my.landscape.1$site.coords %>%
  ggplot(aes(x,y,size = my.landscape.1$site.info$JL)) +
  geom_point()

# Landscape 1 coenoclines
plot.coenoclines(landscape = my.landscape.1,
                 trait.Ef = niche.positions,
                 trait.Ef.sd = niche.breadths)

# Landscape 2 coenoclines -- note the shift of site positions on the x-axis
plot.coenoclines(landscape = my.landscape.2,
                 trait.Ef = niche.positions,
                 trait.Ef.sd = niche.breadths)

# look at change in composition over time
plot.dot.plots(sim.result)

