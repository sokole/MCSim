#####################################################
#####################################################
# TEST 1 -- expect 3 taxa

# -- speciation.limit = 0

#####################################################
set.seed(1234) #set random seed

# make a landscape
my.landscape <- make.landscape(
  site.coords = data.frame(
    x = c(1, 2, 3, 4, 5),
    y = c(1, 3, 1, 5, 2)),
  Ef = c(-.8, -.6, 0, .25, .9),
  m = 0.5,
  JM = 10000)

# niche positions, niche breadths, and relative abundances for three species
niche.positions <-  c(-.5, 0, .5)
niche.breadths <- c(.2, .2, 5)
regional.rel.abund <- c(.8, .1, .1)

# run a simulation with 10 generations
sim.result <- metasim(
  landscape = my.landscape,
  trait.Ef = niche.positions,
  trait.Ef.sd = niche.breadths,
  gamma.abund = regional.rel.abund,
  W.r = 0,
  nu = 0.9,
  speciation.limit = 0,
  n.timestep = 10,
  sim.ID = "my_test_sim",
  output.dir.path = "my_sim_output_directory"
)

# plot coenoclines to view niches
plot.coenoclines(sim.result)

# plot dispersal kernal
plot.standardized.disp.kernel(sim.result)

# plot dot plots
plot.dot.plots(sim.result)

sim.result$J.long$spp %>% table()
sim.result$dat.gamma.t0



#####################################################
#####################################################
# TEST 2 -- expect >= 3 taxa

# -- speciation.limit > 0

#####################################################
set.seed(1234) #set random seed

# make a landscape
my.landscape <- make.landscape(
  site.coords = data.frame(
    x = c(1, 2, 3, 4, 5),
    y = c(1, 3, 1, 5, 2)),
  Ef = c(-.8, -.6, 0, .25, .9),
  m = 0.5,
  JM = 10000)

# niche positions, niche breadths, and relative abundances for three species
niche.positions <-  c(-.5, 0, .5)
niche.breadths <- c(.2, .2, 5)
regional.rel.abund <- c(.8, .1, .1)

# run a simulation with 10 generations
sim.result <- metasim(
  landscape = my.landscape,
  trait.Ef = niche.positions,
  trait.Ef.sd = niche.breadths,
  gamma.abund = regional.rel.abund,
  W.r = 0,
  nu = 0.9,
  speciation.limit = 5,
  n.timestep = 10,
  sim.ID = "my_test_sim",
  output.dir.path = "my_sim_output_directory"
)

# plot coenoclines to view niches
plot.coenoclines(sim.result)

# plot dispersal kernal
plot.standardized.disp.kernel(sim.result)

# plot dot plots
plot.dot.plots(sim.result)

sim.result$J.long$spp %>% table()
sim.result$dat.gamma.t0






#####################################################
#####################################################
# TEST 3 -- expect <= 4 taxa

# -- speciation.limit = 0
# -- Jt.0 has 4 taxa
#####################################################
set.seed(1234) #set random seed


# niche positions, niche breadths, and relative abundances for three species
niche.positions <-  c(-.5, -.1, 0, .1)
niche.breadths <- c(.2, .2, .5, 5)
# regional.rel.abund <- c(.8, .1, .1, .01, .01)

my.Jt.0 <- data.frame(
  a = c(5, 5, 5, 8, 8),
  b = c(1, 1, 1, 10, 10),
  c = c(20, 20, 1, 1, 1),
  d = c(1, 1, 20, 20, 20))

# make a landscape
my.landscape <- make.landscape(
  site.coords = data.frame(
    x = c(1, 2, 3, 4, 5),
    y = c(1, 3, 1, 5, 2)),
  Ef = c(-.8, -.6, 0, .25, .9),
  m = 0.5,
  JL = rowSums(my.Jt.0))

# run a simulation with 10 generations
sim.result <- metasim(
  landscape = my.landscape,
  trait.Ef = niche.positions,
  trait.Ef.sd = niche.breadths,
  J.t0 = my.Jt.0,
  W.r = 0,
  speciation.limit = 0,
  n.timestep = 10,
  sim.ID = "my_test_sim",
  output.dir.path = "my_sim_output_directory"
)

# plot coenoclines to view niches
plot.coenoclines(sim.result)

# plot dispersal kernal
plot.standardized.disp.kernel(sim.result)

# plot dot plots
plot.dot.plots(sim.result)

sim.result$J.long$spp %>% table()
sim.result$dat.gamma.t0








