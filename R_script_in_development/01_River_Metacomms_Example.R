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


# -----------------------------------------------------
# -- make some dot blot plots
# ----------------------------

# -- extract timesteps to plot
J.long <- filter(my_sim_result$J.long, timestep%in%c(1,2,max(timestep)))

# -- group spp counts by time and site to calculate site count totals
J.time.site <- group_by(J.long, timestep, site)
JLs <- summarise(J.time.site, site.totals=sum(count))
J.JLs <- full_join(J.long, JLs, by=c('timestep','site'))

# -- calculate relative abundances (RAs) from counts and site totals, remove
# observations with RAs of 0
J.RAs.long<-J.JLs
J.RAs.long$RA<-J.RAs.long$count/J.RAs.long$site.totals
J.RAs.long<-filter(J.RAs.long, RA > 0)

# -- add environmental data and species names for plotting
J.RAs.long<-mutate(J.RAs.long, 
                   spp.no = as.numeric(as.factor(spp)))

# extract site info from the landscape object
# rank sites locations along the env gradient
site_info <- my_landscape$site.info %>%
  mutate(Ef_rank = rank(Ef))

# extract starting regional species pool from sim object
# rank spp positions along env gradient
spp_info <- my_sim_result$dat.gamma.t0 %>%
  mutate(trait_rank = rank(trait.Ef))

# J.RAs.long$region <- d.site.1D[J.RAs.long$site, 'region']
J.RAs.long$Ef_rank <- site_info[J.RAs.long$site, 'Ef_rank']
J.RAs.long$trait_rank <- spp_info[J.RAs.long$spp,'trait_rank']

# -- make a species characteristic data frame for labeling the spp axis
d.spp.RAs <- as.data.frame(J.RAs.long %>% 
                             group_by(spp) %>% 
                             summarise(max.RA = max(RA),
                                       spp.no = spp.no[1],
                                       trait_rank = trait_rank[1]))

# -- only label species with RAs over 0.40
spp.labels <- dplyr::filter(d.spp.RAs, max.RA > .4)

# -- make plot
p <- ggplot(J.RAs.long, aes(trait_rank, 
                            Ef_rank, 
                            size = RA))

p + geom_point() + 
  facet_grid(. ~ timestep) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 6)) +
  scale_size('Relative\nAbundance', range = c(.5,4)) +
  labs(title = 'Simulated Metacommunity Simulation\nNo. generations ( ---> )',
       x = 'Site ID',
       y = 'Species ID') +
  scale_y_continuous(
    breaks = spp.labels$env.rank,
    labels = spp.labels$spp) +
  theme_bw()


