MCSim
=====
#### An R package for metacommunity simulations and functions to calculate biodiversity metrics  
MCSim is a package for R that can be used to design lottery-based simulations of metacommunities, which can be used to test assumptions about how metacommunity dynamics influence emergent biodiversity patterns. See my paper in Ecological Modelling ([link](http://www.sciencedirect.com/science/article/pii/S0304380014004918)), which uses MCSim [v0.3](https://github.com/sokole/MCSim/releases/tag/v0.3). In [v0.4.2](https://github.com/sokole/MCSim/releases/tag/v0.4.2), users can make use of igraph to design network topology and seed simulations with observed data. Most recently, [v0.4.9](https://github.com/sokole/MCSim/tree/v0.4.9) provides users control over temporal variability in the landscape as well as some new plotting functions. [v0.5.0](https://github.com/sokole/MCSim/releases/tag/v0.5.0) provides plotting functions and an intro vignette ("MCSim-intro"). Please find the current version under development at https://github.com/sokole/MCSim.   

#### Tutorials

* [An introduction to metacommunity simulations (MCSim) for R](https://rpubs.com/sokole/MCSim-intro)
* [Modeling disturbance with metacommunity simulations (MCSim)](https://rpubs.com/sokole/MCSim-disturbance)
* [Exploring species area relationships with MetaCommunity Simulations (MCSim)](https://rpubs.com/sokole/MCSim-SARs) - still in development

#### Useful links and related projects

* [MCSim details, news, and tutorials](https://sites.google.com/site/metacommunitysimulation/)
* [Application of metacommunity simulations at the McMurdo Dry Valleys LTER](http://mcm.lternet.edu/content/metacommunity-dynamics-simulations-diatoms-antarctic-ponds)
* [The LTER Metacommunities Working Group project page](https://sites.google.com/site/ltermetacommunities/home)
* [LTER Metacommunities on github](https://github.com/sokole/ltermetacommunities/)

Please feel free to email me with questions, suggestions, etc.  
#### Eric R. Sokol  
sokole@gmail.com

#### How to install MCSim
```
# Install without vignette(s)
remotes::install_github("sokole/MCSim")

# Install with vignette(s)
remotes::install_github("sokole/MCSim", build_vignettes = TRUE)
```
