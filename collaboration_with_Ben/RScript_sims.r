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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# landscape design common to all sims -- two sets of clumped sites
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# generating pre and post disturbance landscapes
xy_coordinates <- data.frame(
  x = c(1, 1.5, 3, 4, 5),
  y = c(1, 1.5, 1, 6, 2))

plot(xy_coordinates)

