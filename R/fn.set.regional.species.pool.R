#' @title Initialize species pool for fn.metaSIM
#' 
#' @description Initialize species pool for fn.metaSIM
#' 
#' @usage 
#' fn.set.regional.species.pool(alpha.fisher = 1)
#' fn.set.regional.species.pool(gamma.abund = c(.9, .4, .2, .3))

#' @param alpha.fisher Fisher's alpha used to initiate the simulation.
#' @param n.timestep Number of generations (time steps) in the simulation, can affect number of novel species. Default is 0.
#' @param nu Scalar value representing the probability that a novel species 
#' is recruited.  Hubbell's \dQuote{speciation rate}.
#' @param speciation.limit A limit to the number of novel species that can occur 
#' in the simulation, default is 0.  
#' @param JM Total number of individuals in the metacommunity.
#' @param gamma.abund Vector of relative abundances or counts used to set regional species pool RAs
#' @param trait.dispersal Vector of trait scores for dispersal, should be in range 0 to 1, or NULL 
#' @param trait.dispersal.median Median value for species' dispersal trait values.
#' @param trait.dispersal.range Range in variation allowed for species' trait 
#' dispersal values.
#' @param Ef.min Min limit for environmental gradient
#' @param Ef.max Max limit for environmental gradient
#' @param trait.Ef Vector or value of niche positions for species in regional pool
#' @param trait.Ef.sd Vector or value of niche breadths
#' @param taxa.list.prefix A character string to use as a prefix in species' names.
#' @param taxa.list.prefix Text appended as prefix to species names to prevent errors for column names when site by species matrices are created. Default is "spp".
#' 
#' @details Internal function called by \link{fn.metaSIM}.
#' 
#' @export
#' 
fn.set.regional.species.pool <- function ( alpha.fisher = 1,
                                           n.timestep = 0, 
                                           nu = 0.001, 
                                           speciation.limit = 0, 
                                           JM = 1000, 
                                           
                                           gamma.abund = NULL,
                                           
                                           trait.dispersal = NULL, # a vector of dispersal traits
                                           trait.dispersal.median = 1,# or set a median and range of dispersals
                                           trait.dispersal.range = 0, 
                                           
                                           Ef.min = 0, 
                                           Ef.max = 1,
                                           
                                           trait.Ef = NULL,
                                           trait.Ef.sd = NULL,
                                           
                                           taxa.list = NULL,
                                           taxa.list.prefix = "spp"){
  
  # -- set regional RAs if given, otherwise sample dist based on alpha.fisher
  if(is.array(gamma.abund)|is.vector(gamma.abund)){
    alpha.fisher <- NA
  } else {
    # -- default alpha.fisher to 1 when no gamma.abund is provided
    if(is.na(alpha.fisher)) {alpha.fisher <- 1}
    gamma.abund <- fisher.ecosystem(N = JM, nmax = JM, alpha = alpha.fisher)
  }
  
  # -- add columsn for potential new species
  if (is.na(speciation.limit)) {
    n.new.spp <- length(gamma.abund)
  } else {
    n.new.spp <- nu * n.timestep * JM
    if (n.new.spp > speciation.limit) 
      n.new.spp <- speciation.limit
  }
  
  # -- Rel Abundances
  regional.RA <- c(gamma.abund/sum(gamma.abund), rep(0, n.new.spp))
  total.richness <- length(regional.RA)
  
  # -- dispersal traits -- 1 is max dispersal, 0 is min
  while(length(trait.dispersal)<total.richness){
    if(!(is.array(trait.dispersal)|is.vector(trait.dispersal))){
      trait.dispersal <- runif(total.richness, 
                               trait.dispersal.median - trait.dispersal.range, 
                               trait.dispersal.median + trait.dispersal.range)
      trait.dispersal[trait.dispersal > 1] <- 1
      trait.dispersal[trait.dispersal < 0] <- 0
    } else {
      trait.dispersal<-c(trait.dispersal,
                         runif(total.richness-length(trait.dispersal), 
                               trait.dispersal.median - trait.dispersal.range, 
                               trait.dispersal.median + trait.dispersal.range))
      trait.dispersal[trait.dispersal > 1] <- 1
      trait.dispersal[trait.dispersal < 0] <- 0
    }
  }
  
  # -- assign habitat preference traits
  while(length(trait.Ef)<total.richness){
    if(!(is.array(trait.Ef)|is.vector(trait.Ef))){
      trait.Ef <- runif(total.richness, Ef.min, Ef.max)
    } else {
      trait.Ef <- c(trait.Ef,
                    runif(total.richness-length(trait.Ef), 
                          Ef.min, Ef.max))
    }
  }
  
  # -- assign niche breadths
  while(length(trait.Ef.sd)<total.richness){
    if(!(is.array(trait.Ef.sd)|is.vector(trait.Ef.sd))){ #what to do if trait.Ef.sd array is length 0, make NCM
      
      # -- default is neutral community model
      trait.range <- max(trait.Ef)-min(trait.Ef)
      if(!trait.range>0) trait.range <- 1 #make non-zero trait range
      trait.Ef.sd <- rep(1000*trait.range, total.richness) #make neutral community
      
    } else {
      trait.Ef.sd<-c(trait.Ef.sd,
                     sample(x = c(trait.Ef.sd, trait.Ef.sd),
                            size = total.richness-length(trait.Ef.sd), 
                            replace = TRUE))
    }
  }
  
  # -- taxa.list
  while(length(taxa.list)<total.richness){
    if(!(is.array(taxa.list)|is.vector(taxa.list))){
      taxa.list <- paste(taxa.list.prefix, c(1:total.richness), 
                         sep = "")
    } else {
      taxa.list <- c(taxa.list,
                     paste(taxa.list.prefix, 
                           c((length(taxa.list)+1):total.richness), 
                           sep = ""))
    }
  }
  
  if(length(taxa.list)>total.richness){taxa.list<-taxa.list[1:total.richness]}
  
  # -- return dataframe
  return(data.frame(row.names = taxa.list, 
                    taxa.list = taxa.list, 
                    trait.dispersal = trait.dispersal[1:length(taxa.list)], 
                    trait.Ef = trait.Ef[1:length(taxa.list)], 
                    trait.Ef.sd = trait.Ef.sd[1:length(taxa.list)],
                    regional.RA = c(regional.RA), 
                    stringsAsFactors = FALSE))
}