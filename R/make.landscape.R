#' @title Make a simulation landscape
#' 
#' @aliases fn.make.landscape
#' 
#' @description Define the attributes of a MCSim landscape, including number of sites, area, carrying capacity, and local immigration rates.
#' 
#' @usage make.landscape(JM = 10000, m = 0.1)
#' make.landscape(site.coords = c(1:10), m = 0.1, JM = 10000) 
#' 
#' @param site.coords A data.frame of site coordinates. Can be 1, 2, or more dimensions
#' @param dist.mat Alternative to site.coords. Can be a distance matrix or a network map from the igraph package
#' @param site.info A data frame with site information
#' @param JL Scalar or vector number of individuals at each site, overrides JM
#' @param JM Total number of individuals to include in a MCSim simulation.
#' @param m Immigration rate paramter, from Hubbells neutral model. Overrides I.rate.m2.
#' @param I.rate.m2 Alternative to m, immigration rate in number of individuals / m2 / timestep. Default is 1.
#' @param area.m2 Area of each site in m2. Default is 1.
#' @param Ef.specificity Vector of specificity values for environmental filters at each site. If 0 (default), site habitat value is modeled as a single point along an environmental gradient. If > 0, a site's habitat is modeled as a normal curve around a point on an environmental gradient.
#' @param Ef Vector of habitat scores for each site.
#' @param guess.site.coords Binary. If TRUE, Uses PCoA to extract site coordinates if given a distance matrix or network map. Useful to make a map to display sites. Not necessary if igraph input is used because igraph has a function to plot a network map. Default is FALSE.
#' @param list.of.stuff A list that can be used to store other landscape attributes in the landscape object. Useful for storing igraph properties when igraph is used. 
#' 
#' @details There are two steps to creating a metacommunity simulation in MCSim:
#' 1. Make a "landscape" -- The landscape is the “game board” on which the simulation plays out, and it is created using the fn.make.landscape function.
#' 2. Run the simulation -- Once the landscape is created, you can pass the landscape object to fn.metaSIM along with parameter settings that define the rules for how metacommunity dynamics will play out in the metacommunity simulation. Note that the current version of MCSim is zero sum, which means there will always be JM individuals in the simulation during each generation.
#' For a tutorial, see \url{http://rpubs.com/sokole/159425}
#' 
#' @export
#' 

make.landscape<-function(
  # -------------------------------
  # -------------------------------
  # -------------------------------
  # -- data frame inputs
  site.info = data.frame(),
  
  # -- need one or the other, or they have to match. Priority given to dist if they don't
  site.coords = NULL,
  dist.mat = data.frame(),
  
  # -- will fill in info if none given
  
  # -- default metacommunity parameters if none given
  site.ID = NULL,
  JL = NULL,
  JM = NULL,
  m = NULL,
  I.rate.m2 = NULL,
  area.m2 = 1,
  Ef.specificity = 0, # 0 is point specificity along env. gradient
  Ef = 0,
  guess.site.coords = FALSE,
  list.of.stuff = NA,
  ...
){
  # -------------------------------
  # -------------------------------
  # -------------------------------
  
  # -- ESCAPE VAR
  get.the.f.out <- FALSE
  
  # -- extracting variables from the site.info data frame if they are provided in there by the user
  vector_argument_list = c('site.ID','JL','JM','m','I.rate.m2','area.m2','Ef.specificity','Ef')
  
  # extract site.coords from site.info if available
  if(is.null(site.coords) & nrow(dist.mat) == 0){
    try({site.coords <- site.info[ , !names(site.info) %in% vector_argument_list]})
  }
  
  if(is.null(JL) & 'JL' %in% names(site.info)) JL <- site.info$JL
  if(is.null(JM) & 'JM' %in% names(site.info)) JM <- site.info$JM
  if(is.null(m) & 'm' %in% names(site.info)) m <- site.info$m
  if(is.null(I.rate.m2) & 'I.rate.m2' %in% names(site.info)) I.rate.m2 <- site.info$I.rate.m2
  
  # override defaults if these vars are included in site.info
  if('area.m2' %in% names(site.info)) area.m2 <- site.info$area.m2
  if('Ef.specificity' %in% names(site.info)) Ef.specificity <- site.info$Ef.specificity
  if('Ef' %in% names(site.info)) Ef <- site.info$Ef
  
  try({site.coords<-as.data.frame(site.coords)},
      silent = TRUE)
  
  # -------------------------------
  # -- Make geo landscape data.frame
  # -------------------------------
  dist.mat<-as.data.frame(as.matrix(dist.mat)) #make into a data.frame
  
  # --  check info
  if(nrow(dist.mat)>0){
    dat.geo.dist.out<-dist.mat
    if(nrow(site.coords)!=nrow(dist.mat) & guess.site.coords){
      dat.geo.out<-data.frame(cmdscale(dat.geo.dist.out))
      print('I assigned dat.geo for you')
    }else if(nrow(site.coords)!=nrow(dist.mat)){
      dat.geo.out<-data.frame(site.label=c(1:nrow(dat.geo.dist.out)))
      print('Rock on')
    }else if(nrow(site.coords)==nrow(dist.mat)){
      dat.geo.out<-site.coords
    }
  }else if(nrow(site.coords)>0){
    dat.geo.out<-site.coords
    dat.geo.dist.out<-data.frame(as.matrix(dist(site.coords)))
    print('success!!')
  }else if(nrow(site.info)>0){
    get.the.f.out<-TRUE
    print('no geo info!')
  }else{
    get.the.f.out<-TRUE
    print('no landscape info!')
  }
  
  n.obs<-nrow(dat.geo.dist.out)
  
  if(!get.the.f.out){
    # -- turn area into a vector if it is not a vector, otherwise, it remains the same, get's fed up if it's the wrong length
    area.m2<-data.frame(
      dummy=c(1:n.obs),
      area.m2=area.m2)$area.m2
    
    # -------------------------------
    # -- calculate assemblage sizes at sites, JL 
    # -------------------------------
    if( is.null(JL) | !(length(JL)%in%c(1, n.obs)) ){
      JL.wts <- area.m2 / sum(area.m2)
      JL.wts <- JL.wts/sum(JL.wts)
      JL <- round(JL.wts * JM,0)
    }else{
      JL <- data.frame(
        dummy=c(1:n.obs),
        JL=JL)$JL
      JM <- sum(JL)
    }
    
    # -------------------------------
    # -- calculate immigration at sites, IL 
    # -------------------------------
    
    if(!is.null(m) & length(m)%in%c(1,n.obs)){
      m.site <- data.frame(
        dummy=c(1:n.obs),
        m=m)$m
      I.site <- m.site * (JL - 1) / (1 - m.site)
      I.site <- ifelse(is.infinite(I.site), 10*JL, I.site)
    }else if(!is.null(I.rate.m2)){
      I.site <- I.rate.m2 * area.m2
      I.site <- round(I.site,0)
      m.site <- I.site/(I.site + JL - 1)
    }else{
      print('aaaaggghhh, need immigration rate')
    }
    
    
    # -------------------------------
    # -- dat with info
    # -------------------------------
    if(is.null(site.ID)){
      if('site.ID' %in% names(site.info)){
        site.ID <- site.info$site.ID
      }else{
        site.ID <- c(1:nrow(dat.geo.out))
      }
    }

    
    dat.info.default<-data.frame(
      site.ID = site.ID,
      area.m2 = area.m2,
      JL = JL,
      Ef = Ef,
      Ef.specificity = Ef.specificity,
      IL = I.site,
      m = m.site
    )
    
    if(nrow(site.info)>0){
      dat.info.out<-site.info  
      # -- check to see if specific vars need to be filled in
      for(i.var in c('site.ID','area.m2','JL','Ef','Ef.specificity','IL','m')){
        if(!i.var%in%names(dat.info.out)) dat.info.out[,i.var] <- dat.info.default[,i.var]
      }
    }else{
      dat.info.out<-dat.info.default
    }
    
    return(
      list(site.info=dat.info.out,
           site.coords=dat.geo.out,
           dist.mat=dat.geo.dist.out,
           list.of.stuff=list.of.stuff)
    )
  }else{
    print('no landscape for you!')
  }
}

# define alias functions
#' @export
fn.make.landscape <- make.landscape
  