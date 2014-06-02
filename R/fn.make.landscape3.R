#' fn.make.landscape3
#' 
#' @title Make a customizable landscape for fn.metaSIM3
#' 
#' @description Function to make a customizable landscape for \link{fn.metaSIM3}.
#' 
#' @param x vector of x coordinates for sites/ponds
#' @param y vector of y coordinates for sites/ponds
#' @param management vector indicating site/pond management, or another grouping variable
#' @param area_m2 vector of pond surface areas (m2) 
#' @param Ef vector of values for site/pond environmental filter
#' @param Ef.specificity.type NOT USED, set to 0 as a default, environmental filter 
#' type/strength, positive numeric value, or 'H1', 'H2', or 'neutral'
#' @param JM total metacommunity size
#' @param I.rate.m2 areal immigration rate
#' @param mgmt.JL.increase NOT USED, default is 1 (no effect), factor indicating the 
#' percent increase in population size at managed sites.
#' 
#' @seealso \link{fn.metaSIM3}
#' 
#' @export
#' 
fn.make.landscape3<-function(
  # -- vectors
  x = d.geo$x,
  y = d.geo$y,
  management = d.siteinfo$managed,
  area_m2 = d.env.raw$area_m,
  Ef = Ef.obs,
  Ef.specificity.type = 0, # positive numeric value, or 'H1', 'H2', or 'neutral'
  # -- scalar values
  JM = JM.est,
  I.rate.m2 = 1,
  mgmt.JL.increase = .25){
  
  require(vegan)

  # -------------------------------
  # -- calculate assemblage sizes at sites, JL 
  # -- JL influenced by management
  # -------------------------------
  JL.wts.E <- area_m2 / sum(area_m2)
  JL.wts <- (JL.wts.E * (1 + mgmt.JL.increase*management))
  JL.wts<-JL.wts/sum(JL.wts)
  JL <- round(JL.wts * JM,0)
  
  # -------------------------------
  # -- calculate immigration at sites, IL 
  # -------------------------------
  I.site <- I.rate.m2 * area_m2
  I.site <- round(I.site,0)
  m.site <- I.site/(I.site + JL - 1)
  
  # -------------------------------
  # -- assign Ef.specificity at sites 
  # -------------------------------
  if(is.na(Ef.specificity.type)){
    Ef.specificity<-rep(0, length(x))
  }else if(is.numeric(Ef.specificity.type)){
    Ef.specificity<-rep(Ef.specificity.type, length(x))
  }else if(length(grep('H1',Ef.specificity.type))>0){
    Ef.specificity<-rep(0, length(x))
  }else if(length(grep('H2',Ef.specificity.type))>0){
    Ef.specificity<-management
  }else if(length(agrep('neutral',Ef.specificity.type))>0){
    Ef.specificity<-rep(2, length(x))
  }else{
    Ef.specificity<-rep(0, length(x))
  }
  
  # -------------------------------
  # -- Make landscape data.frame
  # -------------------------------
  dat.landscape<-data.frame(
    x = x,
    y = y,
    management = management,
    area_m2 = area_m2,
    JL = JL,
    Ef = Ef,
    Ef.specificity = Ef.specificity,
    IL = I.site,
    m = m.site
  )
  
  # -------------------------------
  # -- make spatial variables using pcnm
  # -------------------------------
  geo.dist.landscape <- dist(dat.landscape[,c('x','y')])
  mod.pcnm<-pcnm(geo.dist.landscape)
  dat.pcnm<-data.frame(scores(mod.pcnm))
  dat.landscape <- data.frame(
    dat.landscape,
    dat.pcnm)
  
  
  # -------------------------------
  # -- return data as landscape.list
  # -------------------------------
  landscape.list <- list(dat = dat.landscape, 
                         Ef.scale = NA, 
                         Ef.pcnm = NA, 
                         Ef.specificity = Ef.specificity.type, 
                         IL.scale = NA, 
                         IL.pcnm = NA, 
                         IL.intensity = I.rate.m2, 
                         JL.scale = NA, 
                         JL.pcnm = NA, 
                         JM = JM, 
                         site.list = row.names(dat.landscape), 
                         pcnm.list = names(dat.pcnm), 
                         geo.dist = geo.dist.landscape)
  
  return(landscape.list)
}