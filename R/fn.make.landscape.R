#' fn.make.landscape
#' 
#' @title Make a square landscape for fn.metaSIM
#' 
#' @description Function to make a square landscape for \link{fn.metaSIM}.
#' 
#' @usage 
#' fn.make.landscape(n.sites = 100, JL.scale = 1, Ef.scale = -1, 
#' Ef.specificity = 0.2, Ef.specificity.scale = 1, IL.scale = 2, 
#' IL.intensity = 5, JM = 1000, JL.min = 1, ...)
#' 
#' @param n.sites Number of sites in the landscape.
#' @param JL.scale Numeric value describing spatial heterogeneity in \code{JL}.  
#' See \link{fn.metaSIM}.
#' @param Ef.scale Numeric value describing spatial heterogeneity in \code{Ef}.  
#' See \link{fn.metaSIM}.
#' @param Ef.specificity Specificity of environmental filter.  See \link{fn.metaSIM}.
#' @param Ef.specificity.scale Numeric value describing the spatial heterogeneity 
#' in \code{Ef.specificity}.  Currently can only take on values of -1, 1, or NA.  
#' See \link{fn.metaSIM}.  
#' @param IL.scale Numeric value describing spatial heterogeneity in \code{IL}.  
#' See \link{fn.metaSIM}.
#' @param IL.intensity Number of individuals in the immigrant pool.  This value, 
#' along with \code{JL}, determines \code{m} at each site.
#' @param JM Total number of individuals in the metacommunity.
#' @param JL.min Minimum number of individuals an assemblage is allowed to have.
#' @param \dots other values that are passed on to internal functions.
#' 
#' @seealso \link{fn.metaSIM}
#' 
#' @export
#' 
fn.make.landscape<-function(
  n.sites = 100, 
  JL.scale = 1, 
  Ef.scale = -1, 
  Ef.specificity = 0.2,
  Ef.specificity.scale = 1,
  IL.scale = 2, 
  IL.intensity = 5, 
  JM = 1000, 
  JL.min = 1, ...) 
{
  require(vegan)
  dat.landscape <- data.frame(row.names = paste("site", 1:n.sites, 
                                                sep = ""), expand.grid(x = 1:sqrt(n.sites), y = 1:sqrt(n.sites)))
  site.list <- as.character(row.names(dat.landscape))
  geo.dist.landscape <- dist(dat.landscape[, c("x", "y")])
  mod.pcnm.landscape <- pcnm(geo.dist.landscape)
  pcnm.list <- names(as.data.frame(mod.pcnm.landscape$vectors))
  dat.landscape <- cbind(dat.landscape[site.list, ], mod.pcnm.landscape$vectors[site.list, 
                                                                                ])
  if (is.na(JL.scale)) {
    dat.landscape$JL <- 1
    JL.pcnm <- "homogenous"
  }else {
    if (JL.scale < 0) {
      dat.landscape$JL <- sample(x = c(1, 10^JL.scale), 
                                 size = nrow(dat.landscape), replace = TRUE)
      JL.pcnm <- "async"
    }else if (JL.scale >= 1) {
      dat.landscape$JL <- ifelse(dat.landscape$x > mean(dat.landscape$x), 
                                 10^(-1 * JL.scale), 1)
      JL.pcnm <- "sync"
    }else {
      pcnm.var <- trunc((length(pcnm.list) - 1) * (1 - 
                                                     JL.scale), 0) + 1
      JL.temp <- dat.landscape[, pcnm.list[pcnm.var]]
      dat.landscape$JL <- ifelse(JL.temp > 0, JL.temp, 
                                 0)
      JL.pcnm <- pcnm.var
      rm(JL.temp, pcnm.var)
    }
  }
  JL <- dat.landscape$JL
  JL <- round(JL * ((JM - (JL.min * length(JL)))/(sum(JL))), 
              0) + JL.min
  dat.landscape$JL <- JL
  
  
  
  homogenizing.Ef<-FALSE
  if (Ef.specificity < 0) {
    dat.landscape$Ef.specificity <- 0
    homogenizing.Ef <- TRUE
  }else{
    if (is.na(Ef.specificity.scale)){
      dat.landscape$Ef.specificity<-Ef.specificity
    }else if(is.numeric(Ef.specificity.scale)){
      if(Ef.specificity.scale>=1){
        dat.landscape$Ef.specificity<-ifelse(dat.landscape$x > mean(dat.landscape$x),0,Ef.specificity)
        #syncronous Ef.specificty
      }else if(Ef.specificity.scale<0){
        dat.landscape$Ef.specificity<-sample(c(0, Ef.specificity), size = nrow(dat.landscape),replace = TRUE)
        #async Ef.specificity
      }
    }
  }

  
  
  if (is.na(Ef.scale)) {
    dat.landscape$Ef <- 0
    Ef.pcnm <- "homogenous"
  }else if (is.numeric(Ef.scale)) {
    if (Ef.scale >= 1) {
      dat.landscape$Ef <- ifelse(dat.landscape$x > mean(dat.landscape$x),0, 1)
      Ef.pcnm <- "sync"
    }else if (Ef.scale < 0) {
      dat.landscape$Ef <- sample(c(0, 1), size = nrow(dat.landscape),replace = TRUE)
      Ef.pcnm <- "async"
    }else {
      if (homogenizing.Ef) {
        pcnm.var <- trunc((length(pcnm.list) - 1) * (1 - 
                                                       Ef.scale), 0) + 1
        dat.landscape$Ef <- dat.landscape[, pcnm.list[pcnm.var]]
        dat.landscape[dat.landscape$Ef < 0, "Ef"] <- 0
        if (max(dat.landscape$Ef) > 0) 
          dat.landscape$Ef <- dat.landscape$Ef/max(dat.landscape$Ef)
        Ef.pcnm <- pcnm.var
        rm(pcnm.var)
      }else {
        pcnm.var <- trunc((length(pcnm.list) - 1) * (1 - 
                                                       Ef.scale), 0) + 1
        dat.landscape$Ef <- dat.landscape[, pcnm.list[pcnm.var]]
        dat.landscape$Ef <- dat.landscape$Ef - min(dat.landscape$Ef)
        if (max(dat.landscape$Ef) > 0) 
          dat.landscape$Ef <- dat.landscape$Ef/max(dat.landscape$Ef)
        Ef.pcnm <- pcnm.var
        rm(pcnm.var)
      }
    }
  }else {
    dat.landscape$Ef <- ifelse(dat.landscape$x > mean(dat.landscape$x),0, 1)
    Ef.pcnm <- "sync"
  }
  
  
  if (is.na(IL.scale)) {
    dat.landscape$IL <- IL.intensity
    IL.pcnm <- "homogenous"
  }else {
    if (IL.scale < 0) {
      IL.temp <- trunc(IL.intensity * sample(x = c(1, 10^IL.scale), 
                                             size = nrow(dat.landscape), replace = TRUE))
      IL.temp[IL.temp < 1] <- 1
      dat.landscape$IL <- IL.temp
      IL.pcnm <- "async"
      rm(IL.temp)
    }else if (IL.scale >= 1) {
      IL.temp <- ifelse(dat.landscape$x > mean(dat.landscape$x), 
                        IL.intensity * (10^(-1 * IL.scale)), IL.intensity)
      IL.temp[IL.temp < 1] <- 1
      dat.landscape$IL <- IL.temp
      IL.pcnm <- "sync"
      rm(IL.temp)
    }else {
      pcnm.var <- trunc((length(pcnm.list) - 1) * (1 - 
                                                     IL.scale), 0) + 1
      IL.temp <- dat.landscape[, pcnm.list[pcnm.var]]
      IL.temp <- IL.temp - min(IL.temp)
      if (max(IL.temp) != 0) 
        IL.temp <- IL.temp/max(IL.temp)
      IL.temp <- ifelse(IL.temp > 0, IL.temp, 1)
      IL.temp <- round(IL.intensity * IL.temp, 0)
      IL.pcnm <- pcnm.var
      dat.landscape$IL <- IL.temp
      rm(IL.temp, pcnm.var)
    }
  }
  m.temp <- dat.landscape$IL/(dat.landscape$IL + dat.landscape$JL - 
                                1)
  dat.landscape$m <- m.temp
  landscape.list <- list(dat = dat.landscape, Ef.scale = Ef.scale, 
                         Ef.pcnm = Ef.pcnm, Ef.specificity = Ef.specificity, IL.scale = IL.scale, 
                         IL.pcnm = IL.pcnm, IL.intensity = IL.intensity, JL.scale = JL.scale, 
                         JL.pcnm = JL.pcnm, JM = JM, site.list = site.list, pcnm.list = pcnm.list, 
                         geo.dist = geo.dist.landscape)
  return(landscape.list)
}