#' fn.plot.sim
#' 
#' @title Plots to visualize a simulation
#' 
#' @description NOT CURRENTLY USED... A function to plot select parameters 
#' and diversity outcomes for the first and last timesteps in a simulation.
#' 
#' @usage 
#' fn.plot.sim(sim.result, dist.type = "horn", 
#' plot.timestep = c(1, sim.result$n.timestep), ...)
#' 
#' @param sim.result}{A result from \link{fn.metaSIM}.
#' @param dist.type}{Distance metric, from \link[vegan]{vegdist} used to 
#' plot an ordination of community composition.
#' @param plot.timestep}{Time steps to include in the plot.  Default is to 
#' plot the initial and final time steps. Must be a vector of no more than 
#' 2 integers.
#' @param \dots}{Other plotting parameters passed onto \link{plot} and \link[vegan]{ordiplot}.
#' 
fn.plot.sim <- function(
  sim.result,
  dist.type="horn",
  plot.timestep=c(1,sim.result$n.timestep),
  ...
){
  # -----------------------------------------------------------------
  # -- Data for plotting results
  # -----------------------------------------------------------------
  
  # -- select as percentiles for timesteps to send to ordination plotting funciton
  # -- does not make plots in this function
  timestep.labels<-paste("t",plot.timestep,sep="")
  
  J.plotting<-data.frame(
    groupID=rep(timestep.labels,each=length(sim.result$x)),
    do.call(rbind,sim.result$J[plot.timestep]),
    stringsAsFactors =FALSE)
  
  landscape<-sim.result$landscape
  taxa.list<-sim.result$taxa.list
  
  # -----------------------------------------------------------------
  # -- Plot landscape
  # -----------------------------------------------------------------
  par(mfrow=c(2,3),mex=.5,cex=.7)
  with(landscape$dat,{
    response.var<-JL
    ordisurf(cbind(x,y), response.var,
             pch=21,
             cex=.5+1.5*abs(response.var)/max(abs(response.var)),
             bg=ifelse(response.var<0,2,NA),
             main="JL")
    response.var<-Ef
    ordisurf(cbind(x,y), response.var,
             pch=21,
             cex=.5+2*response.var,
             main="Ef")
    response.var<-m
    ordisurf(cbind(x,y), response.var,
             pch=21,
             cex=.5+1.5*abs(response.var)/max(abs(response.var)),
             bg=ifelse(response.var<0,2,NA),
             main="m")
  })
  
  # -----------------------------------------------------------------
  # -- Plot Metacommunity
  # -----------------------------------------------------------------
  
  # -- rank abundance curves
  dat.rankabund<-do.call(rbind,
                         by(data=ifelse(J.plotting[,taxa.list]>0,1,0),
                            INDICES=J.plotting$groupID,FUN=colSums))
  
  show.taxa<-colSums(dat.rankabund)>0
  
  par(las=2)
  barplot(dat.rankabund[,show.taxa],beside=TRUE,col=c(1,2),main="Rank Occurrence")
  
  # -- trait diversity
  dat.traits<-data.frame(
    row.names=sim.result$taxa.list,
    taxa.list=taxa.list,
    trait.dispersal=sim.result$trait.dispersal,
    trait.Ef=sim.result$trait.Ef
  )[show.taxa,]
  
  # -- plot trait vals
  barplot(t(as.matrix(dat.traits[,c("trait.dispersal","trait.Ef")])),
          beside=TRUE,col=c(1,NA),
          main="Trait Scores")
  legend(x="bottomright",
         legend=c("Dispersal","Ef affinity"),
         fill=c(1,NA),
         bty="o",cex=.6)
  
  # -- make PCOA model
  mod.pcoa<-cmdscale(vegdist(decostand(J.plotting[,taxa.list],"log"),dist.type))
  
  # -- plot PCOA for t0 and selected quantile
  mod.plot<-ordiplot(mod.pcoa,type="n")
  
  #   points(mod.plot,what="sites",select=J.plotting$groupID==timestep.labels[1],col=1,pch=1)
  #   points(mod.plot,what="sites",select=J.plotting$groupID==timestep.labels[2],col=2,pch=2)
  #   
  color.vec<-ifelse(J.plotting$groupID==timestep.labels[1],1,2)
  pch.vec<-ifelse(J.plotting$groupID==timestep.labels[1],1,2)
  cex.vec<-1.5*(landscape$dat$JL/max(landscape$dat$JL))+.5
  
  points(mod.plot,what="sites",col=color.vec,pch=pch.vec,cex=cex.vec)
  
  
  ordihull(mod.pcoa,
           groups=J.plotting$groupID,
           show.groups=timestep.labels[1],
           display="sites",col=1,lty=1)
  ordihull(mod.pcoa,
           groups=J.plotting$groupID,
           show.groups=timestep.labels[2],
           display="sites",col=2,lty=2)
  
  legend(x="bottomleft",
         legend=timestep.labels,
         pch=c(1:2),col=c(1:2),lty=c(1:2),
         bty="n",cex=.75)
  
  dat.plot<-data.frame(
    J.plotting,
    m=sim.result$m,
    Ef=sim.result$Ef)
  
  mod.envfit<-envfit(mod.pcoa,
                     env=data.frame(
                       m=dat.plot$m,
                       Ef=dat.plot$Ef),
                     permutations=0,
                     na.rm=TRUE)
  plot(mod.envfit)
}
