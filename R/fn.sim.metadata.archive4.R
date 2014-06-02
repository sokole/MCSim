#' fn.sim.metadata.archive4
#' 
#' @title Archive your simulation
#' 
#' @description Saves metadata to a .csv file in an output directory, default 
#' is SIM_OUTPUT. Optionally saves the simulation output as a .rda file.  
#' These files can be big.
#' 
#' @usage 
#' fn.sim.metadata.archive4(sim.result = sim.result, var.dir = "SIM_OUTPUT",
#' var.filetype = ".rda", save.sim = TRUE, q.order = q.order, 
#' keep.timesteps = c(1, 10))
#' 
#' @param sim.result An object produced by \link{fn.metaSIM}, metacommunity simulation 
#' output.
#' @param var.dir Directory where simulation outputs will be saved.
#' @param var.filetype File type used to save simulation outputs.  Currently must 
#' be .rda.
#' @param save.sim Boolean (\code{TRUE} or \code{FALSE}), should the simulation output 
#' be saved as a .rda file?  If \code{FALSE}, only the metadata file from the simulation 
#' will be saved.
#' @param q.order Order \code{q} to use to calculate diversity statistics.  Can be a 
#' single integer or a vector of integers.
#' @param keep.timesteps Timesteps to report in the metadata file.  
#' 
#' @export
#'                    
fn.sim.metadata.archive4 <- function(
  sim.result=sim.result,
  var.dir="SIM_OUTPUT",
  var.filetype=".rda",
  save.sim=TRUE,
  q.order=q.order,
  keep.timesteps=c(1,10)){
  
  require(plyr)
  calc.comm.stats<-!is.na(keep.timesteps)

  # -- check for output directory, create if it does not exit
  if(!file.exists(var.dir)) dir.create(var.dir)
  
  # -- create sim.result name using parameter settings and timestamp
  sim.result.name<-sim.result$sim.result.name
  sim.result.filename<-paste(var.dir,"/",sim.result.name,var.filetype,sep="")
  
  # -- extract metadata from sim.result file
  n.timestep<-sim.result$n.timestep
  
  dat.trait<-data.frame(
    taxa.list=sim.result$taxa.list,
    trait.dispersal=sim.result$trait.dispersal,
    trait.Ef=sim.result$trait.Ef,
    stringsAsFactors=FALSE)
  
  # -- diversity stats for order q, may be calculated for more than 1 value of q
  dat.divpart<-list()
  dat.varpart<-list()
  
  if(calc.comm.stats){
    for(i.q in q.order){
    
        # diversity partition for order i.q
        dat.divpart[[as.character(i.q)]]<-ldply(
          .data=sim.result$J[keep.timesteps],
          .fun=fn.diversity.partition,
          q.order=i.q)
    
		# variation partition for order i.q
        dat.varpart[[as.character(i.q)]]<-ldply(
          .data=sim.result$J[keep.timesteps],
          .fun=fn.variation.partition,
          E=data.frame(Ef=sim.result$Ef),
          S=sim.result$landscape$dat[,sim.result$landscape$pcnm.list],
          q.order=i.q)
    }
  }

  sim.result.metadata<-data.frame(
    row.names=paste(sim.result.name,"t",keep.timesteps,sep="_"),
    scenario.ID=sim.result$scenario.ID,
    sim.ID=sim.result.name,
    timestep=keep.timesteps,
    sim.result.filename=sim.result.filename,
    n.sites=length(sim.result$x),
    n.timestep=sim.result$n.timestep,
    alpha.fisher=sim.result$alpha.fisher,
    nu.sim=sim.result$nu.sim,
    JM=sim.result$landscape$JM,
    JL.mean=mean(sim.result$JL),
    JL.sd=sd(sim.result$JL),
    JL.scale=sim.result$landscape$JL.scale,
    JL.pcnm=sim.result$landscape$JL.pcnm,
    m.mean=mean(sim.result$landscape$dat$m),
    m.sd=sd(sim.result$landscape$dat$m),
    IL.mean=mean(sim.result$landscape$dat$IL),
    IL.sd=sd(sim.result$landscape$dat$IL),
    IL.scale=sim.result$landscape$IL.scale,
    IL.pcnm=sim.result$landscape$IL.pcnm,
    SWM.slope=sim.result$SWM.slope,
    Ef.mean=mean(sim.result$landscape$dat$Ef),
    Ef.sd=sd(sim.result$landscape$dat$Ef),
    Ef.scale=sim.result$landscape$Ef.scale,
    Ef.pcnm=sim.result$landscape$Ef.pcnm,
    Ef.specificity=sim.result$Ef.specificity,
    Tr.disp.mean=mean(dat.trait$trait.dispersal),
    Tr.disp.sd=sd(dat.trait$trait.dispersal),
    Niche.breadth=sim.result$trait.Ef.sd,
    stringsAsFactors=FALSE)

  if(calc.comm.stats){
	  # -- turn list into table
	  divpart.mat<-do.call(cbind,dat.divpart)
	  names(divpart.mat)<-unlist(lapply(dat.divpart,names))
	  
	  varpart.mat<-do.call(cbind,dat.varpart)
	  name.mat<-expand.grid(
		names(dat.varpart[[1]]),
		paste("q",q.order,sep=""))
	  names(varpart.mat)<-apply(name.mat,1,paste,collapse=".")
	  
	  # -- combine parameters and results into one table
	  sim.result.metadata<-cbind(sim.result.metadata,divpart.mat,varpart.mat)
  }
  
  # -- write metadata to a table, save as csv
  # save parameter settings and sim.result name
  filname.sim.metadata<-paste(var.dir,"/sim.metadata_",sim.result$scenario.ID,".csv",sep="")
  
  if(file.exists(filname.sim.metadata)){
    dat.sim.metadata<-read.csv(filname.sim.metadata,
                               row.names=1,header=TRUE)
    if(!sim.result.name%in%row.names(dat.sim.metadata)){
      dat.sim.metadata<-rbind(dat.sim.metadata,sim.result.metadata)
      write.csv(dat.sim.metadata,filname.sim.metadata)
    }
  }else{
    write.csv(sim.result.metadata,filname.sim.metadata)
  }
  
  # -- save sim.result as a .rda file
  if(save.sim) save(sim.result,file=sim.result.filename)
}
