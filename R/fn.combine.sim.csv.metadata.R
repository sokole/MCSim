#' fn.combine.sim.csv.metadata
#' 
#' @title fn.combine.sim.csv.metadata
#' 
#' @description Function from the \pkg{untb} package, 
#' see \link[untb]{fisher.ecosystem}  -- used for internal 
#' calculations in \link{fn.metaSIM}.
#' 
#' @param out.dir.name Name of the target directory in which simulation 
#' results are saved, including the .csv metadata files associated with each simulation.  
#' @param file.name currently must be \code{sim.metadata_}.  This is the file name used 
#' by \code{fn.metaSIM} for the simulation metadata files.
#' @param n.sims The number of simulations saved to the target directory.
#' 
#' @export
#' 
fn.combine.sim.csv.metadata <- function(
  out.dir.name=out.dir.name,
  file.name="sim.metadata_",
  n.sims=nrow(data.simparameters)){
  require(plyr)
  if((length(n.sims)==1)&(is.numeric(n.sims))){
	scenario.list<-c(1:n.sims)
  }else{
	scenario.list<-n.sims
  }
  dat.list<-alply(
    .data=scenario.list,
    .margins=1,
    .fun=function(i.sim,out.dir.name,file.name){
      filename.temp<-paste(out.dir.name,"/",file.name,i.sim,".csv",sep="")
      if(file.exists(filename.temp)){
        return(read.csv(file=filename.temp,
                        stringsAsFactors=FALSE,row.names=1))
      }else{
        return("no sim output")
      }
    },
    out.dir.name=out.dir.name,
    file.name=file.name)
  
  expected.cols<-max(unlist(lapply(dat.list,length)))
  
  if(expected.cols>0){
    dat.list.keepers<-dat.list[lapply(dat.list,length)==expected.cols]
    dat.all<-do.call(rbind,dat.list.keepers)
    write.csv(dat.all,file=paste(out.dir.name,"/sim.metadata.csv",sep=""))
  }else{
    print("warning: no metadata or results")
  }
}
