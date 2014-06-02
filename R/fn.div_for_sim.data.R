#' fn.div_for_sim.data
#' 
#' @title Simulation diversity partitioning
#' 
#' @description Read metacommunity simulation outcomes from a specified directory, then
#' calculate diversity partitions and write output to a .csv file in the same directory.
#' 
#' @param i.sim identify simulation scenario to read
#' @param out.dir.name identify directory from which to read simulation outcomes
#' @param managed A vector identifying sites in a metacommunity simulation as "managed"
#' or "unmanaged"
#' 
#' @export
#' 
fn.div_for_sim.data<-function(
  i.sim,
  out.dir.name,
  managed = site.management){
  
  #' ---------------------------------------------------------------------------------
  #' -- check to see if diversities have already been calculated and stored in out.dir.name
  #' ---------------------------------------------------------------------------------
  filenames.in.dir<-list.files(out.dir.name)
  div.metric.name.check<-paste("div.metrics.sim.data_",i.sim,".csv",sep="")
  div.metric.name.temp<-paste(out.dir.name,"/","div.metrics.sim.data_",i.sim,".csv",sep="")

  while(!div.metric.name.check%in%filenames.in.dir){
    #' ---------------------------------------------------------------------------------
    #' -- load output data from simulation
    #' ---------------------------------------------------------------------------------
    file.name.temp<-paste(out.dir.name,"/","sim.data_",i.sim,".rda",sep="")
    load(file.name.temp)
    
    dat.comm<-sim.result.list$dat.comm
    dat.siteinfo<-sim.result.list$dat.landscape
        
    dat.div.summary<-data.frame(stringsAsFactors=FALSE)
    monodominant<-!(sum(colSums(dat.comm)>0))>1
    
    #' get managed vector from siteinfo if it's not provided
    if(length(managed)!=nrow(dat.comm)){
      managed<-dat.siteinfo[,agrep('manage',names(dat.siteinfo))]
    }
    
    #' ---------------------------------------------------------------------------------
    #' ---------------------------------------------------------------------------------
    #' -- divpart q = 0
    dat.divpart.temp<-do.call(
      rbind,
      by(dat.comm,managed,fn.diversity.partition,q.order=0))
    names(dat.divpart.temp)<-c("alpha","beta","gamma")
    
    dat.divpart.all<-fn.diversity.partition(dat.comm,q.order=0)
    names(dat.divpart.all)<-c("alpha","beta","gamma")
    dat.divpart.all<-data.frame(hypothesis=out.dir.name,
               scenario.ID=i.sim,
               management.type='all',
               dat.divpart.all,
               stringsAsFactors=FALSE)
      
    diff.divpart<-dat.divpart.temp[2,]-dat.divpart.temp[1,]
    
    change.divpart.q0<-data.frame(hypothesis=out.dir.name,
                                  scenario.ID=i.sim,
                                  management.type='diff',
                                  diff.divpart,
                                  stringsAsFactors=FALSE)
    
    dat.div.summary<-rbind(
      data.frame(hypothesis=out.dir.name,
                 scenario.ID=i.sim,
                 management.type=ifelse(
                   row.names(dat.divpart.temp)==0,
                   'unmanaged','managed'),
                 rbind(dat.divpart.temp),
                 stringsAsFactors=FALSE),
      dat.divpart.all,
      change.divpart.q0)
   
    row.names(dat.div.summary)<-NULL
    
    write.csv(dat.div.summary,file=div.metric.name.temp)
    filenames.in.dir<-list.files(out.dir.name)
  }
}