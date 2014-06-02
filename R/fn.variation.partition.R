#' fn.variation.partition
#' 
#' @title Variation partitioning using dbRDA
#' 
#' @description Variation partitioning, similar to \link[vegan]{varpart} in 
#' the \pkg{vegan} package, but uses \link[vegan]{capscale} (dbRDA) to 
#' calculate [a], [b], [c], and [d] components instead of RDA.  Current 
#' implementation only calculates species dissimilarities (of any order q) based 
#' on Jost 2006, Jost 2007, and Chao et al. 2012 formula for relative turnover 
#' (pairwise "true" beta diversities).
#' 
#' @usage 
#' fn.variation.partition(dat.comm, E, S, q.order = 1)
#'  
#' @param dat.comm community matrix
#' @param E A data frame in which columns are environemntal vars -- currently, the 
#' simulation only produces one environmental gradient.
#' @param S A data frame of spatial variables (i.e., PCNM eigenvectors produced 
#' by \link[vegan]{pcnm} in the \pkg{vegan} package.
#' @param q.order Order q to use to calculate among site dissimilarities in community 
#' composition.  See \link{fn.true.dist}.
#' 
#' @references 
#' Borcard, D., and P. Legendre. 2002. All-scale spatial analysis of ecological data 
#' by means of principal coordinates of neighbour matrices. Ecological Modelling 
#' 153:51--68.
#' 
#' Peres-Neto, P. R., P. Legendre, S. Dray, and D. Borcard. 2006. Variation 
#' partitioning of species data matrices - estimation and comparison of fractions. 
#' Ecology 87:2614--2625.
#' 
#' @export
#' 
fn.variation.partition <- function(
  dat.comm,
  E,
  S,
  q.order=1){

  require(vegan)
  
  # -- remove unobserved species
  dat.comm<-dat.comm[,colSums(dat.comm)>0]
  
  # -- Variation partitioning
  S.formula<-"ALL"
  
  # -- calc dists based on CqN overlap
  if(q.order==0){
    dat.comm.dist<-vegdist(dat.comm,"bray",binary=TRUE) #same result as fn.true.dist, but vegdist is faster
  }else if(q.order==2){
    dat.comm.dist<-vegdist(dat.comm,"horn")
  }else{
    dat.comm.dist<-as.dist(fn.true.dist(d=dat.comm,q.order=q.order)) #same result as fn.true.dist, but vegdist is faster
  }
  
  # -- check for postivit eigenvalues in comm matrix
  eig.comm<-eigen(dat.comm.dist)
  eig.check<-sum(eig.comm$values[eig.comm$values>0])>0
  
  # -- select significant PCNM vars
  fn.sig.dbRDA<-function(X=S[,1],Y=dat.comm.dist){
    require(vegan)
    mod.anova<-anova(capscale(Y~X,na.action="na.omit"))
    return(mod.anova["Model","Pr(>F)"])
  }
  
  if(eig.check){
    pcnm.pvals<-apply(X=S,
                      MARGIN=2,
                      FUN=fn.sig.dbRDA,
                      Y=dat.comm.dist)
    
    pcnm.pvals[is.na(pcnm.pvals)]<-1 #set NAs as non-significant
    
    pcnm.sig.list<-names(pcnm.pvals[pcnm.pvals<=0.05])
    if(length(pcnm.sig.list)==0) pcnm.sig.list<-names(pcnm.pvals[pcnm.pvals==min(pcnm.pvals)])
    
    if(nrow(S)>length(pcnm.sig.list)) S.formula<-paste(pcnm.sig.list,collapse="+")
    
    if(length(pcnm.sig.list)>1){
      S<-S[,pcnm.sig.list]
    }else{
      S<-data.frame(x=S[,pcnm.sig.list])
      names(S)<-pcnm.sig.list
    }
    
    # -- calculate dbRDA models
    mod.ab<-capscale(dat.comm.dist~.,data=E,na.action="na.omit")
    mod.bc<-capscale(dat.comm.dist~.,data=S,na.action="na.omit")
    mod.abc<-capscale(dat.comm.dist~.,data=cbind(E,S),na.action="na.omit")
    
    ab<-RsquareAdj(mod.ab)$adj.r.squared
    #   p.ab<-anova(mod.ab)["Model","Pr(>F)"]
    bc<-RsquareAdj(mod.bc)$adj.r.squared
    #   p.bc<-anova(mod.bc)["Model","Pr(>F)"]
    abc<-RsquareAdj(mod.abc)$adj.r.squared
    #   p.abc<-anova(mod.abc)["Model","Pr(>F)"]
    
    a<-abc-bc
    c<-abc-ab
    b<-ab+bc-abc
    d<-1-abc
    
    var.part.results<-data.frame(a=a,b=b,c=c,d=d,
                                 E=ab,
                                 S=bc,
                                 ES=abc,
                                 S.vars=S.formula)
  }else{
    var.part.results<-data.frame(a=NA,b=NA,c=NA,d=NA,
                                 E=NA,
                                 S=NA,
                                 ES=NA,
                                 S.vars=NA)
  }
  
  return(var.part.results)
}
