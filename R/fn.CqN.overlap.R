#' fn.CqN.overlap
#' 
#' @title Community overlap of order q
#' 
#' @description Function to calculate community overlap among two or more sites, as 
#' described by Jost et al. in Ch 6 in Biological Diversity (2011).
#' 
#' @param d Site by species matrix or data frame of abundances or relative abundances
#' @param q.order Order q used to calculate overlap.  0 is presence/absence (more 
#' sensitive to rare taxa), 1 is unbiased, and 2 is more sensitive to dominant taxa. 
#' 
#' @references Jost, L., A. Chao, and R. L. Chazdon. 2011. Compositional similarity and 
#' beta diversity. Pages 68--84 in A. E. Magurran and B. J. McGill, editors. Biological 
#' Diversity: Frontiers in Measurement and Assessment. Oxford University Press, Oxford, UK.
#' 
#' @export
#' 
fn.CqN.overlap <- function(
  d,
  q.order=1){
  if(q.order==0){ #eq for C0N, Jost et al. 2011, pp74 in Ch6 of Biological Diversity
    N<-nrow(d)
    S.bar<-mean(rowSums(d>0))
    S.total<-sum(colSums(d)>0)
    C<-(N - (S.total/S.bar) ) / (N - 1)
  }else if(q.order==1){ #eq for C1N, eq 6.4 on pp74
    p<-as.matrix(d/rowSums(d))
    N<-nrow(d)    
    p.ratio.sum<-apply(X=p,
                       FUN=function(x){
                         x.inv<-1/x
                         I.inv<-matrix(data=1,nrow=length(x),ncol=length(x))
                         diag(I.inv)<-0
                         x.sum<-colSums((x%o%x.inv)*I.inv)
                         x.sum[!is.finite(x.sum)]<-0
                         return(x.sum)
                       },
                       MARGIN=2)
    C<-(1/log(N)) * sum( (p/N)*log(1+p.ratio.sum) )
  }else{
    p<-data.frame(d/rowSums(d))
    #     p<-p[,colSums(p)!=0]
    N<-nrow(p)
    numerator.sum<-(1/( (N^q.order) - N)) * sum(((colSums(p))^q.order) - (colSums(p^q.order)))
    denominator.sum<-(1/N)*( sum(colSums(p^q.order)) )
    C<-numerator.sum/denominator.sum
    
  }
  return(C)
}
