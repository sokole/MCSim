#' count
#' 
#' @title count
#' @aliases as.count
#' @aliases is.count
#' 
#' @description Function from the \pkg{untb} package, see \link[untb]{count},
#' used for internal calculations in \link{fn.metaSIM}.
#' 
#' @param a species count data, see \link[untb]{count} documentation
#' @param add see \link[untb]{count} documentation
#' 
#' @references Hankin, R. K. 2007. Introducing untb, an R package for simulating 
#' ecological drift under the unified neutral theory of biodiversity. Journal of 
#' Statistical Software 22:1--15. 
#' \url{http://cran.r-project.org/web/packages/untb/index.html}
#' 
#' @seealso \code{\link[untb]{count}}
#' 
as.count <- function (a, add = ""){
  if (is.count(a)) {
    out <- a
  }
  else if (is.data.frame(a)) {
    if (nrow(a) != 1) {
      stop("data frame supplied: must have only one row")
    }
    out <- as.numeric(a)
    names(out) <- colnames(a)
    out <- count(out[out > 0])
  }
  else if (is.table(a)) {
    out <- count(a)
  }
  else {
    out <- count(table(a))
  }
  if (length(out) > 0) {
    names(out) <- paste(add, names(out), sep = "")
  }
  return(out)
}
