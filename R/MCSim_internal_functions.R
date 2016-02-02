# --------------------------------------------------------
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

count <- function(a){
  if (length(a) > 0) {
    out <- sort(as.table(a), decreasing = TRUE)
  }
  else {
    out <- integer(0)
  }
  class(out) <- c("count", "table")
  return(out)
}

is.count <- function(a){
  inherits(a, "count")
}
# --------------------------------------------------------
#' fisher.ecosystem
#' 
#' @title fisher.ecosystem
#' @aliases fisher.ecosystem
#' 
#' @description Function from the \pkg{untb} package, 
#' see \link[untb]{fisher.ecosystem}  -- used for internal 
#' calculations in \link{fn.metaSIM}.
#' 
#' @param N Number of individuals
#' @param S Number of species
#' @param nmax Max number of species abundance classes
#' @param alpha Fisher's alpha. If not supplied it will be calculated from \code{N} and \code{S}
#' @param c Rare species advantage term
#' 
#' @references Hankin, R. K. 2007. Introducing untb, an R package for simulating 
#' ecological drift under the unified neutral theory of biodiversity. Journal of 
#' Statistical Software 22:1--15. 
#' \url{http://cran.r-project.org/web/packages/untb/index.html}
#' 
#' @seealso \code{\link[untb]{fisher.ecosystem}}
#' 
fisher.ecosystem <- function(N, S, nmax, alpha = 2, c = 0){
  x <- N/(N + alpha)
  j <- 1:nmax
  jj <- rpois(n = nmax, lambda = alpha * x^j/(j + c))
  out <- as.count(rep(1:sum(jj), rep(j, jj)))
  names(out) <- paste("sp", names(out), sep = ".")
  return(out)
}

# --------------------------------------------------------
