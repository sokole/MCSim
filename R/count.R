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
