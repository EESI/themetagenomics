jsd <- function(p,q) {

  m <- .5*(p + q)
  div <- .5*sum(p*log(p/m)) + .5*sum(q*log(q/m))
  dis <- sqrt(div)

  return(dis)
}

norm10 <- function(x) (x-min(x))/(max(x)-min(x))
