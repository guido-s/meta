kentau <- function(x, y, correct=FALSE, keep.data=FALSE){
  ##
  ## Check:
  ##
  if(length(x) != length(y))
    stop("length of argument x and y must be equal")
  ##
  sel <- !is.na(x) & !is.na(y)
  if (length(x) != sum(sel))
    warning(paste(length(x) - sum(sel),
                  "observation(s) dropped due to missing values"))
  ##
  x <- x[sel]
  y <- y[sel]
  n <- length(x)
  ##
  ks <- .C("kenscore",
           kenscore = as.double(0),
           x = as.double(x),
           y = as.double(y),
           n = as.integer(n),
           PACKAGE="meta")$kenscore
  ##
  ## Calculate S and s.e(S) according to
  ## Stata, release 5, Reference P-Z, p.239-240
  ##
  ## see also Kendall, Gibbons (1990), Rank Correlation Methods
  ## p. 66-68
  ##
  t <- rle(sort(x))$lengths
  u <- rle(sort(y))$lengths
  ##
  ##
  se.ks <- sqrt(1/18 * (n*(n-1)*(2*n+5) -
                        sum(t*(t-1)*(2*t+5)) -
                        sum(u*(u-1)*(2*u+5))) +
                1/(9*n*(n-1)*(n-2)) *
                sum(t*(t-1)*(t-2))*
                sum(u*(u-1)*(u-2)) +
                1/(2*n*(n-1)) *
                sum(t*(t-1)) *
                sum(u*(u-1)))
  ##
  if (as.logical(correct) &
      any(c(length(unique(x)), length(unique(y)))==2))
    warning(paste("Continuity corrected statistic may be inappropriate,\n",
                  "see Kendall, Gibbons (1990), Rank Correlation Methods, ",
                  "p.67.\n ", sep=""))
  ##
  statistic <- (ks - sign(ks) * as.logical(correct)) / se.ks
  p.value <- 2*pnorm(abs(statistic), lower.tail=FALSE)
  ##
  N <- 0.5 * n * (n-1)
  N1 <- N-sum(t*(t-1)/2)
  N2 <- N-sum(u*(u-1)/2)
  ##
  res <- list(tau.a=ks/N,
              tau.b=ks/(sqrt(N1)*sqrt(N2)),
              ks=ks - sign(ks) * as.logical(correct),
              se.ks=se.ks,
              statistic=statistic,
              p.value=p.value,
              correction=as.logical(correct))
  ##
  if (keep.data){
    res$x <- x
    res$y <- y
  }
  ##
  res
}
