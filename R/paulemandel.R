paulemandel <- function(TE, seTE,
                        tol=.Machine$double.eps^0.25,
                        maxiter=25, ...){
  
  ##
  ## Mandel-Paule algorithm
  ## Based on R function mpaule.default from R package metRology 
  ## Author: S. Cowen <simon.cowen@lgc.co.uk> with amendments by
  ##         S.L.R. Ellison <s.ellison@lgc.co.uk>
  ##

  sel <- !is.na(TE) & !is.na(seTE)
  TE <- TE[sel]
  seTE <- seTE[sel]
  ##
  if (length(TE)==0)
    stop("Mandel-Paule estimate not defined as all studies have missing values in argument TE and/or seTE.")
  
  if (tol >= 1.0 )
    stop("Tolerance must be smaller than 1.0.")
  
  ## Guarantees TE.random exists in case no iterations are run
  TE.random <- NA
  
  tau2 <- variance.TE <- dv <- var(TE)
  ##
  n.iter <- 0
  converged <- 0L
  ##
  while (n.iter < maxiter && abs(dv) > tol*variance.TE){
    n.iter <- n.iter + 1
    converged <- 0L
    ##
    w.random <- 1 / (seTE^2 + tau2)
    w.random[is.na(w.random)] <- 0
    TE.random <- weighted.mean(TE, w.random)
    ##
    F <- sum(w.random*(TE - TE.random)^2) - (length(TE) - 1)
    dv <- F / sum(w.random^2 * (TE - TE.random)^2)
    tau2 <- tau2 + dv
    ##
    if (tau2 < 0){
      tau2 <- 0.0
      w.random <- 1 / (seTE^2 + tau2)
      w.random[is.na(w.random)] <- 0
      TE.random <- weighted.mean(TE, w.random)
      converged <- 2L
    }
  }
  ##
  seTE.random <- 1/sqrt(sum(w.random))
  
  if(abs(dv) >= tol*variance.TE){
    warning("Maximum iterations reached; Mandel-Paule algorithm may not have converged.")
  }
  else{
    ## Not changed if already set to 2L
    if (converged==0L)
      converged <- 1L
  }
  
  w.random.all <- rep(0, length(sel))
  w.random.all[sel] <- w.random
  
  res <- list(TE.random=TE.random, seTE.random=seTE.random,
              w.random=w.random.all, tau=sqrt(tau2),
              n.iter=n.iter, converged=converged, tol=tol*variance.TE)
  
  res
}
