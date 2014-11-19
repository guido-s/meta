hetcalc <- function(TE, seTE, method.tau, TE.tau, byvar){
  
  Ccalc <- function(x){
    res <- (sum(x  , na.rm=TRUE) -
            sum(x^2, na.rm=TRUE)/
            sum(x  , na.rm=TRUE))
    ##
    res
  }
  
  
  if (!missing(byvar)){
    bylevs <- bylevs(byvar)
    ##
    if (method.tau=="DL"){
      res.w <- matrix(NA, ncol=3, nrow=length(bylevs))
      j <- 0
      ##
      for (i in bylevs){
        j <- j+1
        sel <- byvar == i
        ##
        if (all(is.na(TE[sel])))
          stop("No data available for byvar = ", i)
        ##
        m1 <- metagen(TE[sel], seTE[sel], method.tau="DL")
        ##
        res.w[j,] <- c(m1$Q, m1$k-1, m1$C)
      }
      ##
      Q    <- sum(res.w[,1])
      df.Q <- sum(res.w[,2])
      Cval <- sum(res.w[,3])
      ##
      if (round(Q, digits=18)<=df.Q)
        tau2 <- 0
      else
        tau2 <- (Q-df.Q)/Cval
      ##
      se.tau2 <- NULL
    }
    else{
      if (is.numeric(byvar))
        byvar <- as.factor(byvar)
      ##
      mf1 <- metafor::rma.uni(yi=TE, vi=seTE^2, method=method.tau, mods=~byvar)
      ##
      Q    <- mf1$QE
      df.Q <- mf1$k - mf1$p
      Cval <- NA
      ##
      tau2 <- mf1$tau2
      ##
      se.tau2 <- mf1$se.tau2
    }
  }
  else{
    if (method.tau=="DL" | method.tau=="PM"){
      w.fixed <- 1/seTE^2
      w.fixed[is.na(w.fixed)] <- 0
      ##
      TE.fixed   <- weighted.mean(TE, w.fixed, na.rm=TRUE)
      ##
      if (is.null(TE.tau))
        Q <- sum(w.fixed * (TE - TE.fixed)^2, na.rm=TRUE)
      else
        Q <- sum(w.fixed * (TE - TE.tau  )^2, na.rm=TRUE)
      ##
      df.Q <- sum(!is.na(seTE)) - 1
      Cval <- Ccalc(w.fixed)
      ##
      if (round(Q, digits=18)<=df.Q)
        tau2 <- 0
      else
        tau2 <- (Q-df.Q)/Cval
      ##
      se.tau2 <- NULL
    }
    else{
      mf2 <- metafor::rma.uni(yi=TE, vi=seTE^2, method=method.tau)
      Q    <- mf2$QE
      df.Q <- mf2$k - mf2$p
      Cval <- NA
      ##
      tau2 <- mf2$tau2
      ##
      se.tau2 <- mf2$se.tau2
    }
  }
  
  
  res <- list(tau=sqrt(tau2),
              se.tau2=se.tau2,
              Q=Q,
              df.Q=df.Q,
              Cval=Cval)
  ##
  res
}
