catmeth <- function(method, method.tau=NULL,
                    sm="", k.all,
                    hakn=FALSE, metaprop=FALSE,
                    trimfill=FALSE,
                    tau.common=FALSE,
                    tau.preset=NULL,
                    metabin=FALSE,
                    metainc=FALSE,
                    sparse=FALSE,
                    incr=NULL,
                    allincr=FALSE,
                    addincr=FALSE,
                    MH.exact=FALSE,
                    method.ci=NULL,
                    metacont=FALSE,
                    pooledvar=FALSE,
                    method.smd,
                    sd.glass,
                    exact.smd=FALSE){
  
  if  (sm=="PFT")
    sm.details <- "\n- Freeman-Tukey double arcsine transformation"
  else if (sm=="PAS")
    sm.details <- "\n- Arcsine transformation"
  else if (sm=="PLN")
    sm.details <- "\n- Log transformation"
  else if (sm=="PLOGIT")
    sm.details <- "\n- Logit transformation"
  else if (sm=="PRAW")
    sm.details <- "\n- Untransformed proportions"
  else if (sm=="ZCOR")
    sm.details <- "\n- Fisher's z transformation of correlations"
  else if (sm=="COR")
    sm.details <- "\n- Untransformed correlations"
  else
    sm.details <- ""
  ##
  if (metaprop && !is.null(method.ci)){
    if  (method.ci=="CP")
      method.ci.details <- "\n- Clopper-Pearson confidence interval for individual studies"
    else if (method.ci=="WS")
      method.ci.details <- "\n- Wilson Score confidence interval for individual studies"
    else if (method.ci=="WSCC")
      method.ci.details <- "\n- Wilson Score confidence interval with continuity correction for individual studies"
    else if (method.ci=="AC")
      method.ci.details <- "\n- Agresti-Coull confidence interval for individual studies"
    else if (method.ci=="SA")
      method.ci.details <- "\n- Simple approximation confidence interval for individual studies"
    else if (method.ci=="SACC")
      method.ci.details <- "\n- Simple approximation confidence interval with continuity correction for individual studies"
    else if (method.ci=="NAsm")
      method.ci.details <- "\n- Normal approximation confidence interval for individual studies"
    ##
    sm.details <- paste(sm.details, method.ci.details, sep="")
  }
  ##
  if (metacont && sm=="SMD" && !is.null(method.smd)){
    if  (method.smd=="Hedges")
      if (exact.smd)
        method.details <- "\n- Hedges' g (bias corrected standardised mean difference; using exact formulae)"
      else
        method.details <- "\n- Hedges' g (bias corrected standardised mean difference)"
    else if  (method.smd=="Cohen")
      if (exact.smd)
        method.details <- "\n- Cohen's d (standardised mean difference; using exact formulae)"
      else
        method.details <- "\n- Cohen's d (standardised mean difference)"
    else if  (method.smd=="Glass"){
      if (!is.null(sd.glass) && sd.glass=="control")
        method.details <- "\n- Glass' delta (standardised mean difference; based on control group)"
      else if (!is.null(sd.glass) && sd.glass=="experimental")
        method.details <- "\n- Glass' delta (standardised mean difference; based on experimental group)"
    }
    ##
    sm.details <- paste(sm.details, method.details, sep="")
  }
  ##
  if (metabin | metainc | metaprop){
    if (!(sm=="ASD" | method=="Peto")){
      if (addincr){
        if (incr=="TACC")
          sm.details <- paste(sm.details,
                              "\n- Treatment arm continuity correction in all studies",
                              sep="")
        else if (incr!=0)
          sm.details <- paste(sm.details,
                              "\n- Continuity correction of ", round(incr, 4),
                              " in all studies", sep="")
        else
          sm.details <- ""
      }
      else if (sparse){
        if (allincr==FALSE)
          if (incr=="TACC")
            sm.details <- paste(sm.details,
                                "\n- Treatment arm continuity correction in studies with zero cell frequencies",
                                sep="")
          else if (incr!=0)
            sm.details <- paste(sm.details,
                                "\n- Continuity correction of ", round(incr, 4),
                                " in studies with zero cell frequencies", sep="")
          else
            sm.details <- ""
        else
          if (incr=="TACC")
            sm.details <- paste(sm.details,
                                "\n- Treatment arm continuity correction in all studies",
                                sep="")
          else if (incr!=0)
            sm.details <- paste(sm.details,
                                "\n- Continuity correction of ", round(incr, 4),
                                " in all studies", sep="")
          else
            sm.details <- ""
      }
    }
  }
  
  
  lab.method.details <- ""
  ##
  if (is.null(method.tau))
    lab.method.tau <- ""
  else {
    if (!is.null(tau.preset)){
      tau2 <- tau.preset^2
      if (tau2 > 0 & tau2 < 0.0001)
        tau2 <- paste("tau^2", format.tau(tau2))
      else
        tau2 <- paste("tau^2 = ",
                      ifelse(tau2==0,
                             "0",
                             format(round(tau2, 4), 4, nsmall=4, scientific=FALSE)),
                      sep="")
      ##
      lab.method.tau <- paste("\n- Preset between-study variance: ",
                              tau2, sep="")
      ##
      lab.method.details <- lab.method.tau
    }
    else{
      i.lab.method.tau <- charmatch(method.tau,
                                    c("DL", "PM", "REML", "ML", "HS", "SJ", "HE", "EB"),
                                    nomatch = NA)
      ##
      lab.method.tau <- c("\n- DerSimonian-Laird estimator for tau^2",
                          "\n- Paule-Mandel estimator for tau^2",
                          "\n- restricted maximum-likelihood estimator for tau^2",
                          "\n- maximum-likelihood estimator for tau^2",
                          "\n- Hunter-Schmidt estimator for tau^2",
                          "\n- Sidik-Jonkman estimator for tau^2",
                          "\n- Hedges estimator for tau^2",
                          "\n- empirical Bayes estimator for tau^2")[i.lab.method.tau]
      ##
      if (tau.common)
        lab.method.tau <- paste(lab.method.tau, " (assuming common tau^2 in subgroups)", sep="")
      ##
      if (hakn)
        lab.hakn <- "\n- Hartung-Knapp adjustment for random effects model"
      else
        lab.hakn <- ""
      ##      
      lab.method.details <- paste(lab.method.tau, lab.hakn, sep="")
    }
  }
  ##
  imeth <- charmatch(method,
                     c("MH", "Peto", "Inverse", "Cochran"),
                     nomatch = NA)
  ##
  if ((metabin|metainc) & imeth==1 & (sparse | addincr))
    if (MH.exact | metainc)
      lab.method.details <- paste(" (without continuity correction)",
                                  lab.method.details, sep="")
  ##
  if (metacont && !is.null(pooledvar) && pooledvar)
    lab.method.details <- paste(" (with pooled variance for individual studies)",
                                lab.method.details, sep="")
  ##
  method <- c("\n- Mantel-Haenszel method",
              "\n- Peto method",
              "\n- Inverse variance method",
              "\n- Cochran method")[imeth]
  method <- paste(method, lab.method.details, sep="")
  ##
  if (k.all > 1){
    cat(paste("\nDetails on meta-analytical method:", method, sep=""))
    if (trimfill)
      cat("\n- Trim-and-fill method to adjust for funnel plot asymmetry")
  }
  else
    cat(paste("\nDetails:", method, sep=""))
  ##
  if (sm.details!="")
    cat(sm.details,
        "\n",
        sep="")
  else
    cat("\n")
  
  invisible(NULL)
}
