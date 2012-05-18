catmeth <- function(method, method.tau=NULL,
                    sm="", k.all,
                    hakn=FALSE, metaprop=FALSE,
                    trimfill=FALSE){
  
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
  if (metaprop)
    sm.details <- paste(sm.details,
                        "\n- Exact binomial confidence intervals for individual studies",
                        sep="")
  
  lab.method.details <- ""
  ##
  if (is.null(method.tau))
    lab.method.tau <- ""
  else {
    i.lab.method.tau <- charmatch(method.tau,
                                  c("DL", "REML", "ML", "HS", "SJ", "HE", "EB"),
                                  nomatch = NA)
    ##
    lab.method.tau <- c("\n- DerSimonian-Laird estimator for tau^2",
                        "\n- restricted maximum-likelihood estimator for tau^2",
                        "\n- maximum-likelihood estimator for tau^2",
                        "\n- Hunter-Schmidt estimator for tau^2",
                        "\n- Sidik-Jonkman estimator for tau^2",
                        "\n- Hedges estimator for tau^2",
                        "\n- empirical Bayes estimator for tau^2")[i.lab.method.tau]
    ##
    if (hakn)
      lab.hakn <- "\n- Hartung-Knapp adjustment for random effects model"
    else
      lab.hakn <- ""
    ##      
    lab.method.details <- paste(lab.method.tau, lab.hakn, sep="")
  }
  ##
  method <- ifelse(method=="MH",
                   "\n- Mantel-Haenszel method",
                   ifelse(method=="Peto",
                          paste("\n- Peto method", lab.method.details, sep=""),
                          ifelse(method=="Inverse",
                                 paste("\n- Inverse variance method",
                                       lab.method.details,
                                       sep=""),
                                 method)))
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
