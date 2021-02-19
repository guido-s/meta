linregcore <- function(TE, seTE, covar = NULL,
                       model = "lm", method.tau = "DL") {
  
  if (is.null(covar))
    predictor <- "sei"
  else
    predictor <- "ni"
  
  rma1 <- suppressWarnings(rma.uni(TE, sei = seTE, ni = covar,
                                   method = method.tau, test = "t"))
  ##
  reg <- suppressWarnings(regtest(rma1,
                                  ni = covar, predictor = predictor,
                                  model = model))
  creg <- suppressWarnings(coef(reg$fit))
  ##
  if (is.matrix(creg))
    res <- list(intercept = creg[1, 1], se.intercept = creg[1, 2],
                slope = creg[2, 1], se.slope = creg[2, 2])
  else
    res <- list(intercept = reg$fit$beta[1], se.intercept = reg$fit$se[1],
                slope = reg$fit$beta[2], se.slope = reg$fit$se[2])
  ##
  res$statistic <- res$slope / res$se.slope
  res$df <- reg$dfs
  res$pval <- 2 * pt(abs(res$statistic), df = res$df, lower.tail = FALSE)
  ##
  if (model == "rma")
    res$tau <- sqrt(reg$fit$tau2)
  else
    res$tau <- reg$fit$sigma
  ##
  res$TE <- TE
  res$seTE <- seTE
  res$covar <- covar
  res$predictor <- predictor
  res$model <- model
  res$method.tau <- method.tau
  ##
  names(res$statistic) <- "t"
  
  
  res
}
