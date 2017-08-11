trimfill.default <- function(x, seTE, left = NULL, ma.fixed = TRUE,
                             type = "L", n.iter.max = 50,
                             sm = NULL, studlab = NULL,
                             level = 0.95, level.comb = level,
                             comb.fixed = FALSE, comb.random = TRUE,
                             hakn = FALSE,
                             method.tau = "DL",
                             prediction = FALSE, level.predict = level,
                             backtransf = TRUE, pscale = 1,
                             irscale = 1, irunit = "person-years",
                             silent = TRUE, ...) {
  
  
  ##
  ##
  ## (1) Check essential arguments
  ##
  ##
  k.All <- length(x)
  ##
  chknumeric(x)
  chknumeric(seTE)
  ##
  fun <- "trimfill"
  chklength(seTE, k.All, fun)
  ##
  if (!is.null(studlab))
    chklength(studlab, k.All, fun)
  else
    studlab <- seq(along = x)
  ##
  if (is.null(sm)) sm <- ""
  
  
  ##
  ##
  ## (2) Do meta-analysis
  ##
  ##
  m <- metagen(x, seTE, studlab = studlab, sm = sm)
  
  
  ##
  ##
  ## (3) Run trim-and-fill method
  ##
  ##
  res <- trimfill(m, left = left, ma.fixed = ma.fixed,
                  type = type, n.iter.max = n.iter.max,
                  level = level, level.comb = level.comb,
                  comb.fixed = comb.fixed, comb.random = TRUE,
                  hakn = hakn,
                  method.tau = method.tau,
                  prediction = prediction, level.predict = level.predict,
                  backtransf = backtransf, pscale = pscale,
                  irscale = irscale, irunit = irunit,
                  silent = silent, ...)
  
  
  res
}
