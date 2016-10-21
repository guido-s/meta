Rb <- function(seTE, seTE.random, tau2, Q, df.Q, level) {
  ##
  ## Calculate Rb
  ## Crippa et al. (2016), Statistics in Medicine, 35 3661-75.
  ##
  
  
  if (is.na(df.Q) || df.Q == 0)
    return(list(TE = NaN, lower = NaN, upper = NaN))
  ##
  k <- df.Q + 1
  varTE <- seTE^2
  w.fixed <- 1 / varTE
  
  
  ##
  ## Equation (4) in Crippa et al. (2016)
  ##
  ## Rb <- 1 / k * sum(tau2 / (varTE + tau2))
  Rb <- tau2 / (k * seTE.random^2)
  
  
  ##
  ## Appendix D in Crippa et al. (2016)
  ##
  S <- function(n, w) sum(w^n)
  ##
  a <- varTE * (S(1, w.fixed) - S(2, w.fixed) / S(1, w.fixed))
  ##
  seQ <- sqrt(2 * (k - 1) +
              4 * (S(1, w.fixed) - S(2, w.fixed) / S(1, w.fixed)) * tau2 +
              2 * (S(2, w.fixed) - 2 * S(3, w.fixed) / S(1, w.fixed) +
                   S(2, w.fixed)^2 / S(1, w.fixed)^2) * tau2^2
              )
  ##  
  seRb <- sqrt((1 / k * sum(a / (Q + a - (k - 1))^2))^2 * seQ^2)
  
  cint <- ci(Rb, seRb, level)
  ##
  cint$lower <- ifelse(cint$lower < 0, 0, cint$lower)
  cint$upper <- ifelse(cint$upper > 1, 1, cint$upper)
  
  res <- list(TE = cint$TE, lower = cint$lower, upper = cint$upper)
  ##
  res
}
