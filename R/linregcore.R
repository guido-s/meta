linregcore <- function(x, y, w=NULL){
  ##
  ## core function for method.bias linreg and mm
  ##
  ##
  ## Check:
  ##
  if(length(x) != length(y))
    stop("length of argument x and y must be equal")
  ##
  if (!is.null(w)){
    if(length(x) != length(w))
      stop("length of argument x and w must be equal")
  }
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
  if (is.null(w))
    w <- rep(1, n)
  else
    w <- w[sel]
  ##
  W <- diag(w)
  ##
  if (n > 2){
    X <- cbind(rep(1, n), x)
    XWX <- solve(t(X) %*% W %*% X)
    ##
    coefs <- XWX %*% t(X) %*% W %*% y
    MSE.w <- 1/(n-2) * sum(w*(y - X %*% coefs)^2)
    df <- n - 2
    se.coefs <- sqrt(MSE.w * diag(XWX))
    ##
    intercept <- coefs[1]
    se.intercept <- se.coefs[1]
    slope <- coefs[2]
    se.slope <- se.coefs[2]
  }
  else {
    intercept <- NA
    se.intercept <- NA
    slope <- NA
    se.slope <- NA
  }
  ##
  res <- list(intercept=intercept,
              se.intercept=se.intercept,
              slope=slope,
              se.slope=se.slope,
              df=df,
              MSE.w=MSE.w)
  res
}
