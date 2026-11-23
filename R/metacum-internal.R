#' @method summary metacum
#' @export

summary.metacum <- function(object, ...)
  object

select_interval_method <- function(meths, adjusts,
                                   method.interval, missing.method.interval,
                                   adhoc.hakn, missing.adhoc.hakn) {
  name1 <- deparse(substitute(method.interval))
  name2 <- deparse(substitute(adhoc.hakn))
  #
  if (missing.method.interval)
    method.interval <- meths[1]
  else
    method.interval <- setchar(method.interval, meths, length = 1, name = name1)
  #
  if (missing.adhoc.hakn)
    adhoc.hakn <- adjusts[which(meths == method.interval)[1]]
  else {
    adhoc.hakn <-
      setchar(adhoc.hakn, adjusts[meths == method.interval],
              length = 1, name = name2)
  }
  #
  selected <- meths == method.interval & adjusts == adhoc.hakn
  #
  res <- which(selected)[1]
  #
  res
}
