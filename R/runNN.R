runNN <- function(func, args) {
  args <- args[!sapply(args, is.null)]
  do.call(func, args)
}
