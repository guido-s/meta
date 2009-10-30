addvar <- function(x, y, varname,
                   by.x="studlab", by.y=by.x){
  ##
  if (!inherits(x, "meta"))
    stop("Argument 'x' must be an object of class \"meta\"")
  ##
  if (!is.character(varname))
    stop("parameter 'varname' must be a character")
  if (length(varname)!=1)
    stop("parameter 'varname' must be of length 1")
  ##
  if (length(x[[varname]])==0){
    tdata <- as.data.frame(x)
    tdata$orderorder <- 1:dim(tdata)[1]
    ##
    mdata <- merge(tdata, y,
                   by.x=by.x, by.y=by.y,
                   all.x=TRUE, all.y=FALSE,
                   sort=FALSE)
    ##
    res <- as.vector(unlist(mdata[varname]))[mdata$orderorder]
  }
  else
    stop(paste("variable '", varname, "' already exists in '",
               deparse(substitute(x)), "'", sep=""))
  ##
  res
}
