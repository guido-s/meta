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
  if (length(x[[by.x]])==0)
    stop(paste("Variable '", by.x, "' not included in '",
               deparse(substitute(x)), "'", sep=""))
  ##
  if (length(y[[by.y]])==0)
    stop(paste("Variable '", by.y, "' not included in '",
               deparse(substitute(y)), "'", sep=""))
  ##
  if (length(x[[varname]])==0){
    tdata <- as.data.frame(x, stringsAsFactors=FALSE)
    if (
        (length(unique(tdata[[by.x]])) != length(tdata[[by.x]])) |
        (length(unique(y[[by.x]])) != length(y[[by.x]]))
        ){
      if (length(tdata[[by.x]]) != length(y[[by.y]])){
        if (by.x!=by.y)
          warning(paste("Duplicate entries in variable '", by.x,
                        "' or '", by.y,
                        "' and different length of these variables",
                        sep=""))
        else
          warning(paste("Duplicate entries in variables '", by.x,
                        "' and different length of these variables",
                        sep=""))
        res <- NULL
      }
      else{
        if (any(tdata[[by.x]] != y[[by.y]])){
          if (by.x!=by.y)
            warning(paste("Duplicate entries in variable '", by.x,
                          "' or '", by.y,
                          "' and values of these variables do not match",
                          sep=""))
          else
            warning(paste("Duplicate entries in variables '", by.x,
                          "' and values of these variables do not match",
                          sep=""))
        res <- NULL
        }
        else
          res <- y[[varname]]
      }
    }
    else{
      tdata$orderorder <- 1:dim(tdata)[1]
      ##
      mdata <- merge(tdata, y,
                     by.x=by.x, by.y=by.y,
                     all.x=TRUE, all.y=FALSE,
                     sort=FALSE)
      ##
      res <- as.vector(unlist(mdata[varname]))[mdata$orderorder]
    }
  }
  else
    stop(paste("Variable '", varname, "' already exists in '",
               deparse(substitute(x)), "'", sep=""))
  ##
  res
}
