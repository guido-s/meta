catarg <- function(x, newline = TRUE, end = "") {
  xname <- x
  x <- gsub(" ", "", x)
  ##
  if (newline)
    cat("- ")
  ##
  if (is.null(.settings[[x]]))
    cat(paste0(xname, ' = NULL', end, '\n'))
  else if (is.character(.settings[[x]]))
    cat(paste0(xname, ' = "', .settings[[x]], '"', end, '\n'))
  else
    cat(paste0(xname, ' = ', .settings[[x]], end, "\n"))
  invisible(NULL)
}


specificSettings <- function(args, new, setting, quietly = FALSE) {
  isnull.old <- as.vector(unlist(lapply(.settings[args], is.null)))
  ischar.old <- as.vector(unlist(lapply(.settings[args], is.character)))
  old <- as.vector(unlist(.settings[args]))
  ##
  ischar.new <- as.vector(unlist(lapply(new, is.character)))
  new <- as.vector(unlist(new))
  ##
  label.old <- ifelse(isnull.old, "NULL",
               ifelse(ischar.old, paste0("\"", old, "\""), old))
  label.new <- ifelse(ischar.new, paste0("\"", new, "\""), new)
  ##
  sel <- new != old
  if (any(sel)) {
    tdata <- data.frame(argument = c("Argument",
                                     "--------",
                                     args[sel]),
                        space1 = rep("  ", along = c(1:2, sel)),
                        new = c("New value",
                                "---------",
                                label.new[sel]),
                        space2 = rep("  ", along = c(1:2, sel)),
                        previous = c("Previous value",
                                     "--------------",
                                     label.old[sel]))
    
    names(tdata) <- c("--------", "", "---------",
                      "", "--------------")
    ##
    if (!quietly) {
      cat(paste0("\n** Use ", setting, " (R package meta) **\n\n"))
      prmatrix(tdata, quote = FALSE, right = FALSE,
               rowlab = rep_len("", 2 + sum(sel)))
    }
    ##
    for (i in seq(along = args)) {
      new.i <- new[i]
      if (!ischar.new[i]) {
        if (new.i %in% c("TRUE", "FALSE"))
          new.i <- as.logical(new.i)
        else
          new.i <- as.numeric(new.i)
      }
      setOption(args[i], new.i)
    }
  }
  else {
    if (!quietly) {
      if (substring(setting, 1, 1) == "s")
        setting <- paste0("S", substring(setting, 2))
      cat(paste0("\n** ", setting, " already in used (R package meta). **\n\n"))
    } 
  }
}
