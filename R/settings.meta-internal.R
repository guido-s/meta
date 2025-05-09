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
  ##
  invisible(NULL)
}


specificSettings <- function(args, new, setting, quietly = FALSE) {
  isnull.old <- as.vector(unlist(lapply(.settings[args], is.null)))
  ischar.old <- as.vector(unlist(lapply(.settings[args], is.character)))
  old <- rep("character", length(isnull.old))
  old[!isnull.old] <- as.vector(unlist(.settings[args]))
  ##
  isnull.new <- as.vector(unlist(lapply(new, is.null)))
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
      cat(paste0("\n** ", setting, " already in use (R package meta). **\n\n"))
    } 
  }
}


setcharacter <- function(argname, args, set = NULL, length = 1,
                         NULL.ok = FALSE, ignore.other = FALSE,
                         logical.ok = FALSE) {
  id <- argid(names(args), argname)
  ##
  if (!is.na(id)) {
    val <- args[[id]]
    ##
    if (NULL.ok & is.null(val)) {
      setOption(argname, val)
      return(invisible(NULL))
    }
    ##
    if (logical.ok & is.logical(val)) {
      setOption(argname, val)
      return(invisible(NULL))
    }
    ##
    if (!is.character(val) & ignore.other)
      return(invisible(id))
    ##
    if (!is.null(set))
      val <- setchar(val, set, name = argname)
    else
      chkchar(val, length = length, name = argname)
    ##
    setOption(argname, val)
  }
  ##
  invisible(id)
}


setcolor <- function(argname, args) {
  id <- argid(names(args), argname)
  ##
  if (!is.na(id)) {
    val <- args[[id]]
    chkcolor(val, name = argname)
    setOption(argname, val)
  }
  ##
  invisible(id)
}


setlevel <- function(argname, args) {
  id <- argid(names(args), argname)
  ##
  if (!is.na(id)) {
    val <- args[[id]]
    chklevel(val, name = argname)
    setOption(argname, val)
  }
  ##
  invisible(id)
}


setlogical <- function(argname, args, NULL.ok = FALSE,
                       ignore.other = FALSE) {
  id <- argid(names(args), argname)
  ##
  if (!is.na(id)) {
    val <- args[[id]]
    ##
    if (NULL.ok & is.null(val)) {
      setOption(argname, val)
      return(invisible(NULL))
    }
    ##
    if (!is.logical(val) & ignore.other)
      return(invisible(id))
    ##
    chklogical(val, name = argname)
    setOption(argname, val)
  }
  ##
  invisible(id)
}


setnumeric <- function(argname, args, NULL.ok = FALSE) {
  id <- argid(names(args), argname)
  ##
  if (!is.na(id)) {
    val <- args[[id]]
    ##
    if (NULL.ok & is.null(val)) {
      setOption(argname, val)
      return(invisible(NULL))
    }
    ##
    chknumeric(val, min = 0, length = 1, name = argname)
    setOption(argname, val)
  }
  ##
  invisible(id)
}
