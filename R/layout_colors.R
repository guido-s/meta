layout_colors <- function(layout, classes = NULL) {
  bmj <- layout == "BMJ"
  revman5 <- layout == "RevMan5"
  jama <- layout == "JAMA"
  #
  # Colour schemes for layouts:
  # - color1 - vertical line for common effect or random effects model
  # - color2 - diamond for meta-analysis results
  # - color3 - outer lines of diamonds
  # - color4 - square or circle colour
  # - color5 - outer lines of squares or circles
  # - color6 - line color (axes, reference line, header lines)
  # - color7 - prediction interval
  # - color8 - outer lines of prediction intervals
  # - color9 - subgroups
  # - color10 - color within squares or circles
  #
  if (bmj) {
    color1 <- "#6b58a6"
    color2 <- "#6b58a6"
    color3 <- "white"
    color4 <- "#6b58a6"
    color5 <- "white"
    color6 <- "#a7a9ac"
  }
  else if (jama) {
    color1 <- "black"
    color2 <- "lightblue"
    color3 <- "black"
    color4 <- "darkblue"
    color5 <- "darkblue"
    color6 <- "black"
  }
  else if (revman5) {
    color1 <- "black"
    color2 <- "black"
    color3 <- "black"
    #
    if (!is.null(classes)) {
      metacont <- "metacont" %in% classes
      metamean <- "metamean" %in% classes
      metabin <- "metabin" %in% classes
      #
      if (metacont | metamean)
        color4 <- color5 <- "green"
      else if (metabin)
        color4 <- color5 <- "blue"
      else
        color4 <- color5 <- "red"
    }
    else {
      color4 <- color5 <- "gray"
    }
    #
    color6 <- "black"
  }
  else {
    color1 <- "black"
    color2 <- "gray"
    color3 <- "black"
    color4 <- "gray"
    color5 <- "gray"
    color6 <- "black"
  }
  #
  color7 <- "red"
  color8 <- "black"
  color9 <- "black"
  color10 <- "white"
  color11 <- "royalblue"
  
  res <- c(color1, color2, color3, color4, color5, color6,
           color7, color8, color9, color10, color11)
  res
}
