#' Forest plot to display the result of a leave-one-out meta-analysis
#' 
#' @description
#' Draws a forest plot showing the results of omitting each study in turn,
#' using the grid graphics system.
#' 
#' @aliases forest.metainf
#' 
#' @param x An object of class \code{\link{metainf}}.
#' @param layout A character string specifying the layout of the
#'   forest plot (see \code{\link{forest.meta}}).
#' @param prediction A logical indicating whether prediction
#'   intervals should be printed.
#' @param overall A logical indicating whether overall results should be
#'   shown.
#' @param leftcols A character vector specifying columns to be printed on the
#'   left side of the forest plot or a logical value; see
#'   \code{\link{forest.meta}}.
#' @param rightcols A character vector specifying columns to be printed on the
#'   right side of the forest plot or a logical value; see
#'   \code{\link{forest.meta}}.
#' @param leftlabs A character vector specifying labels for columns on the left
#'   side of the forest plot (see Details).
#' @param rightlabs A character vector specifying labels for columns on the
#'   right side of the forest plot (see Details).
#' @param lab.NA A character string to label missing values.
#' @param col A named list specifying colours used in the forest
#'   plot; see \code{\link{forest-formatting}}.
#' @param print A named list specifying which information to print in the forest
#'   plot; see \code{\link{forest-formatting}}.
#' @param type A character string or vector specifying how to
#'   plot treatment effects and confidence intervals for meta-analysis results.
#' @param digits A named list specifying digits used for displayed results; a
#'   scalar value is also accepted. See \code{\link{forest-formatting}}.
#' @param details A logical specifying whether details on statistical
#'   methods should be printed.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in forest plots. If \code{backtransf = TRUE},
#'   results for \code{sm = "OR"} are presented as odds ratios rather
#'   than log odds ratios, for example.
#' @param big.mark A character used as thousands separator.
#' @param main A title for the forest plot.
#' @param filename File name to store forest plot.
#' @param path Path of the directory to save plot to: `path` and `filename` are
#'   combined to create the fully qualified file name. Defaults to the working
#'   directory.
#' @param device Device to use. Can either be a device function (e.g.,
#'   \code{\link{png}}), or one of \code{"eps"}, \code{"ps"},
#'   \code{"tex"} (pictex), \code{"pdf"}, \code{"jpeg"}, \code{"tiff"},
#'   \code{"png"}, \code{"bmp"}, \code{"svg"} or \code{"wmf"} (Windows only).
#'   By default, the device is guessed based on the `filename` extension.
#' @param device.args List with additional graphical parameters passed on
#'   to graphics function.
#' @param units Units of the `width` and `height` of the file containing the
#'   forest plot. One of \code{"in"} (for inches, default), \code{"cm"},
#'   \code{"mm"}, or \code{"px"} (for pixels), can be abbreviated.
#' @param height Height of the graphics file.
#' @param width Width of the graphics file.
#' @param dpi Plot resolution.
#' @param dev.off A logical to specify whether current graphics device
#'   should be shut down.
#' @param new A logical value indicating whether a new figure should
#'   be printed in an existing graphics window.
#' @param \dots Additional graphical arguments (passed on to
#'   \code{\link{forest.meta}}).
#' 
#' @details
#' Each row of the forest plot represents a meta-analysis after omitting one
#' study from the original meta-analysis. Results are shown in the order of
#' the omitted studies. Common-effect or random-effects estimates, prediction
#' intervals, heterogeneity statistics, and study-level information can be
#' selected with the corresponding arguments.
#'
#' A forest plot, also called a confidence interval plot, is drawn in the
#' active graphics window using the grid graphics system. To save the plot,
#' specify a file name with argument \code{filename}; the graphics device is
#' selected from the file extension by default (see \code{\link{forest.meta}}
#' for examples).
#' 
#' The named list \code{print} determines which result columns are shown in
#' the forest plot. The argument \code{layout} determines the default
#' placement of these columns on the left or right side. The columns and their
#' placement are determined internally and cannot be changed by the user.
#'
#' Arguments \code{leftcols} and \code{rightcols} can be used to specify
#' columns printed on the left or right side, respectively. Arguments
#' \code{leftlabs} and \code{rightlabs} can be used to provide headings for
#' these columns. Use \code{NA} for a column heading to retain its default
#' label.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{forest.metacum}}, \code{\link{forest.meta}},
#'   \code{\link{metainf}}, \code{\link{settings.meta}}
#' 
#' @keywords hplot
#' 
#' @examples
#' data(Fleiss1993bin)
#' ma <- metabin(d.asp, n.asp, d.plac, n.plac,
#'   data = Fleiss1993bin, studlab = study, sm = "RR", method = "I")
#' ma
#' metainf(ma)
#' metainf(ma, pooled = "random")
#' 
#' forest(metainf(ma))
#' forest(metainf(ma, pooled = "random"))
#' forest(metainf(ma, pooled = "random", prediction = TRUE))
#'
#' @method forest metainf
#' @export

forest.metainf <- function(x,
                           #
                           layout = gs("layout"),
                           #
                           prediction = x$prediction,
                           overall = x$overall,
                           #
                           leftcols = NULL, rightcols = NULL,
                           leftlabs = NULL, rightlabs = NULL,
                           #
                           lab.NA = gs("lab.NA"),
                           #
                           col = NULL,
                           print = NULL,
                           type = "square",
                           digits = NULL,
                           #
                           details = gs("forest.details"),
                           #
                           backtransf = x$backtransf,
                           big.mark = gs("big.mark"),
                           #
                           main = "Leave-One-Out Meta-Analysis",
                           #
                           # Output
                           #
                           filename = NULL,
                           path = NULL,
                           device = NULL, device.args = NULL,
                           units = "in",
                           height = NULL, width = NULL,
                           dpi = 300,
                           dev.off = NULL,
                           new = TRUE,
                           #
                           ...) {
  
  chkclass(x, "metainf")
  x <- updateversion(x)
  
  invisible(
    forest.metacum(x,
                   #
                   layout = layout,
                   #
                   prediction = prediction,
                   overall = overall,
                   #
                   leftcols = leftcols, rightcols = rightcols,
                   leftlabs = leftlabs, rightlabs = rightlabs,
                   #
                   lab.NA = lab.NA,
                   #
                   col = col,
                   print = print,
                   type = type,
                   digits = digits,
                   #
                   details = details,
                   #
                   backtransf = backtransf,
                   big.mark = big.mark,
                   #
                   main = main,
                   #
                   # Output
                   #
                   filename = filename,
                   path = path,
                   device = device, device.args = device.args,
                   units = units,
                   height = height, width = width,
                   dpi = dpi,
                   dev.off = dev.off,
                   new = new,
                   #
                   ...)
  )
}
