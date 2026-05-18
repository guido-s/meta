remove_html <- function(x) {
  x <- gsub("^<p>", "", x)
  x <- gsub("</p>$", "", x)
  x <- gsub("^Quote: ", "", x)
  x <- gsub("\"", '"', x)
  #
  x <- gsub("&gt;", ">", x)
  x <- gsub("&lt;", "<", x)
  x <- gsub("&nbsp;", " ", x)
  #
  x
}
