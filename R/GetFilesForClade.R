#' Gets file list for given clade
#'
#' Gets file list for given clade from HTML (graemetllouyd.com).
#'
#' @param HTMLendings Vector of strings giving matrix ending for desired HTML.
#'
#' @return A vector of file names intended for use with the \code{Metatree} function.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Get all ichthyosaur file names
#' # from graemetlloyd.com:
#' GetFilesForClade("matricht.html")
#'
#' @export GetFilesForClade
GetFilesForClade <- function(HTMLendings) {
  
  # Build full URL for each ending:
  URLs <- paste("http://www.graemetlloyd.com/", HTMLendings, sep = "")
  
  # Extract files to export:
  FilesToExport <- sort(unlist(lapply(as.list(URLs), function(x) {y <- readLines(x); unlist(lapply(strsplit(y[grep(".xml", y)], "href=\"xml/|.xml"), '[[', 2))})))
  
  # Return files to export:
  return(sort(FilesToExport))
  
}
