#' Gets file list for given clade
#'
#' @description
#'
#' Gets available file list for a given clade from the web site graemetlloyd.com.
#'
#' @param HTMLendings Vector of strings giving matrix ending for desired HTML.
#'
#' @details
#'
#' A list of published cladistic analyses is available from \href{http://www.graemetlloyd.com/matr.html}{graemetlloyd.com}, but not all of these have available data sets. This function serves as a tool to grab a list of the file names available for a specific clade, for example \href{http://www.graemetlloyd.com/matricht.html}{ichthyopterygians} (see example code below).
#'
#' The intended purpose of this function is to isolate filenames for a given clade so only these are handed as input to the \code{Metatree} function instead of using the full set, many of which might be irrelevant for a given target clade.
#'
#' @return
#'
#' A vector of available file names (without appendix).
#'
#' @author
#'
#' Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Get all ichthyosaur file names from graemetlloyd.com:
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
