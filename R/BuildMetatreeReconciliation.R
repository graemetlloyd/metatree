#' Build metatree taxonomic reconcilaiotn
#'
#' @description
#'
#' Given a list of one or more names builds metatree taxonomic reconcilaiton line.
#'
#' @param Taxa A list of taxa exactly as written in the Paleobiology Database, except with underscores within names instead of spaces.
#' @param SortNames Logical indicating whether names should also be sorted (TRUE) or not (FALSE; the default).
#'
#' @details
#'
#' Really intended as a tool for the author, this helps automate the process of building XML lines for multiple-taxon reconciliations.
#'
#' @return
#'
#' Nothing is returned. Instead the requisite line is sent to the screen ready for copy-pasting into an XML file.
#'
#' @author
#'
#' Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Lloyd, G. T., Bapst, D. W., Friedman, M. and Davis, K. E., 2016. Probabilistic divergence time estimation without branch lengths: dating the origins of dinosaurs, avian flight, and crown birds. \emph{Biology Letters}, \bold{12}, 20160609.
#'
#' @examples
#' 
#' # Generate a set of two example trees:
#' Taxa <- c("Elephas_maximus", "Loxodonta_africana", "Loxodonta_cyclotis")
#'
#' # Convert to MRP and show just matrix:
#' BuildMetatreeReconciliation(Taxa)
#'
#' @export BuildMetatreeReconciliation
BuildMetatreeReconciliation <- function(Taxa, SortNames = FALSE) {
  
  # If sorting taxa then sort names:
  if(SortNames) Taxa <- sort(Taxa)
  
  # Output XML test for taxon reconciliaiton:
  cat(paste("<List recon_name=\"", paste(gsub(" ", "_", Taxa), collapse = ","), "\" recon_no=\"", paste(unlist(lapply(apply(PaleobiologyDBTaxaQuerier(taxon_nos = "1", taxon_names = Taxa), 1, list), function(x) {x <- unlist(x); gsub("txn:|var:", "", x[which(!is.na(x[1:2]))[1]])})), collapse = ";"), "\">", sep = ""))
  
}

