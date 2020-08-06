#' Prunes an MRP matrix of any redundant characters
#'
#' @description
#'
#' Prunes an MRP matrix of any redundant (duplicate or autapomorphic) characters following Pisani et al. (2002).
#'
#' @param
#'
#' MRPMatrix An MRP matrix in the format imported by \code{Claddis::read_nexus_matrix}. Must contain only 0 or 1 states.
#'
#' @details
#'
#' Matrix representations must be handled carefully when taxa (rows) or characters (columns) are being pruned or duplicated, as they are by the \link{Metatree} function. Specifically, as noted by Pisani et al. (2002; their Figure 1) if a row gets pruned it will often create duplicate characters that, if unchecked, will inadvertently upweight some relationships incorrectly.
#'
#' To illustrate this, consider the following scenario as proposed by Pisani et al. (2002). A matrix representation for a tree with six tips (labelled A-F) is encoded with four characters (nodes or bipartitions, labelled 1-4):
#'
#' \preformatted{  1 2 3 4
#' A 0 0 0 0
#' B 0 0 0 1
#' C 0 0 1 1
#' D 0 1 1 1
#' E 1 1 1 1
#' F 1 1 1 1}
#'
#' Then, following removal of taxon D (e.g., because taxonomic reconciliation suggests it is a nomen dubium) the matrix will now look like this:
#'
#' \preformatted{  1 2 3 4
#' A 0 0 0 0
#' B 0 0 0 1
#' C 0 0 1 1
#' E 1 1 1 1
#' F 1 1 1 1}
#'
#' We now see that characters 1 and 2 are identical, effectively upweighting this relationship in the matrix representation.
#'
#' This can be countered by removing the offending character (column), to give us the fairer representation:
#'
#' \preformatted{  1 3 4
#' A 0 0 0
#' B 0 0 1
#' C 0 1 1
#' E 1 1 1
#' F 1 1 1}
#'
#' This function performs a check and correction for such duplicate characters, but also removes those that become autapomorphic in a similar way. For example, if we also removed taxon F:
#'
#' \preformatted{  1 3 4
#' A 0 0 0
#' B 0 0 1
#' C 0 1 1
#' E 1 1 1}
#'
#' Character 1 is now redundant, not because it is duplicated but because it carries no parsimony-informative information. NB: in practical terms this wouldn't effect the results, but removing characters like this can dramatically reduce data set size and hence minimise memory issues downstream. In this case, collapsing the data to just two characters:
#'
#' \preformatted{  3 4
#' A 0 0
#' B 0 1
#' C 1 1
#' E 1 1}
#'
#' Note that primarily this function is intended for internal use by \link{Metatree}, but may have some value to a user as a stand alone function.
#'
#' @return
#'
#' A modified MRP matrix with redundant (duplicate or autapomorphic) characters removed.
#'
#' @author
#'
#' Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Pisani, D., Yates, A. M., Langer, M. C. and Benton, M. J., 2002. A genus-level supertree of the Dinosauria. \emph{Proceedings of the Royal Society of London B}, \bold{269}, 915-921.
#'
#' @seealso
#'
#' For another means by which MRP matrices require collapsing see \link{CollapseDuplicateTaxonMRP}.
#'
#' @examples
#'
#' # Build the Pisani et al. (2002) example matrix:
#' PisaniExample <- Claddis::build_cladistic_matrix(matrix(as.character(c(0, 0,
#'   0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)),
#'   nrow = 6, dimnames = list(LETTERS[1:6], c()), byrow = TRUE))
#'
#' # Confirm it matches the Pisani example by inspection:
#' PisaniExample$matrix_1$matrix
#'
#' # Remove taxon D
#' PisaniExample <- Claddis::prune_cladistic_matrix(PisaniExample, taxa2prune = c("D"))
#'
#' # Confirm that first two characters are redundant by inspection:
#' PisaniExample$matrix_1$matrix
#'
#' # Apply pruning alogirhtm to data:
#' PisaniExample <- PisaniMRPPrune(PisaniExample)
#'
#' # Confirm that first two characters are now merged by inspection:
#' PisaniExample$matrix_1$matrix
#'
#' @export PisaniMRPPrune
PisaniMRPPrune <- function(MRPMatrix) {
  
  # Check there are no step matrices:
  if(any(!is.null(MRPMatrix$Topper$StepMatrices))) stop("Function not meant for use with step matrices.")
  
  # Check is only a single block and stop and warn user if not:
  if(length(MRPMatrix) > 2) stop("MRP matrices should consist of a single block.")
  
  # Check is only zeroes and ones and stop and warn user if not:
  if(length(setdiff(unique(as.vector(MRPMatrix$matrix_1$matrix)), c("0", "1"))) > 0) stop("MRP matrix must consist on only zeroes and ones.")
  
  # Check for zero weight characters and stop and warn user if found:
  if(any(MRPMatrix$matrix_1$character_weights == 0)) stop("Function not intended for zero weight characters.")
  
  # Get any duplicated characters:
  DuplicatedCharacters <- which(duplicated(apply(MRPMatrix$matrix_1$matrix, 2, paste, collapse = "")))
  
  # Get any constant characters:
  ConstantCharacters <- which(apply(MRPMatrix$matrix_1$matrix, 2, function(x) length(unique(x))) < 2)
  
  # Get any autapomorphic characters:
  AutapomorphicCharacters <- which(apply(MRPMatrix$matrix_1$matrix, 2, function(x) sum(as.numeric(x))) == 1)
  
  # Join to make characters to prune:
  CharactersToPrune <- sort(unique(c(DuplicatedCharacters, ConstantCharacters, AutapomorphicCharacters)))
  
  # Special case of all characters being pruned:
  if(length(CharactersToPrune) == ncol(MRPMatrix$matrix_1$matrix)) {
    
    # Make matrix empty (zero by zero):
    MRPMatrix$matrix_1$matrix <- matrix(nrow = 0, ncol = 0)
    
    # Make ordering an empty character vector:
    MRPMatrix$matrix_1$ordering <- vector(mode = "character")
    
    # Make all otehr values an empty numeric vector:
    MRPMatrix$matrix_1$character_weights <- MRPMatrix$matrix_1$minimum_values <- MRPMatrix$matrix_1$maximum_values <- vector(mode = "numeric")
  
  # As longas characters to prune is fewer than total number of characters:
  } else {
    
    # If characters to prune are found then prune them:
    if(length(CharactersToPrune) > 0) MRPMatrix <- Claddis::prune_cladistic_matrix(MRPMatrix, characters2prune = CharactersToPrune)
  
  }
  
  # Return MRP matrix:
  return(MRPMatrix)

}
