#' Prunes an MRP matrix of redundant characters
#'
#' Prunes an MRP matrix of redundant characters (Pisani et al. 2002)
#'
#' @param MRPMatrix An MRP matrix in the format imported by \code{Claddis::ReadMorphNexus}.
#'
#' @return A modified matrix with any now redundant characters removed.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Pisani et al. 2002
#'
#' @examples
#'
#' # Nothing yet.
#'
#' @export PisaniMRPPrune
PisaniMRPPrune <- function(MRPMatrix) {
  
  # Check there are no step matrices:
  if(any(!is.null(MRPMatrix$Topper$StepMatrices))) stop("Function not meant for use with step matrices.")
  
  # Check is only a single block and stop and warn user if not:
  if(length(MRPMatrix) > 2) stop("MRP matrices should consist of a single block.")
  
  # Check is only zeroes and ones and stop and warn user if not:
  if(length(setdiff(unique(as.vector(MRPMatrix$Matrix_1$Matrix)), c("0", "1"))) > 0) stop("MRP matrix must consist on only zeroes and ones.")
  
  # Check for zero weight characters and stop and warn user if found:
  if(any(MRPMatrix$Matrix_1$Weights == 0)) stop("Function not intended for zero weight characters.")
  
  # Get any duplicated characters:
  DuplicatedCharacters <- which(duplicated(apply(MRPMatrix$Matrix_1$Matrix, 2, paste, collapse = "")))
  
  # Get any constant characters:
  ConstantCharacters <- which(apply(MRPMatrix$Matrix_1$Matrix, 2, function(x) length(unique(x))) < 2)
  
  # Join to make characters to prune:
  CharactersToPrune <- sort(unique(c(DuplicatedCharacters, ConstantCharacters)))
  
  # Special case of all characters being pruned:
  if(length(CharactersToPrune) == ncol(MRPMatrix$Matrix_1$Matrix)) {
    
    # Make matrix empty (zero by zero):
    MRPMatrix$Matrix_1$Matrix <- matrix(nrow = 0, ncol = 0)
    
    # Make ordering an empty character vector:
    MRPMatrix$Matrix_1$Ordering <- vector(mode = "character")
    
    # Make all otehr values an empty numeric vector:
    MRPMatrix$Matrix_1$Weights <- MRPMatrix$Matrix_1$MinVals <- MRPMatrix$Matrix_1$MaxVals <- vector(mode = "numeric")
  
  # As longas characters to prune is fewer than total number of characters:
  } else {
    
    # If characters to prune are found then prune them:
    if(length(CharactersToPrune) > 0) MRPMatrix <- MatrixPruner(MRPMatrix, characters2prune = CharactersToPrune)
  
  }
  
  # Return MRP matrix:
  return(MRPMatrix)

}
