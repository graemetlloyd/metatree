#' Collapse MRP matrix
#'
#' Collapses an MRP matrix to unique characters
#'
#' Collapses an MRP matrix to just unique characters following the protocol laid out in Pisani et al. (2002).
#'
#' @param clad.matrix An MRP matrix in \link{ReadMorphNexus} format.
#'
#' @return An MRP matrix in \link{ReadMorphNexus} format where all characters are unique.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Pisani et al. (2002)
#'
#' @examples
#'
#' # Nothing yet.
#'
#' @export MRPCollapse
MRPCollapse <- function(clad.matrix) {
  
  # Get number of characters:
  NChars <- sum(unlist(lapply(lapply(clad.matrix[2:length(clad.matrix)], '[[', "Matrix"), ncol)))
  
  # Find any constant characters:
  constantcharacters <- which(unlist(lapply(lapply(split(clad.matrix$Matrix_1$Matrix, rep(1:ncol(clad.matrix$Matrix_1$Matrix), each = nrow(clad.matrix$Matrix_1$Matrix))), unique), length)) < 2)
  
  # Prune these if found (unless all characters are constant):
  if(length(constantcharacters) > 0 && length(constantcharacters) < NChars) clad.matrix <- MatrixPruner(clad.matrix, characters2prune = constantcharacters)
  
  # If all characters are constant prune without replacing matrices with NULL:
  if(length(constantcharacters) == NChars) clad.matrix[2:length(clad.matrix)] <- lapply(clad.matrix[2:length(clad.matrix)], function(x) {x$Matrix <- x$Matrix[, -(1:ncol(x)), drop = FALSE]; return(x)})
  
  # Only continue if matrix still has at least one column:
  if(ncol(clad.matrix$Matrix_1$Matrix) > 0) {
    
    # Find any zero-weight characters:
    zeroweightcharacters <- which(clad.matrix$Matrix_1$Weights == 0)
    
    # Prune these if found:
    if(length(zeroweightcharacters) > 0) clad.matrix <- MatrixPruner(clad.matrix, characters2prune = zeroweightcharacters)
    
    # Only continue if matrix still has at least one column:
    if(ncol(clad.matrix$Matrix_1$Matrix) > 0) {
      
      # Check for step matrices and stop if found:
      if(!is.null(clad.matrix$Topper$StepMatrices)) stop("Function not designed to work with step matrices.")
      
      # Check for non-binary characters and stop if found:
      if(length(unique(as.vector(clad.matrix$Matrix_1$Matrix))) != 2) stop("Function not designed for non-binary data.")
      
      # Store taxon names:
      taxon.names <- rownames(clad.matrix$Matrix_1$Matrix)
      
      # Collapse cladistic matrix to just unique characters:
      clad.matrix$Matrix_1$Matrix <- matrix(unlist(strsplit(sort(apply(clad.matrix$Matrix_1$Matrix, 2, paste, collapse = ""))[!duplicated(sort(apply(clad.matrix$Matrix_1$Matrix, 2, paste, collapse = "")))], "")), nrow = length(taxon.names), dimnames = list(taxon.names, c()))
      
      # Overwrite ordering to new length:
      clad.matrix$Matrix_1$Ordering <- rep("unord", ncol(clad.matrix$Matrix_1$Matrix))
      
      # Overwrite weights to new length:
      clad.matrix$Matrix_1$Weights <- rep(1, ncol(clad.matrix$Matrix_1$Matrix))
      
      # Overwrite max vals to new length:
      clad.matrix$Matrix_1$MaxVals <- rep(1, ncol(clad.matrix$Matrix_1$Matrix))
      
      # Overwrite min vals to new length:
      clad.matrix$Matrix_1$MinVals <- rep(0, ncol(clad.matrix$Matrix_1$Matrix))
      
    }
    
  }
  
  # Return collapsed matrix:
  return(clad.matrix)
  
}
