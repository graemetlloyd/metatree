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
  
  # Find any constant characters:
  constantcharacters <- which(unlist(lapply(lapply(split(clad.matrix$matrix, rep(1:ncol(clad.matrix$matrix), each = nrow(clad.matrix$matrix))), unique), length)) < 2)
  
  # Prune these if found:
  if(length(constantcharacters) > 0) clad.matrix <- MatrixPruner(clad.matrix, characters2prune = constantcharacters)
  
  # Only continue if matrix still has at least one column:
  if(ncol(clad.matrix$matrix) > 0) {
    
    # Find any zero-weight characters:
    zeroweightcharacters <- which(clad.matrix$weights == 0)
    
    # Prune these if found:
    if(length(zeroweightcharacters) > 0) clad.matrix <- MatrixPruner(clad.matrix, characters2prune = zeroweightcharacters)
    
    # Only continue if matrix still has at least one column:
    if(ncol(clad.matrix$matrix) > 0) {
      
      # Check for step matrices and stop if found:
      if(!is.null(clad.matrix$step.matrices[[1]][1])) stop("Function not designed to work with step matrices.")
      
      # Check for non-binary characters and stop if found:
      if(length(unique(as.vector(clad.matrix$matrix))) != 2) stop("Function not designed for non-binary data.")
      
      # Store taxon names:
      taxon.names <- rownames(clad.matrix$matrix)
      
      # Collapse cladistic matrix to just unique characters:
      clad.matrix$matrix <- matrix(unlist(strsplit(sort(apply(clad.matrix$matrix, 2, paste, collapse = ""))[!duplicated(sort(apply(clad.matrix$matrix, 2, paste, collapse = "")))], "")), nrow = length(taxon.names), dimnames = list(taxon.names, c()))
      
      # Overwrite ordering to new length:
      clad.matrix$ordering <- rep("unord", ncol(clad.matrix$matrix))
      
      # Overwrite weights to new length:
      clad.matrix$weights <- rep(1, ncol(clad.matrix$matrix))
      
      # Overwrite max vals to new length:
      clad.matrix$max.vals <- rep(1, ncol(clad.matrix$matrix))
      
      # Overwrite min vals to new length:
      clad.matrix$min.vals <- rep(0, ncol(clad.matrix$matrix))
      
    }
    
  }
  
  # Return collapsed matrix:
  return(clad.matrix)
  
}
