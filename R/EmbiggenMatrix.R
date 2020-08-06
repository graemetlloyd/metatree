#' Embiggen a cladistic matrix
#'
#' @description
#'
#' Given a cladistic matrix, will duplicate each block of characters N times.
#'
#' @param CladisticMatrix A cladistic matrix in the format imported by \code{Claddis::read_nexus_matrix}.
#' @param N The number of times to duplicate the data.
#'
#' @details
#'
#' Given a cladistic matrix in the format imported by \code{Claddis::read_nexus_matrix} will duplicate each block N times. Thus if a block contains 10 characters and N is 5 the resulting block will have 50 characters total.
#'
#' This function exists for internal use in the \link{Metatree} function as a means of upweighting MRP data beyond the limits of TNT (maximum weight 1000) by duplicating blocks. However, it is also made available for individual use here.
#'
#' @return
#'
#' A cladistic matrix in the format imported by \code{Claddis::read_nexus_matrix}.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Build an example matrix of six species and four characters:
#' ExampleMatrix <- Claddis::build_cladistic_matrix(matrix(as.character(c(0, 0,
#'   0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)),
#'   nrow = 6, dimnames = list(LETTERS[1:6], c()), byrow = TRUE))
#'
#' # Inspect starting matrix:
#' ExampleMatrix$matrix_1$matrix
#'
#' # Duplicate block five times (20 characters):
#' ExampleMatrix <- EmbiggenMatrix(ExampleMatrix, 5)
#'
#' # Inspect matrix to check it now contains 20 characters:
#' ExampleMatrix$matrix_1$matrix
#'
#' @export EmbiggenMatrix
EmbiggenMatrix <- function(CladisticMatrix, N) {
  
  # TO DO:
  # - Add data checks
  
  # Duplicate main block N times:
  CladisticMatrix[2:length(CladisticMatrix)] <- lapply(CladisticMatrix[2:length(CladisticMatrix)], function(x) {
    
    # Isolate matrix block:
    MatrixBlock <- x$matrix
    
    # Keep adding block to end of matrix N times:
    for(i in 2:N) x$matrix <- cbind(x$matrix, MatrixBlock)
    
    # Now update ordering:
    x$ordering <- rep(x$ordering, times = N)
    
    # Now update weights:
    x$character_weights <- rep(x$character_weights, times = N)
    
    # Now update minimum values:
    x$minimum_values <- rep(x$minimum_values, times = N)
    
    # Now update maximum values:
    x$maximum_values <- rep(x$maximum_values, times = N)
    
    # Return duplicated block:
    return(x)
    
  })
  
  # Return embiggened matrix:
  return(CladisticMatrix)
  
}
