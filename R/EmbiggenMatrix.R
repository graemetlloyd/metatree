#' Embiggen a cladistic matrix
#'
#' @description
#'
#' Given a cladistic matrix, will duplicate each block of characters N times.
#'
#' @param CladisticMatrix A cladistic matrix in the format imported by \code{Claddis::ReadMorphNexus}.
#' @param N The number of times to duplicate the data.
#'
#' @details
#'
#' Given a cladistic matrix in the format imported by \code{Claddis::ReadMorphNexus} will duplicate each block N times. Thus if a block contains 10 characters and N is 5 the resulting block will have 50 characters total.
#'
#' This function exists for internal use in the \link{Metatree} function as a means of upweighting MRP data beyond the limits of TNT (maximum weight 1000) by duplicating blocks. However, it is also made available for individual use here.
#'
#' @return
#'
#' A cladistic matrix in the format imported by \code{Claddis::ReadMorphNexus}.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Build an example matrix of six species and four characters:
#' ExampleMatrix <- Claddis::MakeMorphMatrix(matrix(as.character(c(0, 0,
#'   0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)),
#'   nrow = 6, dimnames = list(LETTERS[1:6], c()), byrow = TRUE))
#'
#' # Inspect starting matrix:
#' ExampleMatrix$Matrix_1$Matrix
#'
#' # Duplicate block five times (20 characters):
#' ExampleMatrix <- EmbiggenMatrix(ExampleMatrix, 5)
#'
#' # Inspect matrix to check it now contains 20 characters:
#' ExampleMatrix$Matrix_1$Matrix
#'
#' @export EmbiggenMatrix
EmbiggenMatrix <- function(CladisticMatrix, N) {
  
  # TO DO:
  # - Add data checks
  
  # Duplicate main block N times:
  CladisticMatrix[2:length(CladisticMatrix)] <- lapply(CladisticMatrix[2:length(CladisticMatrix)], function(x) {
    
    # Isolate matrix block:
    MatrixBlock <- x$Matrix
    
    # Keep adding block to end of matrix N times:
    for(i in 2:N) x$Matrix <- cbind(x$Matrix, MatrixBlock)
    
    # Now update ordering:
    x$Ordering <- rep(x$Ordering, times = N)
    
    # Now update weights:
    x$Weights <- rep(x$Weights, times = N)
    
    # Now update minimum values:
    x$MinVals <- rep(x$MinVals, times = N)
    
    # Now update maximum values:
    x$MaxVals <- rep(x$MaxVals, times = N)
    
    # Return duplicated block:
    return(x)
    
  })
  
  # Return embiggened matrix:
  return(CladisticMatrix)
  
}
