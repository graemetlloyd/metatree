#' Embiggen a cladistic matrix
#' 
#' Given a cladistic matrix duplicates the block N times.
#' 
#' Given a cladistic matrix duplicates the block N times.
#' 
#' @param CladisticMatrix A cladistic matrix in the format imported by \code{Claddis::ReadMorphNexus}.
#' @param N The number of times to duplicated the data.
#'
#' @return A cladistic matrix in the format imported by \code{Claddis::ReadMorphNexus}.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Nothing yet.
#'
#' @export EmbiggenMatrix
EmbiggenMatrix <- function(CladisticMatrix, N) {
  
  CladisticMatrix[2:length(CladisticMatrix)] <- lapply(CladisticMatrix[2:length(CladisticMatrix)], function(x) {
    
    
    MatrixBlock <- x$Matrix
    
    for(i in 2:N) x$Matrix <- cbind(x$Matrix, MatrixBlock)
    
    x$Ordering <- rep(x$Ordering, times = N)
    x$Weights <- rep(x$Weights, times = N)
    x$MinVals <- rep(x$MinVals, times = N)
    x$MaxVals <- rep(x$MaxVals, times = N)
    
    return(x)
    
  })
  
  return(CladisticMatrix)

}
