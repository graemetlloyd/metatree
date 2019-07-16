#' Gives position of the kth item in a list
#' 
#' Gives position of the kth item in a list of vectors.
#' 
#' @param Position The kth value for which the position is desired.
#' @param ListDimensions The dimensions (vector length) of each part of the list in order.
#'
#' @return The position, in tersm of both list number ($ListNumber) and list position ($ListPosition), of the kth value.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Find the position of the 13th value:.
#' ListPosition2ListCoordinates(13, c(1,3,2,6,1,3))
#'
#' @export ListPosition2ListCoordinates
ListPosition2ListCoordinates <- function(Position, ListDimensions) {
  
  # Add zero to beginning of list dimensions (makes following simpler):
  ListDimensions <- c(0, ListDimensions)
  
  # Build positions matrix (from-to numbers):
  PositionMatrix <- cbind((cumsum(ListDimensions) + 1)[1:(length(ListDimensions) - 1)], cumsum(ListDimensions)[2:length(ListDimensions)])
  
  # Build list of indices (1:N):
  ListIndices <- lapply(apply(PositionMatrix, 1, as.list), function(x) {y <- unlist(x); y[1]:y[2]})
  
  # Find the list number for the desired position:
  ListNumber <- which(unlist(lapply(ListIndices, function(x) sum(x == Position))) == 1)
  
  # Find the list position for the desired position:
  ListPosition <- which(ListIndices[[ListNumber]] == Position)
  
  # Build output:
  Output <- list(ListNumber, ListPosition)
  
  # Add naems to output:
  names(Output) <- c("ListNumber", "ListPosition")
  
  # Return output:
  return(Output)
  
}
