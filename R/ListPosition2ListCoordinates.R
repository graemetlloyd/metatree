#' Gives the exact position of the kth item in a list
#' 
#' @description
#'
#' Gives the exact position of the kth item in a list comprised of vectors.
#' 
#' @param Position The kth value for which the position is desired.
#' @param ListDimensions The dimensions (vector length) of each part of the list in order.
#'
#' @details
#'
#' The R language has a great variety of ways that data can be stored, but sometimes it can be hard to isolate the required information. For example, imagine you store data in a list composed of vectors of length 1, 3, 2, 6, 1, and 3, but you just want to know where you will find the 13th value. This function will give "coordinates" for that value corresponding to the list number (5 in the example) and vector position in that list (1 in the example).
#'
#' @return
#'
#' The position, in terms of both list number ($ListNumber) and list position ($ListPosition), of the kth value.
#'
#' @author
#'
#' Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Find the position of the 13th value:.
#' ListPosition2ListCoordinates(13, c(1, 3, 2, 6, 1, 3))
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
  Output <- list(ListNumber = ListNumber, ListPosition = ListPosition)
  
  # Return output:
  return(Output)
  
}
