#' Counts the changes in a series of time bins
#' 
#' Given a vector of dates for a series of time bins and another for the times when a character change occurred will return the total number of changes in each bin.
#' 
#' Calculates the total number of evolutionary changes in a series of time bins. This is intended as an internal function for rate calculations, but could be used for other purposes (e.g., counting any point events in a series of time bins).
#' 
#' @param filename Text.
#'
#' @return A vector giving the number of changes for each time bin. Names indicate the maximum and minimum (bottom and top) values for each time bin.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#' 
#' # Create a random dataset of 100 changes:
#' change.times <- runif(100, 0, 100)
#' 
#' # Create time bins:
#' time.bins <- seq(100, 0, length.out=11)
#' 
#' # Get N changes for each bin:
#' ChangesInBins(change.times, time.bins)
#' 
#' @export Trees2MPTsAndStrict
Trees2MPTsAndStrict <- function(filename) {
    
  # Read in tnt treefile
  X <- scan(file = filename, what = "", sep = "\n", quiet = TRUE)
  
  # Get just the trees
  X <- X[grep("\\(", X)]

  # Remove spaces
  X <- gsub(" ", "", X)

  # Get strict
  strict <- X[length(X)]
    
  # Get mpts
  mpts <- X[1:(length(X) - 1)]
    
  # Make output variable
  result <- list(mpts = mpts, strict = strict)
    
  # Return
  return(result)
  
}
