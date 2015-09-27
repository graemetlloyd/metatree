#' Counts the changes in a series of time bins
#' 
#' Given a vector of dates for a series of time bins and another for the times when a character change occurred will return the total number of changes in each bin.
#' 
#' Calculates the total number of evolutionary changes in a series of time bins. This is intended as an internal function for rate calculations, but could be used for other purposes (e.g., counting any point events in a series of time bins).
#' 
#' @param change.times A vector of ages in millions of years at which character changes are hypothesised to have occurred.
#' @param time.bins A vector of ages in millions of years of time bin boundaries in old-to-young order.
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
#' @export ChangesInBins
tntmrp.to.nexusmrp <- function(file)
{
    # Read in raw file:
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE) # Read in NEXUS file
    start <- grep("ROOT", X)+1
    end <- grep(";", X[start:length(X)])[1]+start-2
    MATRIX <- X[start:end]
    while(length(grep("  ", MATRIX)) > 0) MATRIX <- gsub("  ", " ", MATRIX)
    characters <- names <- vector(mode="character")
    for(i in 1:length(MATRIX)) {
        names[i] <- strsplit(MATRIX[i], " ")[[1]][1]
        characters[i] <- strsplit(MATRIX[i], " ")[[1]][2]
    }
    char.block <- matrix(nrow=length(names), ncol=nchar(characters[1]))
    for(i in 1:length(MATRIX)) char.block[i, ] <- strsplit(characters[i], "")[[1]]
    char.block <- t(char.block)
    characters <- vector(mode="character")
    for(i in 1:length(char.block[, 1])) characters[i] <- paste(char.block[i, ], collapse="")
    characters <- sort(characters)
    characters <- rle(characters)
    weights <- characters$lengths
    characters <- characters$values
    char.block <- char.block[1:length(characters), ]
    for(i in 1:length(characters)) char.block[i, ] <- strsplit(characters[i], "")[[1]]
    char.block <- t(char.block)
    rownames(char.block) <- names
    # Make into clad.matrix format:
    header <- ""
    ordering <- rep("unord", length(char.block[1, ]))
    max.vals <- rep(1, length(char.block[1, ]))
    min.vals <- rep(0, length(char.block[1, ]))
    result <- list(header, char.block, ordering, weights, max.vals, min.vals)
    names(result) <- c("header", "matrix", "ordering", "weights", "max.vals", "min.vals")
    return(result)
}
