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
#' @export MRPCollapse
# MRP collapse function (removes redundant characters a la Pisani et al. 2002):
MRPCollapse <- function(matrix) {
    
    matrix <- as.matrix(matrix)
    taxon.names <- rownames(matrix)
    transposed.matrix <- t(matrix)
    unq.rows <- vector(mode="numeric")
    for(i in length(transposed.matrix[, 1]):1) {
        unq.rows <- unique(c(unq.rows, paste(transposed.matrix[i, ], collapse="")))
    }
    unq.rows <- sort(unique(unq.rows), decreasing=TRUE)
    if(length(unq.rows) > 0) {
        for(i in length(unq.rows):1) {
            if(sum(as.numeric(strsplit(unq.rows[i], "")[[1]])) < 2) {
                unq.rows <- unq.rows[-i]
            }
        }
        if(length(unq.rows) > 0) {
            transposed.matrix <- vector(mode="character")
            for(i in 1:length(unq.rows)) {
                transposed.matrix <- rbind(transposed.matrix, strsplit(unq.rows[i], "")[[1]])
            }
            if(nrow(as.matrix(transposed.matrix)) > 1) {
                out <- t(transposed.matrix)
                rownames(out) <- taxon.names
            } else {
                out <- as.matrix(transposed.matrix)
                colnames(out) <- taxon.names
                out <- t(out)
            }
        } else {
            out <- "There are no informative characters after collapsing"
        }
    } else {
        out <- "There are no informative characters after collapsing"
    }
    return(out)
}
