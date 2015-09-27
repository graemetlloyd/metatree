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
mrp.tree <- function(tree) {
    mrp <- matrix(0, ncol=(Nnode(tree)-1), nrow=Ntip(tree)) # MRP matrix
    rownames(mrp) <- tree$tip.label # Name rows as taxa
    colnames(mrp) <- as.character(c((Ntip(tree)+2):(Ntip(tree)+Nnode(tree)))) # Name columns as nodes
    for(j in c((Ntip(tree)+2):(Ntip(tree)+Nnode(tree)))) { # For each internal node...
        mrp[tree$tip.label[FindDescendants(j, tree)], as.character(j)] <- 1 # Fill MRP matrix
    }
    mrp <- mrp[sort(rownames(mrp)), ]
    return(mrp)
}
