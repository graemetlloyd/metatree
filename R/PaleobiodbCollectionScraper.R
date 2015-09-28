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
#' @export PaleobioDBCollectionScraper
PaleobioDBCollectionScraper <- function(col.id) {
    
    # Grab the data:
    X <- strsplit(readLines(paste("http://paleobiodb.org/data1.1/colls/single.tsv?id=", col.id, "&show=bin,attr,ref,loc,paleoloc,prot,time,strat,stratext,lith,lithext,geo,rem,ent,entname,crmod", sep=""), warn=T), "\t")
    
    # If data breaks ovr lines:
    if(length(X) > 2) {
        
        # While line breaks exist:
        while(length(X) > 2) {
            
            # Stitch break:
            X[[2]][length(X[[2]])] <- paste(X[[2]][length(X[[2]])], X[[3]][1])
            
            # Place on single line:
            if(length(X[[3]]) > 1) X[[2]] <- c(X[[2]][1:length(X[[2]])], X[[3]][2:length(X[[3]])])
            
            # Remove now redundant line:
            X <- X[-3]
            
        }
        
    }
    
    # If no data is found:
    if(unlist(X)[1] == "THIS REQUEST DID NOT GENERATE ANY OUTPUT RECORDS") {
        
        # Warn as output:
        X <- "Collection not found"
        
        # If data is found:
    } else {
        
        # Make into matrix:
        X <- matrix(unlist(X), ncol=length(X[[1]]), byrow=T)
        
        # Add first row as headers:
        colnames(X) <- X[1, ]
        
        # Remove first row:
        X <- X[-1, ]
        
        # Turn back into matrix:
        X <- matrix(X, nrow=1, dimnames=c("", list(names(X))))
        
    }
    
    # Return output:
    return(X)
    
}
