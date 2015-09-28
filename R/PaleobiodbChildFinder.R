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
#' @export PaleobioDBChildFinder
PaleobioDBChildFinder <- function(taxon.name, taxon.number=NULL) {
    
    # Load gdata library for trim function:
    require(gdata)
    
    # Set as NA to start:
    X <- NA
    
    # If no taxon number is given:
    if(is.null(taxon.number)) {
        
        # First check database for all taxa matching that name:
        taxa.found <- strsplit(readLines(paste("http://paleobiodb.org/data1.1/taxa/auto.tsv?name=", taxon.name, "&limit=10", sep=""), warn=T), "\t")
        
        # Make output into a matrix:
        taxa.found <- matrix(unlist(taxa.found), ncol=6, byrow=T)
        
        # Place column names row as column names:
        colnames(taxa.found) <- taxa.found[1, ]
        
        # Collapse to just homonyms:
        taxa.found <- taxa.found[taxa.found[, "taxon_name"] == taxon.name, ]
        
        # Ensure result is stll a matrix:
        if(!is.matrix(taxa.found)) taxa.found <- matrix(taxa.found, nrow=1, dimnames=c("", list(names(taxa.found))))
        
        # Stop and warn if multiple taxa have that name:
        if(nrow(taxa.found) > 1) stop("Multiple taxa have that name. Check the database manually to verify the one you want and try again using the taxon number.\n")
        
        # Query database using taxon name:
        while(is.na(X[1])) X <- trim(strsplit(gsub("\t", "\t ", readLines(paste("http://paleobiodb.org/data1.1/taxa/list.tsv?name=", gsub(" ", "%20", taxon.name), "&rel=all_children&status=valid&show=attr,app,size,phylo,nav,img,ent,entname,crmod", sep=""), warn=T)), "\t"))
        
        # If taxon number is given:
    } else {
        
        # Query database using taxon number:
        while(is.na(X[1])) X <- trim(strsplit(gsub("\t", "\t ", readLines(paste("http://paleobiodb.org/data1.1/taxa/list.tsv?id=", taxon.number, "&rel=all_children&status=valid&show=attr,app,size,phylo,nav,img,ent,entname,crmod", sep=""), warn=T)), "\t"))
        
    }
    
    # As long as the taxon was found:
    if(unlist(X)[1] != "THIS REQUEST DID NOT GENERATE ANY OUTPUT RECORDS") {
        
        # Make output into matrix:
        X <- matrix(unlist(X), nrow=length(X), byrow=T)
        
        # Add first row as headers:
        colnames(X) <- X[1, ]
        
        # Remove first row:
        X <- X[-1, ]
        
        # If taxon was not found:
    } else {
        
        # Return warning for X
        X <- "Taxon not present in database"
        
    }
    
    # Return output:
    return(X)
    
}
