#' Counts the changes in a series of time bins
#' 
#' Given a vector of dates for a series of time bins and another for the times when a character change occurred will return the total number of changes in each bin.
#' 
#' Calculates the total number of evolutionary changes in a series of time bins. This is intended as an internal function for rate calculations, but could be used for other purposes (e.g., counting any point events in a series of time bins).
#' 
#' @param taxon.name Text.
#' @param taxon.number Text.
#' @param occurrences Text.
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
#' @export PaleobioDBTaxonScraper
PaleobioDBTaxonScraper <- function(taxon.name, taxon.number = NULL, occurrences = FALSE) {
    
    # Load gdata library for trim function:
    require(gdata)
    
    # Set as NA to start:
    X <- NA
    
    # If no taxon number is given:
    if(is.null(taxon.number)) {
        
        # Abbreviated name:
        if(length(grep("\\.", taxon.name)) > 0) { #taxon.name <- "A. patagonicus"
            
            taxa.found <- readLines(paste("http://paleobiodb.org/data1.1/taxa/auto.json?taxon_name=", gsub(" ", "", taxon.name), "&limit=10", sep = ""), warn = TRUE)
            
            taxa.found <- strsplit(gsub("}", "", gsub("},", "", gsub("\\{", "", taxa.found[3:(length(taxa.found) - 2)]))), ",")
            
        }
        
        # Single name:
        if(length(grep(" ", taxon.name)) == 0 && length(grep("\\.", taxon.name)) == 0) { #taxon.name <- "Tyrannosauridae"
            
            taxa.found <- readLines(paste("http://paleobiodb.org/data1.1/taxa/single.json?name=", taxon.name, "&limit=10", sep = ""), warn = TRUE)
            
            taxa.found <- strsplit(gsub("}", "", gsub("},", "", gsub("\\{", "", taxa.found[3:(length(taxa.found) - 2)]))), ",")
            
        }
        
        # Binomial name:
        if(length(grep(" ", taxon.name)) > 0 && length(grep("\\.", taxon.name)) == 0) { #taxon.name <- "Tyrannosaurus rex"
            
            taxa.found <- readLines(paste("http://paleobiodb.org/data1.1/taxa/list.json?name=", strsplit(taxon.name, " ")[[1]][1], "&rel=children", sep = ""), warn = TRUE)
            
            taxa.found <- strsplit(gsub("}", "", gsub("},", "", gsub("\\{", "", taxa.found[3:(length(taxa.found) - 2)]))), ",")
            
            taxa.found <- taxa.found[grep(taxon.name, taxa.found)]
            
            if(length(taxa.found) == 0) taxa.found <- list(character(0))
            
        }
        
        taxon.numbers <- unlist(strsplit(unlist(taxa.found)[grep("\"oid\"", unlist(taxa.found))], ":"))[unlist(strsplit(unlist(taxa.found)[grep("\"oid\"", unlist(taxa.found))], ":")) != "\"oid\""]
        
        taxon.names <- gsub("\"", "", unlist(strsplit(unlist(taxa.found)[grep("\"nam\"", unlist(taxa.found))], ":"))[unlist(strsplit(unlist(taxa.found)[grep("\"nam\"", unlist(taxa.found))], ":")) != "\"nam\""])
        
        taxa.found <- cbind(taxon.numbers, taxon.names)
        
        # Stop and warn if no taxa have that name:
        if(nrow(taxa.found) == 0) stop("Taxon name not present in database (although may be there as a junior synonym).")
        
        # Stop and warn if multiple taxa have that name:
        if(nrow(taxa.found) > 1) stop("Multiple taxa have that name. Check the database manually to verify the one you want and try again using the taxon number.\n")
        
        # Can now get taxon number:
        if(nrow(taxa.found) == 1) taxon.number <- as.numeric(taxa.found[1, 1])
        
    }
    
    # Query database using taxon number:
    while(is.na(X[1])) X <- strsplit(readLines(paste("http://paleobiodb.org/data1.1/taxa/single.tsv?id=", gsub(" ", "%20", taxon.number), "&show=attr,app,size", sep = ""), warn = TRUE), "\t")
    
    # As long as the taxon was found:
    if(unlist(X)[1] != "THIS REQUEST DID NOT GENERATE ANY OUTPUT RECORDS") {
        
        # Make output into matrix:
        X <- matrix(unlist(X), nrow = 2, byrow = T)
        
        # Add first row as headers:
        colnames(X) <- X[1, ]
        
        # Remove first row:
        X <- X[-1, ]
        
        # If occurrences are requested (and there is at least one):
        if(occurrences && X["firstapp_ea"] != "") {
            
            # Set as NA to start:
            Y <- NA
            
            # Query database using taxon number:
            while(is.na(Y[1])) Y <- readLines(paste("http://paleobiodb.org/data1.1/occs/list.tsv?taxon_id=", gsub(" ", "%20", taxon.number), "&show=coords,attr,ident,loc,paleoloc,prot,time,strat,lith,geo", sep=""))
            
            # Check that each value for Y is equal in length (number of tabs, i.e. columns):
            Y.lengths <- unlist(lapply(trim(strsplit(gsub("\t", "\t ", Y), "\t")), length))
            
            # If there zero-value cells:
            if(length(which(Y.lengths == 0)) > 0) {
                
                # Remove them:
                Y <- Y[-which(Y.lengths == 0)]
                
                # Update Y lengths:
                Y.lengths <- unlist(lapply(trim(strsplit(gsub("\t", "\t ", Y), "\t")), length))
                
            }
            
            # If they are unequal:
            if(length(unique(Y.lengths)) > 1) {
                
                # Until they are equal:
                while(length(unique(Y.lengths)) > 1) {
                    
                    # Find first (shortened) row:
                    first.problem.row <- min(which(Y.lengths < max(Y.lengths)))
                    
                    # Case if line is not last:
                    if(first.problem.row < length(Y)) {
                        
                        # Collapse with subsequent row:
                        Y[first.problem.row] <- paste(Y[first.problem.row:(first.problem.row + 1)], collapse="")
                        
                        # Remove subsequent row:
                        Y <- Y[-(first.problem.row + 1)]
                        
                        # Case if last line is too short:
                    } else {
                        
                        # Collapse with previous row:
                        Y[(first.problem.row - 1)] <- paste(Y[(first.problem.row - 1):first.problem.row], collapse="")
                        
                        # Remove row:
                        Y <- Y[-first.problem.row]
                        
                    }
                    
                    # Update lengths:
                    Y.lengths <- unlist(lapply(trim(strsplit(gsub("\t", "\t ", Y), "\t")), length))
                    
                }
                
            }
            
            # Make into matrix:
            Y <- matrix(unlist(trim(strsplit(gsub("\t", "\t ", Y), "\t"))), nrow=length(Y), byrow=T)
            
            # Add column names:
            colnames(Y) <- Y[1, ]
            
            # Remove column name row:
            Y <- Y[-1, ]
            
            # Make sure Y has not been turned into a vector:
            if(!is.matrix(Y)) Y <- matrix(Y, nrow=1, dimnames=c("", list(names(Y))))
            
            # If occurrences were not requested (or do not exist):
        } else {
            
            # Return NULL variable:
            Y <- NULL
            
        }
        
        # If taxon was not found:
    } else {
        
        # Return warning for X
        X <- "Taxon not present in database"
        
        # Set Y as NULL:
        Y <- NULL
        
    }
    
    # Compile output:
    out <- list(X, Y)
    
    # Add names to output:
    names(out) <- c("Taxon.data", "Occurrences")
    
    # Return output:
    return(out)
    
}
