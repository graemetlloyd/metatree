#' Palaeobiology Database Taxon Querier
#' 
#' Given a Paleobiology Database taxon number returns basic information on that taxon.
#' 
#' Uses the Paleobiology Database (\code{paleobiodb.org}) API to query a known taxon number and returns information on its validity, name, and rank. Intended for use in building dynamic taxonomic resolutions when building metatree matrices (see Lloyd et al. 2016).
#' 
#' @param taxon_no The Paleobiology database taxon number.
#' @param taxon_name A taxon name to search for in the database (default left to NULL).
#' @param original Whether or not to return the original (TRUE) or resolved version (FALSE).
#'
#' @return A eight-item list detailing the original taxon number (if relevant), the valid (resolved) taxon number, the taxon name, the taxon rank (Paleobiology Database rank number), the taxon number of the parent of this taxon, the taxon validity (if relevant; returns NA if already valid), the accepted taxon number (if relevant), and the accepted taxon name (if relevant).
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references Lloyd, G. T., Bapst, D. W., Friedman, M. and Davis, K. E., 2016. Probabilistic divergence time estimation without branch lengths: dating the origins of dinosaurs, avian flight, and crown birds. Biology Letters, 12, 20160609.
#'
#' @examples
#' 
#' # Taxon query for Allosaurus fragilis:
#' PaleobiologyDBTaxonQuerier("52962")
#' 
#' @export PaleobiologyDBTaxonQuerier
PaleobiologyDBTaxonQuerier <- function(taxon_no, taxon_name = NULL, original = TRUE) {
    
    # Shows resolved taxon name for a given id:
    resolvedhttpstring <- ifelse(original, paste("https://paleobiodb.org/data1.2/taxa/single.json?id=var:", taxon_no, "&show=parent", sep = ""), paste("https://paleobiodb.org/data1.2/taxa/single.json?id=txn:", taxon_no, "&show=parent", sep = ""))
    
    # Overwwrite taxon number query if using the taxon name instead:
    if(!is.null(taxon_name)) resolvedhttpstring <- paste("https://paleobiodb.org/data1.2/taxa/single.json?name=", gsub(" ", "%20", trim(taxon_name)), "&show=parent", sep = "")
    
    # Set resolved json to NA (used later to check results are coming back from server):
    resolvedjson <- NA
    
    # Set start value for counter (used later to avoid infinite loop):
    counter <- 0
    
    # While server has not been reached (and querying a taxon number):
    while(is.na(resolvedjson[[1]][1]) && is.null(taxon_name)) {
        
        # Attempt to acquire resolved taxon string:
        try(resolvedjson <- readLines(resolvedhttpstring), silent = TRUE)
        
        # If server was not reached:
        if(is.na(resolvedjson[[1]][1])) {
            
            # Update counter to record how many attempts to reach server have been made:
            counter <- counter + 1
            
            # If repeatedly failing to get results stop trying:
            if(counter == 100) stop("Server not responding after 100 straight attempts")
            
            # Wait two seconds before next attempt (also avoids overloading server):
            Sys.sleep(2)
            
        }
        
    }
    
    # If using a taxon name to query:
    if(!is.null(taxon_name)) {
        
        # Ask server for data:
        try(resolvedjson <- readLines(resolvedhttpstring), silent = TRUE)
        
        # If no data returned tell user:
        if(length(resolvedjson) == 1) stop(paste("Could not find record for ", taxon_name, " in database.", sep = ""))
        
    }
    
    # Subfunction to return specific information from json data:
    jsontotext <- function(jsonstring) {
        
        # Subfunction to extract specific parameter from json string:
        ParameterExtraction <- function(jsonstring, parameterstring) {
            
            # Extract specific paramter:
            output <- ifelse(length(grep(parameterstring, jsonstring)) > 0, gsub("\"", "", strsplit(strsplit(jsonstring, parameterstring)[[1]][2], ",")[[1]][1]), NA)
            
            # Return output:
            return(output)
            
        }
        
        # Isolate record line:
        jsonstring <- jsonstring[(grep("\\[", jsonstring) + 1):(grep("\\]", jsonstring) - 1)]
        
        # Retrueve original taxon number (should be same as input!), if found:
        OriginalTaxonNo <- ParameterExtraction(jsonstring, parameterstring = "\"vid\":")
        
        # Retrieve resolved taxon number, if found:
        ResolvedTaxonNo <- ParameterExtraction(jsonstring, parameterstring = "\"oid\":")
        
        # Retrieve taxon name, if found:
        TaxonName <- ParameterExtraction(jsonstring, parameterstring = "\"nam\":")
        
        # Retrieve taxon rank if found:
        TaxonRank <- ParameterExtraction(jsonstring, parameterstring = "\"rnk\":")
        
        # Retrieve parent taxon number, if found:
        ParentTaxonNo <- ParameterExtraction(jsonstring, parameterstring = "\"par\":")
        
        # retrieve taxon validity, if known:
        TaxonValidity <- ParameterExtraction(jsonstring, parameterstring = "\"tdf\":")
        
        # retrieve taxon validity, if known:
        AcceptedNumber <- ParameterExtraction(jsonstring, parameterstring = "\"acc\":")

        # retrieve taxon validity, if known:
        AcceptedName <- ParameterExtraction(jsonstring, parameterstring = "\"acn\":")

        # Compile output:
        output <- list(OriginalTaxonNo, ResolvedTaxonNo, TaxonName, TaxonRank, ParentTaxonNo, TaxonValidity, AcceptedNumber, AcceptedName)
        
        # Add names:
        names(output) <- c("OriginalTaxonNo", "ResolvedTaxonNo", "TaxonName", "TaxonRank", "ParentTaxonNo", "TaxonValidity", "AcceptedNumber", "AcceptedName")
        
        # Return output:
        return(output)
        
    }
    
    # Get output:
    output <- jsontotext(resolvedjson)
    
    # Return output:
    return(output)
    
}

#taxon_no <- "64336" # Example of a nomen dubium: Claosaurus affinis
#taxon_no <- "52962" # Example of a valid taxon: Allosaurus fragilis
#taxon_no <- "38856" # Example of a junior synonym: Eoceratops
#taxon_no <- "52773" # Example of a senior synonym: Chasmosaurinae
#taxon_no <- "285777" # Example of a misspelled taxon name
