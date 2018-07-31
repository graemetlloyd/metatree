#' Palaeobiology Database Children Finder
#'
#' Given a Paleobiology Database taxon number returns basic information on all species-level children.
#'
#' Uses the Paleobiology Database (\code{paleobiodb.org}) API to query a known taxon number (or name) and returns information on the validity, name, and rank of all its species-level children. Intended for use in building dynamic taxonomic resolutions when building metatree matrices (see Lloyd et al. 2016).
#'
#' @param taxon_no The Paleobiology database taxon number.
#' @param taxon_name A taxon name to search for in the database (default left to NULL); overrides taxon_no if used.
#' @param original Whether or not to return the original (TRUE) or resolved version (FALSE).
#'
#' @return A ten-column matrix detailing the original taxon number (if relevant), the valid (resolved) taxon number, the taxon name, the taxon rank (Paleobiology Database rank number), the taxon number of the parent of this taxon, the taxon validity (if relevant; returns NA if already valid), the accepted taxon number (if relevant), the accepted taxon name (if relevant) of all species-level children found, the attribution of the original name as currently entered in the database, and whether ("1") or not ("0") the species is extant.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references Lloyd, G. T., Bapst, D. W., Friedman, M. and Davis, K. E., 2016. Probabilistic divergence time estimation without branch lengths: dating the origins of dinosaurs, avian flight, and crown birds. Biology Letters, 12, 20160609.
#'
#' @examples
#'
#' # Taxon query for Loxommatinae:
#' PaleobiologyDBChildrenFinder("339413")
#'
#' @export PaleobiologyDBChildrenFinder
PaleobiologyDBChildrenFinder <- function(taxon_no, taxon_name = NULL, original = TRUE) {
  
  # Shows resolved taxon name for a given id (regular taxa only):
  resolvedhttpstring <- ifelse(original, paste("https://paleobiodb.org/data1.2/taxa/list.json?id=var:", taxon_no, "&rel=all_children", sep = ""), paste("https://paleobiodb.org/data1.2/taxa/list.json?id=txn:", taxon_no, "&rel=all_children&pres=regular&show=attr,is_extant", sep = ""))
  
  # Overwwrite taxon number query if using the taxon name instead (regular taxa only):
  if(!is.null(taxon_name)) resolvedhttpstring <- paste("https://paleobiodb.org/data1.2/taxa/list.json?name=", gsub(" ", "%20", gdata::trim(taxon_name)), "&rel=all_children&pres=regular&show=attr,is_extant", sep = "")
  
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
  
  # Trim to just taxon rows:
  resolvedjson <- resolvedjson[which(lapply(strsplit(resolvedjson, ""), '[', 1) == "{")][-1]
  
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
    #jsonstring <- jsonstring[(grep("\\[", jsonstring) + 1):(grep("\\]", jsonstring) - 1)]
    
    # Retrieve original taxon number (should be same as input!), if found:
    OriginalTaxonNo <- ParameterExtraction(jsonstring, parameterstring = "\"vid\":")
    
    # Retrieve resolved taxon number, if found:
    ResolvedTaxonNo <- ParameterExtraction(jsonstring, parameterstring = "\"oid\":")
    
    # Retrieve taxon name, if found:
    TaxonName <- ParameterExtraction(jsonstring, parameterstring = "\"nam\":")
    
    # Retrieve taxon rank if found:
    TaxonRank <- ParameterExtraction(jsonstring, parameterstring = "\"rnk\":")
    
    # Retrieve parent taxon number, if found:
    ParentTaxonNo <- ParameterExtraction(jsonstring, parameterstring = "\"par\":")
    
    # Retrieve taxon validity, if known:
    TaxonValidity <- ParameterExtraction(jsonstring, parameterstring = "\"tdf\":")
    
    # Retrieve taxon validity, if known:
    AcceptedNumber <- ParameterExtraction(jsonstring, parameterstring = "\"acc\":")
    
    # Retrieve taxon validity, if known:
    AcceptedName <- ParameterExtraction(jsonstring, parameterstring = "\"acn\":")
    
    # Retrieve taxon validity, if known:
    Attribution <- ParameterExtraction(jsonstring, parameterstring = "\"att\":")
    
    # Retrieve taxon validity, if known:
    Extant <- ParameterExtraction(jsonstring, parameterstring = "\"ext\":")

    # Compile output:
    output <- list(OriginalTaxonNo, ResolvedTaxonNo, TaxonName, TaxonRank, ParentTaxonNo, TaxonValidity, AcceptedNumber, AcceptedName, Attribution, Extant)
    
    # Add names:
    names(output) <- c("OriginalTaxonNo", "ResolvedTaxonNo", "TaxonName", "TaxonRank", "ParentTaxonNo", "TaxonValidity", "AcceptedNumber", "AcceptedName", "Attribution", "Extant")
    
    # Return output:
    return(output)
    
  }
  
  # Format output as list:
  outputlist <- lapply(as.list(resolvedjson), jsontotext)
  
  # Reformat list as matrix:
  outputmatrix <- matrix(unlist(outputlist), nrow = length(resolvedjson), byrow = TRUE, dimnames = list(c(), names(outputlist[[1]])))
  
  # Collapse to just species:
  outputmatrix <- outputmatrix[which(outputmatrix[, "TaxonRank"] == "3"), , drop = FALSE]
  
  # List of types of resolution that require deletion:
  deletes <- c("corrected to", "misspelling of", "objective synonym of", "obsolete variant of", "recombined as", "replaced by", "subjective synonym of", "nomen dubium", "nomen vanum", "nomen nudum", "nomen oblitum", "invalid subgroup of")
  
  # If all taxa are deletes:
  if(length(setdiff(unique(outputmatrix[, "TaxonValidity"]), deletes)) == 0) {
    
    # Make empty outputmatrix:
    outputmatrix <- outputmatrix[-(1:nrow(outputmatrix)), , drop = FALSE]
    
  # If at least some taxa are valid:
  } else {
    
    # If there are any issues with taxon validity:
    if(any(!is.na(outputmatrix[, "TaxonValidity"]))) {
      
      # Find any rows to delete:
      rowstodelete <- which(apply(apply(as.matrix(deletes), 1, '==', outputmatrix[, "TaxonValidity"]), 1, sum) == 1)
      
      # If found, remove them:
      if(length(rowstodelete) > 0) outputmatrix <- outputmatrix[-rowstodelete, , drop = FALSE]
      
    }
    
  }
  
  # Return output:
  return(outputmatrix)
  
}
