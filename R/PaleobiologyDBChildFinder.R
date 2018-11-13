#' Palaeobiology Database Child Finder
#'
#' Given a Paleobiology Database taxon number returns basic information on all species-level children.
#'
#' Uses the Paleobiology Database (\code{paleobiodb.org}) API to query a known taxon number (or name) and returns information on the validity, name, and rank of all its species-level children. Intended for use in building dynamic taxonomic resolutions when building metatree matrices (see Lloyd et al. 2016).
#'
#' @param taxon_nos The Paleobiology database taxon number.
#' @param taxon_names A taxon name to search for in the database (default left to NULL); overrides taxon_nos if used.
#' @param original Whether or not to return the original (TRUE) or resolved version (FALSE).
#' @param interval The beginning and ending geologic periods if only wanting taxa from a specified time window (default is NULL).
#' @param validonly Whether or not to only retunr valid taxa (TRUE) or all taxa (FALSE).
#' @param returnrank Whether or not to only return taxa of a specific rank (e.g., "3" for species, "5" for genera). See Paleobiology Database API for more infomation.
#' @param breaker Size of breaker to use if querying a large number of taxa (reduces load on database of individual queries).
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
#' PaleobiologyDBChildFinder("339413", returnrank = "3")
#'
#' @export PaleobiologyDBChildFinder
PaleobiologyDBChildFinder <- function(taxon_nos, taxon_names = NULL, original = TRUE, validonly = TRUE, returnrank = NULL, breaker = 100) {
  
  # TO DO: ALLOW EXTANT FILTERING AS WELL AS TIME
  
  # List of geologic periods (no stage sor other intervals for now) in geolgic order:
  GeologicPeriodsInOrder <- c("Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian", "Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene")

  # Subfunction to break N numbers into breaker-sized blocks:
  NumberChunker <- function(N, breaker) {
    
    # Get total numebr of chunks:
    NChunks <- ceiling(N / breaker)
    
    # Make initial list of numbers:
    ListOfNumbers <- rep(list(1:breaker), NChunks)
    
    # Update last element of list (which can be less than breaker) into the correct size (if required):
    if((N %% breaker) > 0) ListOfNumbers[[length(ListOfNumbers)]] <- 1:(N %% breaker)
    
    # Add cumulative breaker value to all numbers (if required):
    if(N > breaker) ListOfNumbers <- mapply('+', ListOfNumbers, c(0, cumsum(rep(breaker, NChunks - 1))), SIMPLIFY = FALSE)
    
    # Return list of numbers:
    return(ListOfNumbers)
    
  }
  
  # If querying taxon numbers:
  if(is.null(taxon_names)) {
    
    # Build list of numbers to query:
    NumbersToQuery <- lapply(NumberChunker(N = length(taxon_nos), breaker = breaker), function(x) taxon_nos[x])
    
    # Build HTTP string(s):
    ResolvedHTTPStrings <- lapply(NumbersToQuery, function(x) ifelse(original, paste("https://paleobiodb.org/data1.2/taxa/list.json?id=", paste(paste("var:", x, sep = ""), collapse = ","), "&rel=all_children&pres=regular&show=attr", sep = ""), paste("https://paleobiodb.org/data1.2/taxa/list.json?id=", paste(paste("txn:", x, sep = ""), collapse = ","), "&rel=all_children&pres=regular&show=attr", sep = "")))
    
  }
  
  # If querying taxon names:
  if(!is.null(taxon_names)) {
    
    # Build list of names to query:
    NamesToQuery <- lapply(NumberChunker(N = length(taxon_names), breaker = breaker), function(x) taxon_names[x])
    
    # Build HTTP string(s):
    ResolvedHTTPStrings <- lapply(NamesToQuery, function(x) paste("https://paleobiodb.org/data1.2/taxa/list.json?name=", paste(gsub(" ", "%20", gsub("_", " ", gdata::trim(x))), collapse = ","), "&rel=all_children&pres=regular&show=attr", sep = ""))
    
  }
  
  # If there are intervals supplied:
  if(!is.null(interval)) {
    
    # Check interval inputs are valid:
    if(length(setdiff(interval, GeologicPeriodsInOrder)) > 0) stop(paste("The following are not geologic periods (or are missspelled:", paste(setdiff(interval, GeologicPeriodsInOrder), collapse = ", ")))
    
    # If only one interval duplicate it so the next bit works:
    if(length(interval) == 1) interval <- c(interval, interval)
    
    # Check there aren't more than two periods being specified and stop and warn if found:
    if(length(interval) > 2) stop("Only supply beginning and end periods (i.e., two values) for interval.")
    
    # Get beginning adn ending matches from geologic periods:
    IntervalMatchesInOrder <- sort(match(interval, GeologicPeriodsInOrder))
    
    # Vuild interval string:
    IntervalStringToAdd <- paste("&interval=", paste(unique(GeologicPeriodsInOrder[IntervalMatchesInOrder[1]:IntervalMatchesInOrder[2]]), collapse =","), sep = "")
    
    # Add interval to resolved HTTP strings:
    ResolvedHTTPStrings <- lapply(ResolvedHTTPStrings, function(x) paste(x, IntervalStringToAdd, sep = ""))
    
  }

  # Get resolved json strings for each chunk:
  ResolvedJSON <- lapply(ResolvedHTTPStrings, function(x) {
    
    # Set resolved json to NA (used later to check results are coming back from server):
    resolvedjson <- NA
    
    # Set start value for counter (used later to avoid infinite loop):
    counter <- 0
    
    # While server has not been reached (and querying a taxon number):
    while(is.na(resolvedjson[[1]][1])) {
      
      # Attempt to acquire resolved taxon string:
      try(resolvedjson <- readLines(x), silent = TRUE)
      
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
    # Return resolvedjson string:
    return(resolvedjson)
    
  })
  
  # Extract data from json data:
  Output <- do.call(rbind, lapply(ResolvedJSON, function(x) {
    
    # Subfunction to extract specific parameter from json string:
    ParameterExtraction <- function(jsonstring, parameterstring) {
      
      # Extract specific parameter:
      output <- unlist(lapply(as.list(jsonstring), function(x) ifelse(length(grep(parameterstring, x)) > 0, gsub("\"", "", strsplit(strsplit(x, parameterstring)[[1]][2], ",")[[1]][1]), NA)))
      
      # Return output:
      return(output)
      
    }
    
    # If querying taxon numbers:
    if(is.null(taxon_names)) {
      
      # Find any taxon numbers not in database:
      UnknownTaxonHits <- grep("Unknown taxon", x)
      
      # If found stop and warn user:
      if(length(UnknownTaxonHits) > 0) stop(paste("The following taxon numbers were not found in the database: ", paste(unlist(lapply(strsplit(x[UnknownTaxonHits], split = "Unknown taxon '|'\""), '[', 2)), collapse = ", "), sep = ""))
      
    }
    
    # If querying taxon names:
    if(!is.null(taxon_names)) {
      
      # Find any taxon names not in database:
      UnknownTaxonHits <- grep("The name '", x)
      
      # If found stop and warn user:
      if(length(UnknownTaxonHits) > 0) stop(paste("The following taxon names were not found in the database: ", paste(unlist(lapply(strsplit(x[UnknownTaxonHits], split = "The name '|' did not match"), '[', 2)), collapse = ", "), sep = ""))
      
    }
    
    # Isolate record line:
    jsonstring <- x[(grep("\\[", x) + 1):(grep("\\]", x) - 1)]
    
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
    output <- cbind(OriginalTaxonNo, ResolvedTaxonNo, TaxonName, TaxonRank, ParentTaxonNo, TaxonValidity, AcceptedNumber, AcceptedName, Attribution, Extant)
    
    # Add names:
    colnames(output) <- c("OriginalTaxonNo", "ResolvedTaxonNo", "TaxonName", "TaxonRank", "ParentTaxonNo", "TaxonValidity", "AcceptedNumber", "AcceptedName", "Attribution", "Extant")
    
    # Return output:
    return(output)
    
  }))
  
  # If only requesting valid taxa remove all invalids from output:
  if(validonly) Output <- Output[is.na(Output[, "TaxonValidity"]), , drop = FALSE]
  
  # If requesting a specific rank of the output then filter by that rank:
  if(!is.null(returnrank)) Output <- Output[Output[, "TaxonRank"] == returnrank, , drop = FALSE]
  
  # Return output:
  return(Output)
  
}
