#' Palaeobiology Database Taxa Querier
#'
#' Given a Paleobiology Database taxon number returns basic information on that taxon.
#'
#' Uses the Paleobiology Database (\code{paleobiodb.org}) API to query known taxon numbers and returns information on their validity, name, and rank. Intended for use in building dynamic taxonomic resolutions when constructing metatree matrices (see Lloyd et al. 2016).
#'
#' @param taxon_nos A vector of Paleobiology database taxon numbers to retrieve from the database.
#' @param taxon_names A vector of taxon names to search for in the database (default left to NULL).
#' @param original Whether or not to return the original (TRUE) or resolved version (FALSE).
#' @param interval The beginning and ending geologic periods if only wanting taxa from a specified time window (default is NULL).
#' @param extant What to do with extant taxa, one of: "only" (only return extant taxa), "exclude" (exclude extant taxa), or "include" (make no distinction, the default).
#' @param stopfororphans Whether or not to stop with an Error message for taxa with no parent.
#' @param breaker Size of breaker to use if querying a large number of taxa (reduces load on database of individual queries; default is 100).
#'
#' @return An eight-column matrix detailing the original taxon number (if relevant), the valid (resolved) taxon number, the taxon name, the taxon rank (Paleobiology Database rank number), the taxon number of the parent of each taxon, the taxon validity (if relevant; returns NA if already valid), the accepted taxon number (if relevant), and the accepted taxon name (if relevant).
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references Lloyd, G. T., Bapst, D. W., Friedman, M. and Davis, K. E., 2016. Probabilistic divergence time estimation without branch lengths: dating the origins of dinosaurs, avian flight, and crown birds. Biology Letters, 12, 20160609.
#'
#' @examples
#'
#' # Taxon query for Allosaurus fragilis:
#' PaleobiologyDBTaxaQuerier(taxon_nos = "52962")
#'
#' @export PaleobiologyDBTaxaQuerier
PaleobiologyDBTaxaQuerier <- function(taxon_nos, taxon_names = NULL, original = TRUE, interval = NULL, extant = "include", stopfororphans = TRUE, breaker = 100) {
  
  # List of geologic periods (no stage sor other intervals for now) in geolgic order:
  GeologicPeriodsInOrder <- c("Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian", "Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene")
  
  # Check extant option is a valid chocie and stop and warn user f not
  if(length(setdiff(extant, c("exclude", "include", "only"))) > 0) stop("extant must be one of \"exclude\", \"include\", or \"only\".")
  
  # If extant option is include set text to empty string (will be default in API):
  if(extant == "include") extantoption <- ""
  
  # If extant option is only set text to extant=yes:
  if(extant == "only") extantoption <- "&extant=yes"
  
  # If extant option is exclude set text to extant=no:
  if(extant == "exclude") extantoption <- "&extant=no"
  
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
    ResolvedHTTPStrings <- lapply(NumbersToQuery, function(x) ifelse(original, paste("https://paleobiodb.org/data1.2/taxa/list.json?id=", paste(paste("var:", x, sep = ""), collapse = ","), "&show=parent", sep = ""), paste("https://paleobiodb.org/data1.2/taxa/list.json?id=", paste(paste("txn:", x, sep = ""), collapse = ","), "&show=parent", sep = "")))
    
  }
  
  # If querying taxon names:
  if(!is.null(taxon_names)) {
    
    # Build list of names to query:
    NamesToQuery <- lapply(NumberChunker(N = length(taxon_names), breaker = breaker), function(x) taxon_names[x])
    
    # Build HTTP string(s):
    ResolvedHTTPStrings <- lapply(NamesToQuery, function(x) paste("https://paleobiodb.org/data1.2/taxa/list.json?name=", paste(gsub(" ", "%20", gsub("_", " ", gdata::trim(x))), collapse = ","), "&show=parent", sep = ""))
    
  }
  
  # Add extant option to query:
  ResolvedHTTPStrings <- paste(ResolvedHTTPStrings, extantoption, sep = "")
  
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
    
    # retrieve taxon validity, if known:
    TaxonValidity <- ParameterExtraction(jsonstring, parameterstring = "\"tdf\":")
    
    # retrieve taxon validity, if known:
    AcceptedNumber <- ParameterExtraction(jsonstring, parameterstring = "\"acc\":")
    
    # retrieve taxon validity, if known:
    AcceptedName <- ParameterExtraction(jsonstring, parameterstring = "\"acn\":")
    
    # Compile output:
    output <- cbind(OriginalTaxonNo, ResolvedTaxonNo, TaxonName, TaxonRank, ParentTaxonNo, TaxonValidity, AcceptedNumber, AcceptedName)
    
    # Add names:
    colnames(output) <- c("OriginalTaxonNo", "ResolvedTaxonNo", "TaxonName", "TaxonRank", "ParentTaxonNo", "TaxonValidity", "AcceptedNumber", "AcceptedName")
    
    # Return output:
    return(output)
    
  }))
  
  # If querying based on taxon number(s):
  if(is.null(taxon_names)) {
    
    # Extract taxon numbers found in database:
    TaxonNumbers <- unlist(lapply(lapply(lapply(apply(cbind(gsub("var:", "", Output[, "OriginalTaxonNo"]), gsub("txn:", "", Output[, "ResolvedTaxonNo"])), 1, as.list), unlist), function(x) x[!is.na(x)]), '[', 1))
    
    # Sort output by taxon numbers in the order they were supplied:
    Output <- Output[match(taxon_nos, TaxonNumbers), , drop = FALSE]

  }
  
  # If querying based on taxon name(s):
  if(!is.null(taxon_names)) {
    
    # Sort output by taxon names in order they were supplied:
    Output <- Output[match(gsub("_", " ", taxon_names), Output[, "TaxonName" ]), , drop = FALSE]
    
  }
  
  # If stopping for orphans check for them (excluding life) and stop if found:
  if(stopfororphans) if(any(sort(as.numeric(is.na(Output[, "ParentTaxonNo"])) - as.numeric(Output[, "TaxonName"] == "Life") == 1))) {
    
    # Make vector of orphan tax(a) names:
    OrphanTaxa <- which(as.numeric(is.na(Output[, "ParentTaxonNo"])) - as.numeric(Output[, "TaxonName"] == "Life") == 1)
    
    # Stop and warn user of orphans:
    if(length(OrphanTaxa) > 0) stop(paste("The following taxa are orphans: ", paste(Output[OrphanTaxa, "TaxonName"], collapse = ", "), sep = ""))
    
  }
  
  # Return output:
  return(Output)
  
}

# Example of orphan error:
#PaleobiologyDBTaxonQuerier(taxon_nos = as.character(370001:370100), taxon_names = NULL, original = TRUE, stopfororphans = TRUE)

#Examples of different types of taxa:
#PaleobiologyDBTaxonQuerier(taxon_nos = c("64336", "52962", "38856", "52773", "285777"), taxon_names = NULL, original = TRUE, stopfororphans = TRUE)
