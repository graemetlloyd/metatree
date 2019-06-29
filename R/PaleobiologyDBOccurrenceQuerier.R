#' Palaeobiology Database Occurrence Querier
#'
#' Given a set of Paleobiology Database taxon number(s) returns occurrence information for those tax(a).
#'
#' Uses the Paleobiology Database (\code{paleobiodb.org}) API to query known taxon numbers and returns information on their occurrence as fossils.
#'
#' @param taxon_nos A vector of Paleobiology database taxon number(s) to retrieve from the database.
#' @param original Whether or not to return the original (TRUE) or resolved version (FALSE) of names.
#' @param breaker Size of breaker to use if querying a large number of taxa (reduces load on database of individual queries; default is 100).
#'
#' @return A multi-column matrix with rows as occurrences.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Occurrence query for Allosaurus fragilis:
#' PaleobiologyDBOccurrenceQuerier(taxon_nos = "52962")
#'
#' @export PaleobiologyDBOccurrenceQuerier
PaleobiologyDBOccurrenceQuerier <- function(taxon_nos, original = TRUE, breaker = 100) {
  
  # CHECK FOR EXTANT AND HENCE SET AGES AS ZERO?
  
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
  
  # Build list of numbers to query:
  NumbersToQuery <- lapply(NumberChunker(N = length(taxon_nos), breaker = breaker), function(x) taxon_nos[x])
  
  # Build HTTP string(s):
  ResolvedHTTPStrings <- lapply(NumbersToQuery, function(x) ifelse(original, paste("https://paleobiodb.org/data1.2/occs/list.json?taxon_id=", paste(paste("var:", x, sep = ""), collapse = ","), "&show=coords,paleoloc", sep = ""), paste("https://paleobiodb.org/data1.2/occs/list.json?taxon_id=", paste(paste("txn:", x, sep = ""), collapse = ","), "&show=coords,paleoloc", sep = "")))
  
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
    
    # Find any taxon numbers not in database:
    UnknownTaxonHits <- grep("Unknown taxon", x)
      
    # If found stop and warn user:
    if(length(UnknownTaxonHits) > 0) stop(paste("The following taxon numbers were not found in the database: ", paste(unlist(lapply(strsplit(x[UnknownTaxonHits], split = "Unknown taxon '|'\""), '[', 2)), collapse = ", "), sep = ""))
    
    # Isolate record line:
    jsonstring <- x[(grep("\\[", x) + 1):(grep("\\]", x) - 1)]
    
    # List of parameters to extract (check API documentation for meaning):
    Parameters <- c("cid", "idn", "tna", "oei", "eag", "lag", "lng", "lat", "pln", "pla")
    
    # Compile output:
    output <- do.call(cbind, lapply(as.list(Parameters), function(x) gsub("col:", "", ParameterExtraction(jsonstring, parameterstring = paste("\"", x, "\":", sep = "")))))
    
    # Add names:
    colnames(output) <- c("CollectionNo", "IdentifiedName", "TaxonName", "Age", "MaxMa", "MinMa", "Longitude", "Latitude", "PalaeoLongitude", "PalaeoLatitude")
    
    # Return output:
    return(output)
    
  }))
  
  # If there are any unidentified
  if(any(is.na(Output[, "IdentifiedName"]))) {
    
    # Store rows with NAs for identified name(s):
    Rows <- which(is.na(Output[, "IdentifiedName"]))
    
    # Overwrite NAs with taxon name:
    Output[Rows, "IdentifiedName"] <- Output[Rows, "TaxonName"]
    
  }
  
  # Return output:
  return(Output)
  
}

#taxon_nos <- c("251932", "64409", "66344", "159719", "192929", "182713", "85754", "66909", "174879", "172073", "66661", "55001", "55000", "56379", "133339", "54097", "65490", "96666", "56382", "162588", "370732", "370734", "372124", "55477", "67921", "91403", "378703", "347522", "57454", "243275", "68124", "370986", "242722", "153776", "243373", "119230", "133334", "131589", "57406", "142534", "66571", "53365", "57413", "57412", "56627", "62959", "243277", "378960", "56592", "65100", "56683", "57418", "56596", "64684", "68584", "56585", "327194", "56614", "56586", "66682", "68588", "55644", "322708", "54980", "335603")
#x <- PaleobiologyDBOccurrenceQuerier(taxon_nos)
#x <- matrix(cbind(rep(mean((as.numeric(x[, "MaxMa"]) + as.numeric(x[, "MinMa"])) / 2), nrow(x)), as.numeric(x[, "PalaeoLongitude"]), as.numeric(x[, "PalaeoLatitude"])), ncol = 3, dimnames = list(c(), c("recon_age", "paleolng", "paleolat")))
#maps <- getmap(ma = sort(unique(x[, "recon_age"])), model = "PALEOMAP", do.plot = FALSE)
#mapast(model = "PALEOMAP", data = as.data.frame(x), map = maps)
