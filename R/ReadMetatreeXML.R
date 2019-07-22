#' Reads in a metatree-format XML file
#'
#' Reads in a metatree-format XML file (Lloyd et al. 2016).
#'
#' @param File Path to XML file to read.
#' @param Invisible Logical indicating whether to immediately show entire output (FALSE) or not (TRUE, the default).
#'
#' @return A nested list reflecting the nested XML tags of the input file.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Lloyd, G. T., Bapst, D. W., Friedman, M., and Davis, K. E., 2016. Probabilistic divergence time estimation without branch lengths: dating the origins of dinosaurs, avian flight and crown birds. Biology Letters, 12, 20160609.
#'
#' @examples
#'
#' # Nothing yet.
#'
#' @export ReadMetatreeXML
ReadMetatreeXML <- function(File, Invisible = TRUE) {
  
  # SEPARATE TEST THAT FILE BEGINS WITH ?xml .. (WEIRD TAG AS NOT CLOSE?)
  # IDEALLY WILL ALSO TEST FOR PRESENCE OF APPROPRIATE LINEBREAKS (DO NOT WANT DIFFERENT TAGS ON SAME LINE AT ANY POINT IN FILE)
  
  # OTHER CHECKS TO MAYBE ADD:
  #
  # 1. That taxon count matches number of taxa actually listed.
  # 2. That total number of characters matches external NEXUS file.
  # 3. Other checks from metatree master (spaces in names, numbers etc.).
  # 4. Author name sare formatted correctly.
  # 5. Rferences are formatted correctly overall.

  # Build simple tags database (can be edited in future) to set expected structure to be found in file - does not include generic <List> or <Type> subtags:
  TagsDatabase <- matrix(c("<SourceTree>", "", "<Source>", "<SourceTree>", "<Author>", "<Source>", "<Year>", "<Source>", "<Title>", "<Source>", "<Journal>", "<Source>", "<Volume>", "<Source>", "<Pages>", "<Source>", "<Booktitle>", "<Source>", "<Publisher>", "<Source>", "<City>", "<Source>", "<Editor>", "<Source>", "<Taxa>", "<SourceTree>", "<Characters>", "<SourceTree>", "<Molecular>", "<Characters>", "<Morphological>", "<Characters>", "<Behavioural>", "<Characters>", "<Other>", "<Characters>", "<Analysis>", "<SourceTree>", "<Notes>", "<SourceTree>", "<Filename>", "<SourceTree>", "<Parent>", "<SourceTree>", "<Sibling>", "<SourceTree>"), ncol = 2, byrow = TRUE, dimnames = list(c(), c("Tag", "Nesting")))
  
  # Build list of tags to ignore if foudn insde other tags (e.g., typeface tags like bold or emphasis):
  IgnoreTags <- c("<em>", "</em>")
  
  # Genearte vector of tag paths (will allow easy access of list later):
  TagPaths <- unlist(lapply(as.list(TagsDatabase[, "Tag"]), function(x) {CurrentTag <- x; CurrentTag <- c(CurrentTag, base::unname(TagsDatabase[TagsDatabase[, "Tag"] == CurrentTag, "Nesting"])); while(!any(CurrentTag == "")) CurrentTag <- c(CurrentTag, base::unname(TagsDatabase[TagsDatabase[, "Tag"] == CurrentTag[length(CurrentTag)], "Nesting"])); CurrentTag <- base::setdiff(CurrentTag, ""); paste(rev(base::gsub("<|>", "", CurrentTag)), collapse = "$")}))
  
  # Subfunction to convert two-column matrix of nestings (X nested in Y) to nested list:
  NestingMatrixToNestedList <- function(NestingMatrix) {
    
    # FOLLOWING WILL BREAK IF ORDER OF NESTING IS NOT LEAST TO MOST NESTED IN TAGS DATABASE MATRIX SO MAYBE SORT THIS? ALTHOUGHT HAt WOULD BREAK ORDER IN XML, SO MAYBE GET ORDER FROM XML ITSELF?

    # Generate iitial output as empty list:
    Output <- base::list()
    
    # Add in top level tag(s) to Output:
    Output[[base::gsub("<|>", "", NestingMatrix[NestingMatrix[, "Nesting"] == "", "Tag"])]] <- base::list()
    
    # Find any tags that have other tags nested inside them (effectively excludes <List> and <Type> which can duplicate or are generic):
    NestTags <- base::setdiff(base::unique(NestingMatrix[, "Nesting"]), "")
    
    # Build nested lists of tags ready to add to output:
    NestingTags <- base::lapply(base::as.list(NestTags), function(x) {ListNames <- base::gsub("<|>", "", NestingMatrix[NestingMatrix[, "Nesting"] == x, "Tag"]); y <- base::list(); for(i in ListNames) y[[i]] <- base::list(); z <- base::list(); z[[base::gsub("<|>", "", x)]] <- base::list(); z <- y; z})
    
    # For each nest tag:
    for(i in 1:length(NestTags)) {
      
      # Explicitly
      CurrentTag <- NestTags[i]
      
      # Add nxt nesting tag to current tag:
      CurrentTag <- c(CurrentTag, base::unname(NestingMatrix[NestingMatrix[, "Tag"] == CurrentTag, "Nesting"]))
      
      # Keep adding tags until have hit top-level:
      while(!any(CurrentTag == "")) CurrentTag <- c(CurrentTag, base::unname(NestingMatrix[NestingMatrix[, "Tag"] == CurrentTag[length(CurrentTag)], "Nesting"]))
      
      # Remove dead top-level "tag" (""):
      CurrentTag <- base::setdiff(CurrentTag, "")
      
      # Add nested tags to output:
      base::eval(base::parse(text = paste(paste(c("Output", rev(base::gsub("<|>", "", CurrentTag))), collapse = "$"), " <- NestingTags[[i]]", sep = "")))
      
    }
    
    # Return output:
    return(Output)
    
  }
  
  # Subfunction to grab data from inside a tag:
  RetrieveTagSupplement <- function(Tag, RawText) {
    
    # Make a raw match between the tag with space (implying extra data) and raw text:
    RawMatch <- RawText[base::grep(base::gsub(">", " ", Tag), RawText)]
    
    # As long as a raw match was made:
    if(length(RawMatch) > 0) {
      
      # Split raw match by spaces:
      RawMatch <- strsplit(RawMatch, " ")[[1]]
      
      # Find beginning of tag:
      TagBegins <- base::grep(base::gsub(">", "", Tag), RawMatch)
      
      # Find end of tag:
      TagEnds <- base::grep(">", RawMatch)[base::grep(">", RawMatch) > TagBegins][1]
      
      # Build new raw match from inner tag info:
      RawMatch <- unlist(lapply(as.list(RawMatch[(TagBegins + 1):TagEnds]), function(x) strsplit(x, ">")[[1]][1]))
      
      # Check tag contains equals symbol (required for split below):
      if(length(base::setdiff(1:length(RawMatch), base::grep("=", RawMatch))) > 0) stop(paste("Values inside ", Tag, " tag are missing equals symbol", sep = ""))
      
      # Build output matrix:
      Output <- matrix(unlist(strsplit(base::gsub("\"", "", RawMatch), "=")), ncol = 2, byrow = TRUE, dimnames = list(c(), c("Measure", "Value")))
      
    # If not inside tag data:
    } else {
      
      # Make output NULL:
      Output <- NULL
      
    }
    
    # Return output:
    return(Output)
    
  }
  
  # Build starting output (to fill as file read below):
  Output <- NestingMatrixToNestedList(TagsDatabase)
  
  # Read in raw text version of XML:
  RawText <- base::readLines(File)
  
  # Remove any empty lines:
  RawText <- RawText[base::nchar(RawText) > 0]
  
  # Trim any leading or trailing whitespace:
  RawText <- gdata::trim(RawText)
  
  # Find any missing tags:
  MissingTags <- sort(unlist(lapply(as.list(TagsDatabase[, "Tag"]), function(x) {Matches <- c(base::grep(base::gsub(">", "", x), RawText), base::grep(base::gsub(">", "/>", x), RawText)); ifelse(length(Matches) == 0, x, character(0))})))
  
  # If found stop and warn user:
  if(length(MissingTags) > 0) stop(paste("The following tag(s) are missing from the XML file: ", paste(MissingTags, collapse = ", "), ". Check for typographical errors or omissions and try again.", sep = ""))
  
  # Subfunction to actually read tag data from raw text:
  TagReader <- function(Tag, RawText) {
    
    # Formulate the unused tag:
    UnusedTag <- base::gsub(">", "/>", Tag)
    
    # If unused tag is found:
    if(length(base::grep(UnusedTag, RawText)) > 0) {
      
      # Make output NULL (as tag was not actually used):
      Output <- NULL
      
    # If no unused tag:
    } else {
      
      # Build opening tag:
      OpeningTag <- paste(base::gsub(">", " ", Tag), base::gsub(">", ">", Tag), sep ="|")
      
      # Build closing tag:
      ClosingTag <- base::gsub("<", "</", Tag)
      
      # Find position of opening tag:
      Start <- base::grep(OpeningTag, RawText)
      
      # Find position of closing tag:
      End <- base::grep(ClosingTag, RawText)
      
      # Check start tag is even found and stop and warn user if not:
      if(length(Start) == 0) stop(paste("Missing opening ", Tag, " tag.", sep = ""))
      
      # Check there are not multiple start tags and warn user if found:
      if(length(Start) > 1) stop(paste("Multiple opening ", Tag, " tags.", sep = ""))
      
      # Check end tag is even found and stop and warn user if not:
      if(length(End) == 0) stop(paste("Missing closing ", Tag, " tag.", sep = ""))
      
      # Check there are not multiple end tags and warn user if found:
      if(length(End) > 1) stop(paste("Multiple closing ", Tag, " tags.", sep = ""))
      
      # Collapse text to just the part inside the current tags:
      RawText <- RawText[Start:End]
      
      # Get any data from inside the tag:
      InsideTagData <- RetrieveTagSupplement(Tag, RawText)
      
      # Get basic count of tags in current text (ignoring IgnoreTags tags):
      TagCount <- sum(sort(unlist(strsplit(gsub(paste(IgnoreTags, collapse = "|"), "", RawText), split = ""))) == "<")
      
      # If is only the current tag:
      if(TagCount == 2) {
        
        # Extract tag contents and store as output:
        Output <- base::gsub(paste(".*<", base::gsub("<|>", "", Tag), "*>|<", base::gsub("<|>", "", ClosingTag), ">.*", sep = ""), "", RawText)
        
      # If potential subtags found (assumes either <List> or <Type>):
      } else {
        
        # Find any <List> tag positions:
        ListPositions <- base::grep("<List", RawText)
        
        # Find any <Type> tag positions:
        TypePosition <- base::grep("<Type", RawText)
        
        # Check both are not used and stop and warn user if found:
        if(length(ListPositions) > 0 && length(TypePosition) > 0) stop("Cannot use both <List> and <Type> subtags for the saem nesting tag.")
        
        # If <List> tags are used:
        if(length(ListPositions) > 0) {
          
          # Get lsit names (contents of list tags):
          ListNames <- unlist(lapply(as.list(RawText[ListPositions]), function(x) base::gsub(".*>", "", base::gsub("</List>.*", "", x))))

          # If there is supplemental content to the List tags:
          if(length(base::grep("<List ", RawText)) > 0) {
            
            # Extract that supplement:
            ListSupplement <- lapply(as.list(RawText[ListPositions]), function(x) RetrieveTagSupplement("<List>", x))
            
            # Find intersecting values across each list tag:
            IntersectingValues <- base::Reduce(base::intersect, lapply(ListSupplement, function(x) x[, "Measure"]))
            
            # Find unique values across each list tag:
            UniqueValues <- unique(unlist(lapply(ListSupplement, function(x) x[, "Measure"])))
            
            # Check these match up (if not suggests incosistent use of list tag and hence stop and warn user:
            if(paste(sort(UniqueValues), collapse = "") != paste(sort(IntersectingValues), collapse = "")) stop("Supplemental values inside list tags are not consistent. Check these and try again.")
            
            # Build output into a matrix:
            Output <- cbind(matrix(unlist(lapply(ListSupplement, function(x) x[, "Value"])), ncol = length(UniqueValues), byrow = TRUE, dimnames = list(c(), UniqueValues)), ListNames)
            
            # Add column name to last column:
            colnames(Output)[ncol(Output)] <- "ListValue"
            
          # Case if no supplemental information in list tag:
          } else {
            
            # Build a one-column matrix of list values:
            Output <- matrix(ListNames, ncol = 1, dimnames = list(c(), "ListValue"))
            
          }
          
        }
        
        # If type tags were found:
        if(length(TypePosition) > 0) {
          
          # Built type value matrix:
          Output <- matrix(base::gsub(".*<Type>|</Type>.*", "", RawText[TypePosition]), ncol = 1, dimnames = list(c(), "TypeValue"))
          
        }
        
        # Check for case where neither <List> nor <Type> is used; stop and warn user:
        if(length(ListPositions) == 0 && length(TypePosition) == 0) stop("Unexpected subtags found, but these are either:\n1. Missing from the internal tags database (function requires updating).\n2. Not drawn from the current generic list (only <List> or <Type>).\n3. Not really subtags (but a linebreak should have been employed).\n4. Not really tags but use of the greater than/less than symbols inside a tag.")
        
      }
      
      # Form tag supplement from inside tag data (NULL if none found):
      ifelse(length(InsideTagData) > 0, TagSupplement <- InsideTagData, TagSupplement <- list(NULL))
      
      # Reform output as two item list:
      Output <- list(Output, TagSupplement)
      
      # Add anems to output:
      names(Output) <- c("TagContents", "TagSupplement")

    }
    
    # Return output:
    return(Output)
    
  }
  
  # Form tasg to extract from just nested tags (higher-level tags simply contain these and would create redundancy issues if parsed):
  TagsToExtract <- base::setdiff(TagsDatabase[, "Tag"], TagsDatabase[, "Nesting"])
  
  # Build read tags (will still need to store these but this is faster than a for loop):
  ReadTags <- lapply(as.list(TagsToExtract), function(x) TagReader(x, RawText))
  
  # get matching paths t store the data in ReadTags:
  PathsToStore <- TagPaths[unlist(lapply(as.list(TagsToExtract), function(x) which(TagsDatabase[, "Tag"] == x)))]
  
  # For each tag to extract:
  for(i in 1:length(TagsToExtract)) {
    
    # If a NULL value was returned:
    if(is.null(ReadTags[[i]])) {
      
      # Store NULL in tag value:
      base::eval(base::parse(text = paste("Output$", PathsToStore[i], " <- list(NULL)", sep = "")))
      
    # If a non-NULL value was returned:
    } else {
      
      # Store values in tag:
      base::eval(base::parse(text = paste("Output$", PathsToStore[i], " <- ReadTags[[i]]", sep = "")))
      
    }
    
  }
  
  # Return output (inviisbly or not depending on what is requested):
  if(Invisible) {invisible(Output)} else {return(Output)}

}
