#' Writes a metatree-format XML file
#'
#' Writes out a metatree-format XML file (Lloyd et al. 2016).
#'
#' @param XML The XML file to write. Must be in format imported by \link{ReadMetatreeXML}.
#' @param File Path to XML file to write.
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
#' @export WriteMetatreeXML
WriteMetatreeXML <- function(XML, File) {
  
  # Build simple tags database:
  TagsDatabase <- matrix(c("<SourceTree>", "", "<Source>", "<SourceTree>", "<Author>", "<Source>", "<Year>", "<Source>", "<Title>", "<Source>", "<Journal>", "<Source>", "<Volume>", "<Source>", "<Pages>", "<Source>", "<Booktitle>", "<Source>", "<Publisher>", "<Source>", "<City>", "<Source>", "<Editor>", "<Source>", "<Taxa>", "<SourceTree>", "<Characters>", "<SourceTree>", "<Molecular>", "<Characters>", "<Morphological>", "<Characters>", "<Behavioural>", "<Characters>", "<Other>", "<Characters>", "<Analysis>", "<SourceTree>", "<Notes>", "<SourceTree>", "<Filename>", "<SourceTree>", "<Parent>", "<SourceTree>", "<Sibling>", "<SourceTree>"), ncol = 2, byrow = TRUE, dimnames = list(c(), c("Tag", "Nesting")))
  
  # Genearte vector of tag paths (will allow easy access of list later):
  TagPaths <- unlist(lapply(as.list(TagsDatabase[, "Tag"]), function(x) {CurrentTag <- x; CurrentTag <- c(CurrentTag, base::unname(TagsDatabase[TagsDatabase[, "Tag"] == CurrentTag, "Nesting"])); while(!any(CurrentTag == "")) CurrentTag <- c(CurrentTag, base::unname(TagsDatabase[TagsDatabase[, "Tag"] == CurrentTag[length(CurrentTag)], "Nesting"])); CurrentTag <- base::setdiff(CurrentTag, ""); paste(rev(base::gsub("<|>", "", CurrentTag)), collapse = "$")}))
  
  # Get just the nested tags:
  TagsToBuild <- base::setdiff(TagsDatabase[, "Tag"], TagsDatabase[, "Nesting"])
  
  # Get the opposite (nest tags) in revese ordr (so do most nested first):
  NestTags <- base::rev(base::setdiff(TagsDatabase[, "Tag"], TagsToBuild))
  
  # Create empty output vector:
  Output <- vector(mode = "character")
  
  # For each nested tag:
  for(i in TagsToBuild) {
    
    # Set j as position of tag in database (used below to get right path):
    j <- which(TagsDatabase[, "Tag"] == i)
    
    # Get current tag list:
    CurrentList <- base::eval(base::parse(text = paste("XML$", TagPaths[j], sep = "")))
    
    # If tag was used:
    if(any(names(CurrentList) == "TagContents")) {
      
      # Form close tag:
      CloseTag <- gsub("<", "</", i)
      
      # If there is no tag supplement:
      if(any(is.null(unlist(CurrentList$TagSupplement)))) {
        
        # Form start tag:
        StartTag <- i
        
      # If there is a tag supplement:
      } else {
        
        # Form start tag with supplement included:
        StartTag <- paste(gsub(">", "", i), " ", paste(apply(CurrentList$TagSupplement, 1, function(x) paste(x[1], "=\"", x[2], "\"", sep = "")), sep = " ", collapse = " "), ">", sep = "")
        
      }
      
      # Case if <List> or <Type> subtag used:
      if(is.matrix(CurrentList$TagContents)) {
        
        # Get tag type (List or Type):
        TagType <- gsub("Value", "", colnames(CurrentList$TagContents)[ncol(CurrentList$TagContents)])
        
        # If supplemental tag information exists:
        if(ncol(CurrentList$TagContents) > 1) {
          
          # Get tag supplement values:
          TagSupplementValues <- colnames(CurrentList$TagContents)[-ncol(CurrentList$TagContents)]
          
          # Build full tags and add them to the output:
          Output <- c(Output, paste(StartTag, "\n", paste(paste(paste("\t<", TagType, " ", do.call(paste, lapply(as.list(TagSupplementValues), function(x) paste(x, "=\"", CurrentList$TagContents[, x], "\"", sep = ""))), ">", sep = ""), CurrentList$TagContents[, ncol(CurrentList$TagContents)], paste("</", TagType, ">", sep = ""), sep = ""), collapse = "\n"), "\n", CloseTag, sep = ""))
          
        # If no supplemental tag information exists:
        } else {
          
          #
          Output <- c(Output, paste(StartTag, paste(paste(paste("\t<", TagType, ">", sep = ""), CurrentList$TagContents[, 1], paste("</", TagType, ">", sep = ""), sep = ""), collapse = "\n"), CloseTag, sep = "\n"))
          
        }
        
      # Case if no subtag used:
      } else {
        
        # Add tag to output:
        Output <- c(Output, paste(StartTag, CurrentList$TagContents, CloseTag, sep = ""))
        
      }

    # Case if tag is NULL (unused):
    } else {
      
      # Make closed (/>) tag and store in output:
      Output <- c(Output, gsub(">", "/>", i))
      
    }
    
  }
  
  # Pull out linebreaks (will allow tabs to be added properly later):
  Output <- unlist(strsplit(Output, split = "\n"))
  
  # For each nest tag:
  for(i in NestTags) {
    
    # Find the tags nested inside it:
    NestedTags <- TagsDatabase[TagsDatabase[, "Nesting"] == i, "Tag"]
    
    # Find start index of first tag:
    StartIndex <- grep(gsub(">", "", NestedTags[1]), Output)
    
    # Find end index of last tag:
    EndIndex <- c(grep(gsub("<", "</", NestedTags[length(NestedTags)]), Output), grep(gsub(">", "/>", NestedTags[length(NestedTags)]), Output))
    
    # Add extra tab to these values as they will beindented inside the current tag:
    Output[StartIndex:EndIndex] <- paste("\t", Output[StartIndex:EndIndex], sep = "")
    
    # Add opening nest tag to beginning of block:
    Output[StartIndex] <- paste(i, "\n", Output[StartIndex], sep = "")
    
    # Add closing nest tag to end of block:
    Output[EndIndex] <- paste(Output[EndIndex], "\n", gsub("<", "</", i), sep = "")
    
    # Resplit by linebreaks so ready for next level of nesting:
    Output <- unlist(strsplit(Output, split = "\n"))

  }
  
  # Add XML header to file:
  Output <- c("<?xml version=\"1.0\" standalone=\"yes\"?>", Output)
  
  # Write XML to file:
  write(Output, File)

}

#XML <- ReadMetatreeXML("~/Documents/Homepage/www.graemetlloyd.com/xml/Zhu_et_Ahlberg_2004a.xml")
#File <- "Zhu_et_Ahlberg_2004a.xml"
#WriteMetatreeXML(XML, File)
