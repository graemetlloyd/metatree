#' Reads in a metatree-format XML file
#'
#' @description
#'
#' Reads in a metatree-formatted XML file.
#'
#' @param File Path to XML file to read.
#' @param Invisible Logical indicating whether to print output to screen (FALSE) or not (TRUE, the default).
#'
#' @details
#'
#' \bold{Introduction}
#'
#' The two main file inputs to the \link{Metatree} function are an MRP file, summarising the unique bipartitions found amongst the source trees, and an XML file containing important metadata about the data set. This file must be in a very specific format to be used by the \link{Metatree} function and that format is described in detail here.
#'
#' Note that this is the format the metatree approach (Lloyd et al. 2016) is based on and borrows heavily from the Supertree Toolkit format of Hill and Davis (2014).
#'
#' \bold{Example data set}
#'
#' Multiple examples of this format are available at \href{http://www.graemetlloyd.com/matr.html}{graemetlloyd.com}, but a simple version is also shown below.
#'
#' \preformatted{<?xml version="1.0" standalone="yes"?>
#' <SourceTree>
#'   <Source>
#'     <Author>
#'       <List>Michaux, B.</List>
#'     </Author>
#'     <Year>1989</Year>
#'     <Title>Cladograms can reconstruct phylogenies: an example from the fossil record</Title>
#'     <Journal>Alcheringa</Journal>
#'     <Volume>13</Volume>
#'     <Pages>21-36</Pages>
#'     <Booktitle/>
#'     <Publisher/>
#'     <City/>
#'     <Editor/>
#'   </Source>
#'   <Taxa number="4">
#'     <List recon_name="Ancilla" recon_no="10760">Ancilla</List>
#'     <List recon_name="DELETE" recon_no="-1">Turrancilla</List>
#'     <List recon_name="Ancillista" recon_no="10763">Ancillista</List>
#'     <List recon_name="Amalda" recon_no="10743">Amalda</List>
#'   </Taxa>
#'   <Characters>
#'     <Molecular/>
#'     <Morphological number="11">
#'       <Type>Osteology</Type>
#'     </Morphological>
#'     <Behavioural/>
#'     <Other/>
#'   </Characters>
#'   <Analysis>
#'     <Type>Maximum Parsimony</Type>
#'   </Analysis>
#'   <Notes>Based on reanalysis of the original matrix.</Notes>
#'   <Filename>Michaux_1989aa</Filename>
#'   <Parent/>
#'   <Sibling/>
#' </SourceTree>}
#'
#' The XML format is composed of a series of "tags", typically with opening (e.g., <List>) and closing (e.g., </List>) pairs, with unused tags containg a lsahs at the end (e.g., <Sibling/>).
#'
#' The main tags are visited in order below.
#'
#' \bold{XML tag}
#'
#' \preformatted{<?xml version="1.0" standalone="yes"?>}
#'
#' This simply states the filetype for machine-readable purposes and should not vary at all.
#'
#' \bold{SourceTree tag}
#'
#' \preformatted{<SourceTree> ... </SourceTree>}
#'
#' This tag envelops the whole rest of the file and should always be used.
#'
#' \bold{Source tag}
#'
#' \preformatted{<Source> ... </Source>}
#'
#' This tag contains the information regarding the source reference and contains multiple subtags, some of which will not always be employed (e.g., depending on whether the reference corresponds to a journal article or a book chapter).
#'
#' \emph{Author tag}
#'
#' \preformatted{<Author> ... </Author>}
#'
#' For the authors of the work. Multiple authors are allowed and each should be placed in separate <List> tags.
#'
#' \emph{Year tag}
#'
#' \preformatted{<Year> ... </Year>}
#'
#' The year of publication.
#'
#' \emph{Title tag}
#'
#' \preformatted{<Title> ... </Title>}
#'
#' The title of the work. If a book (but not a book chapter) then use <Booktitle> instead (see below).
#'
#' \emph{Journal tag}
#'
#' \preformatted{<Journal> ... </Journal>}
#'
#' The journal name (if a journal article).
#'
#' \emph{Volume tag}
#'
#' \preformatted{<Volume> ... </Volume>}
#'
#' The volume number (for journal articles).
#'
#' \emph{Pages tag}
#'
#' \preformatted{<Pages> ... </Pages>}
#'
#' The page numbers (or article number for some newer journal types).
#'
#' \emph{Booktitle tag}
#'
#' \preformatted{<Booktitle> ... </Booktitle>}
#'
#' The book title if a book or book chapter.
#'
#' \emph{Publisher tag}
#'
#' \preformatted{<Publisher> ... </Publisher>}
#'
#' The publisher, if a book or book chapter.
#'
#' \emph{City tag}
#'
#' \preformatted{<City> ... </City>}
#'
#' The city of publication, if a book or book chapter.
#'
#' \emph{Editor tag}
#'
#' \preformatted{<Editor> ... </Editor>}
#'
#' The editor, or editors (for book chapters). Like authors there can be multiple so each should be included in a <List> tag.
#'
#' \bold{Taxa tag}
#'
#' \preformatted{<Taxa> ... </Taxa>}
#'
#' This tag contains the information on the taxa (Operational Taxonomic Units; OTUs), including the reconciliation between these and the Paleobiology Database \href{https://paleobiodb.org/}{paleobiodb.org}.
#'
#' Note that the total number of taxa should be included in the \emph{opening} tag, e.g., <Taxa number="12">.
#'
#' Individual taxa should then be included in <List> tags. Note that the name inside the tag MUST match exactly the names in the MRP file and should not be manually edited to avoid this. The \link{Metatree} function will check for this, but the
#'
#' Information on the reconciliation is included inside each opening <List> tag using both the "recon_name" and "recon_no" values (e.g., <List recon_name="Amalda" recon_no="10743">). By default these should be recon_name="DELETE" and recon_no="-1", indicating the OTU has not yet been reconciled with the Paleobiology Database and should thus be pruned during the metatree construction process. (I.e., this is much safer than assigning a random real taxon to each OTU.) In operation the "DELETE" value will lead to these taxa being pruned inside the \link{Metatree} function and the "-1" indicates that the taxon has not been checked yet. A user may wish to set a taxon to DELETE \emph{after} checking, for example because it represents a hypothetical outgroup not a real species. In these case the recon_no should be set to "0" instead. Note that "0" and "1" are used here as reserved values because the Paleobiology Database numbering system begins at "1" and hence any number greater than "0" will (potentially) represent a real taxon.
#'
#' However, for the \link{Metatree} function to work in any meaningful way the majority of OTUs must be reconciled with the Paleobiology Database, including BOTH the taxon name(s) and taxon number(s). The reason for this apparent redundancy is to ensure data integrity and is predicated on multiple considerations. Firstly, names are not unique and can exist multiple times because, for example, they are used for both an animal and a plant, they are used separately to denote higher or lower taxa (allowed within ICZN rules), or simply there is an uncorrected homonym issue in the Paleobiology Database. Thus names should never be used in isolation as many mishaps may befall you if you do.
#'
#' Numbers could theoretically be used in isolation, but two dangers arise here. First, a typographical error is much more easily made without being spotted as human beings are more attuned to spot say, Tyrannosuarus, instead of Tyrannosaurus than 314567, instead of 314657. Second, a name corresponding to a specific taxon number can be updated or edited and cross-validation is thus required to ensure names and numbers match. Again, the \link{Metatree} function will check for these issues as it operates.
#'
#' Many users may not be aware of how to find taxon numbers in the database, but a simple way is to search for the required name and find the corresponding page in the database. For example, the page for \emph{Tyrannosaurus rex} is \href{https://paleobiodb.org/classic/basicTaxonInfo?taxon_no=54833}{here}. Look at the URL in your browser and you should see it ends (or contains) \code{taxon_no=54833}. Thus we could reconcile the OTU "Trex" as follows:
#'
#' \preformatted{<List recon_name="Tyrannosaurus_rex" recon_no="54833">Trex</List>}
#'
#' Note that recon_name is the full name, as it is spelled in the Paleobiology Database, and using underscores where spaces would otherwise exist (i.e., between the genus and the species).
#'
#' This might seem like a laborious process if you have hundreds of tips (or many thousands across a series of input data sets), but without careful manual checking for taxonomic reconciliation any composite analysis will be confounded. (I can point you in the direction of some truly awful examples if you ask me but will not publicly shame the guilty parties here.) Thus here I do not provide, or encourage, attempts to automate this process.
#'
#' Taxonomic reconciliation (matching OTU names to valid and appropriate taxa) is obviously critical and hence some best practice guidelinea are offered here:
#'
#' 1. \bold{Endeavour to use species-level reconciliations at all times}. It has become common practice amongst many palaeontological authors to only write the genus name in phylogenetic analyses, instead of the species actually examined. This can lead to unintended consequences when the contents of those genera shift (e.g., track the history of \emph{Brontosaurus}). Using species-level names avoids this issue, allowing a dynamic taxonomy (i.e., the Paleobiology Database) to take care of these updates for you. If a supra-specific system is used instead then the user is doomed to repeat the same manual checks over and over again. In other words, if you reconcile with a species name at the beginning then you never need to update that reconciliation again.
#'
#' 2. \bold{Never manually perform synonymisation}. In other words, taxa should be entered as the original authors intended and any synonymisations done automatically by calls to a single database (here, the Paleobiology Database). The same goes for nomen dubia and the like. The reason is simply that this is not a sustainable aproach and will require contnued re-examination of XMLs and ultimately errors will creep in, compromising the data. Again, reconcile once correctly when the file is created and it never needs to be revisited.
#'
#' 3. \bold{Multiple taxon reconciliation is possible}. In some cases authors will explicitly code a single OTU from multiple specimens or species. This can be accounted for two by assigning multiple species to that OTU. For example, let's say an author codes both species of \emph{Unenlagia} as a single OTU (Unenlagia). This means there should really be two recon_name values and two recon_no values. This is dealt with by using commas and semicolons respectively, i.e.:
#'
#' \code{<List recon_name="Unenlagia_comahuensis,Unenlagia_paynemili" recon_no="65422;65423">Unenlagia</List>}
#'
#' Note that the order matters here (the numbers and names must be in the correct order). Again, the \link{Metatree} function will check this and warn the user, but it is always less work to get it right the first time.
#'
#' It is also important to remember that in doing multiple taxon reconciliations this way the metatree process will assume the OTU is monophyletic, i.e., that the species involved form a clade. It can create downstream problems if this is not the case so multiple taxon reconciliations should always be performed carefully.
#'
#' Overall the key thing to remember is that taxonomy is not static and that a well designed taxonomic reconciliation process will take this into account: that is the intended aim here.
#'
#' \bold{Characters tag}
#'
#' \preformatted{<Characters> ... </Characters>}
#'
#' This tag contains very limited information on the characters used in the matrix, including total number and type. Currently the \link{Metatree} function does not use this data directly, although it may be used in future to apply some automated filtering (e.g., morphology only, exclude MRP etc.). Thus at present this is not a critical field to fill in, but may become so in future.
#'
#' Four main overall types are currently included (molecular, morphological, behavioural and other) and the number of each should be included in the opening tag (as long as there are at least one), followed by the more specific type(s) in <Type> tags. E.g.:
#'
#' \preformatted{<Morphological number="58">
#'   <Type>Osteology</Type>
#' </Morphological>}
#'
#' \emph{Molecular tag}
#'
#' \preformatted{<Molecular> ... </Molecular>}
#'
#' The number and type (e.g., mtDNA, RAG1 etc.) of molecular characters.
#'
#' \emph{Molecular tag}
#'
#' \preformatted{<Morphological> ... </Morphological>}
#'
#' The number and type (e.g., osteology, dermal etc.) of morphological characters.
#'
#' \emph{Behavioural tag}
#'
#' \preformatted{<Behavioural> ... </Behavioural>}
#'
#' The number and type (e.g., nesting style, diurnality etc.) of behavioural characters.
#'
#' \emph{Other tag}
#'
#' \preformatted{<Behavioural> ... </Behavioural>}
#'
#' The number and type (e.g., MRP, geographic etc.) of any other characters.
#'
#' \bold{Analysis tag}
#'
#' \preformatted{<Analysis> ... </Analysis>}
#'
#' This tag contains inforamtion on the type of analysis performed to generate the MRP data (e.g., parsimony, Bayesian, likelihood). The specific approach should appear in a single <Type> tag. E.g.:
#'
#' \preformatted{<Analysis>
#'   <Type>Maximum Parsimony</Type>
#' </Analysis>}
#'
#' \bold{Notes tag}
#'
#' \preformatted{<Notes> ... </Notes>}
#'
#' Simply any notes the user may want to append to the file.
#'
#' \bold{Filename tag}
#'
#' \preformatted{<Filename> ... </Filename>}
#'
#' The filename used (excluding extension). Ideally this should match across the XML and MRP files and be formatted as a pseudo-citation. For example, if the citation would be Rogers et al. (2012), the filename would be Rogers_etal_2012. Additonally, it is helpful to append a lowercase letter to the end of files in case these names end up being non-unique, i.e., Rogers_etal_2012a for the first citation, Rogers_etal_2012b for the second, and so on. Finally, some references may include multiple analyses and so these can be distinguished by using an additional lowercase letter, e.g., Rogers_etal_2012aa and Rogers_etal_2012ab etc.
#'
#' (NB: The disadvantage of such a system is that it will breakdown as soon as there are 27 or more duplicate citations or data sets from a single reference. In practice this has not become a problem yet, but I might review and change this in future.)
#'
#' \bold{Parent tag}
#'
#' \preformatted{<Parent> ... </Parent>}
#'
#' The filename of any other data set in the sample that can logically be considered the parent of the current data set. This tag is used to deal with the fact that many morphological data sets are not independent (<Parent/> would be used if they were), but are based wholly or primarily on some older data set with little or no modification (e.g., adding a single new row (taxon) or updating some codings). This information is used within \link{Metatree} to prune redundnat data sets and reweight remaining non-indpendent ones and so is critical to the metatree construction process.
#'
#' \bold{Sibling tag}
#'
#' \preformatted{<Sibling> ... </Sibling>}
#'
#' Similar to the parent tag, this tag is used to denote data sets with equal claim to priority, but that do not represent a parent-child relationship. These tend to be rarer in my experience, but do crop up occasionally when, for example, two different coding schemes are applied.
#'
#' @return
#'
#' A nested list reflecting the nested XML tags of the input file.
#'
#' @author
#'
#' Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Hill, J. and Davis, K. E., 2014. The Supertree Toolkit 2: a new and improved software package with a Graphical User Interface for supertree construction. \emph{Biodiversity Data Journal}, \bold{2}, e1053.
#'
#' Lloyd, G. T., Bapst, D. W., Friedman, M., and Davis, K. E., 2016. Probabilistic divergence time estimation without branch lengths: dating the origins of dinosaurs, avian flight and crown birds. \emph{Biology Letters}, \bold{12}, 20160609.
#'
#' @seealso
#'
#' \link{WriteMetatreeXML}.
#'
#' @examples
#'
#' # Example line (that would print to screen):
#' #ReadMetatreeXML("Rogers_etal_2012a.xml", Invisible = FALSE)
#'
#' # (Note that this is commented out as it would only work locally,
#' # but should give the user an idea of the syntax)
#'
#' @export ReadMetatreeXML
ReadMetatreeXML <- function(File, Invisible = TRUE) {
  
  # TO DO:
  
  # SEPARATE TEST THAT FILE BEGINS WITH ?xml .. (WEIRD TAG AS HAS NO PAIRED CLOSING TAG)
  # IDEALLY WILL ALSO TEST FOR PRESENCE OF APPROPRIATE LINEBREAKS (DO NOT WANT DIFFERENT TAGS ON SAME LINE AT ANY POINT IN FILE)
  # ADD MAKE METATREE FUNCTION
  
  # OTHER CHECKS TO MAYBE ADD:
  #
  # 1. That taxon count matches number of taxa actually listed.
  # 2. That total number of characters matches external NEXUS file.
  # 3. Other checks from metatree master (spaces in names, numbers etc.).
  # 4. Author names are formatted correctly.
  # 5. References are formatted correctly overall.

  # Build simple tags database (can be edited in future) to set expected structure to be found in file - does not include generic <List> or <Type> subtags:
  TagsDatabase <- matrix(c("<SourceTree>", "", "<Source>", "<SourceTree>", "<Author>", "<Source>", "<Year>", "<Source>", "<Title>", "<Source>", "<Journal>", "<Source>", "<Volume>", "<Source>", "<Pages>", "<Source>", "<Booktitle>", "<Source>", "<Publisher>", "<Source>", "<City>", "<Source>", "<Editor>", "<Source>", "<Taxa>", "<SourceTree>", "<Characters>", "<SourceTree>", "<Molecular>", "<Characters>", "<Morphological>", "<Characters>", "<Behavioural>", "<Characters>", "<Other>", "<Characters>", "<Analysis>", "<SourceTree>", "<Notes>", "<SourceTree>", "<Filename>", "<SourceTree>", "<Parent>", "<SourceTree>", "<Sibling>", "<SourceTree>"), ncol = 2, byrow = TRUE, dimnames = list(c(), c("Tag", "Nesting")))
  
  # Build list of tags to ignore if foudn insde other tags (e.g., typeface tags like bold or emphasis):
  IgnoreTags <- c("<em>", "</em>")
  
  # Genearte vector of tag paths (will allow easy access of list later):
  TagPaths <- unlist(lapply(as.list(TagsDatabase[, "Tag"]), function(x) {CurrentTag <- x; CurrentTag <- c(CurrentTag, base::unname(TagsDatabase[TagsDatabase[, "Tag"] == CurrentTag, "Nesting"])); while(!any(CurrentTag == "")) CurrentTag <- c(CurrentTag, base::unname(TagsDatabase[TagsDatabase[, "Tag"] == CurrentTag[length(CurrentTag)], "Nesting"])); CurrentTag <- base::setdiff(CurrentTag, ""); paste(rev(base::gsub("<|>", "", CurrentTag)), collapse = "$")}))
  
  # Subfunction to convert two-column matrix of nestings (X nested in Y) to nested list:
  NestingMatrixToNestedList <- function(NestingMatrix) {
    
    # FOLLOWING WILL BREAK IF ORDER OF NESTING IS NOT LEAST TO MOST NESTED IN TAGS DATABASE MATRIX SO MAYBE SORT THIS? ALTHOUGH THAT WOULD BREAK ORDER IN XML, SO MAYBE GET ORDER FROM XML ITSELF?

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
