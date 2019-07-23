#' Builds metatree from source data
#'
#' Builds metatree from source data
#'
#' Use approach laid out in Lloyd et al. (2008; 2016).
#'
#' @param MRPDirectory The directory in which the MRP files are to be read from.
#' @param XMLDirectory The directory in which the XML files are to be read from.
#' @param TargetClade The name of the target clade of the metatree (e.g., "Dinosauria").
#' @param InclusiveDataList A vector of the data sest to include in the metatree. Can be left empty to just read all files in \code{MRPDirectory} abd \code{XMLDirectory}.
#' @param ExclusiveDataList A vector of any data sets to exclude from the metatree. Can be left empty if all data sets in \code{MRPDirectory} abd \code{XMLDirectory} are valid. (Intended to exclude things like oogenera or footprint analyses, other supertree data sets etc.)
#' @param HigherTaxaToCollapse Vector of any higher taxa to collapse (e.g., if you are focused on relationships in a stem-group).
#' @param MissingSpecies What to do with species assigned to the target clade, but not present in the source data. Options are: "exclude" (excludes these missing species), "genus" (include those species in a genus-level polytomy if the genus is sampled in the source data), and "all"
#' @param Interval If restricting the sample to a specific interval of geologic time then use this option (passed to \link{PaleobiologyDBChildFinder} which should be consulted for formatting). Default is NULL (no restriction on ages of tips to be included).
#' @param VeilLine A logical indiicating whther to remove older data sets that do not increase taxonomic coverage (TRUE; the default) or not (FALSE). See Lloyd et al. (2016).
#' @param SpeciesToExclude Vector of any species to be excluded from the final metatree. E.g., Eshanosaurus, Ricardoestesia.
#' @param IncludeSpecimenLevelOTUs A logical indicating whther specimen-level OTUs should (TRUE; the default) or should not (FALSE) be included in the metatree.
#' @param BackboneConstraint A Newick string of a backbone constraint (will enforce topology in final metatree but allows taxa not in topology to fall out inside the constraint). This is not required and the default (NULL) will mean no constraint is applied.
#' @param MonophylyConstraint A Newick string of a monophyly constraint (will enforce topology in final metatree and force taxa not in topology to fall outside the constraint). This is not required and the default (NULL) will mean no constraint is applied.
#'
#' @return TBC.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Lloyd et al. (@016)
#'
#' @examples
#'
#' # Nothing yet.
#'
#' @export Metatree
Metatree <- function(MRPDirectory, XMLDirectory, TargetClade = "", InclusiveDataList = c(), ExclusiveDataList = c(), HigherTaxaToCollapse = c(), MissingSpecies = "exclude", Interval = NULL, VeilLine = TRUE, SpeciesToExclude = c(), IncludeSpecimenLevelOTUs = TRUE, BackboneConstraint = NULL, MonophylyConstraint = NULL) {
  
  MRPDirectory <- "/Users/eargtl/Documents/Homepage/www.graemetlloyd.com/mrp"
  XMLDirectory <- "/Users/eargtl/Documents/Homepage/www.graemetlloyd.com/xml"
  TargetClade <- "Ichthyopterygia"
  InclusiveDataList <- sort(c(GetFilesForClade("matricht.html"), "Bickelmann_etal_2009a", "Caldwell_1996a", "Chen_etal_2014ba", "Chen_etal_2014bb", "deBraga_et_Rieppel_1997a", "Gauthier_etal_1988b", "Laurin_et_Reisz_1995a", "Muller_2004a", "Reisz_etal_2011a", "Rieppel_et_Reisz_1999a", "Rieppel_et_deBraga_1996a"))
  ExclusiveDataList <- c("Averianov_inpressa", "Bravo_et_Gaete_2015a", "Brocklehurst_etal_2013a", "Brocklehurst_etal_2015aa", "Brocklehurst_etal_2015ab", "Brocklehurst_etal_2015ac", "Brocklehurst_etal_2015ad", "Brocklehurst_etal_2015ae", "Brocklehurst_etal_2015af", "Bronzati_etal_2012a", "Bronzati_etal_2015ab", "Brusatte_etal_2009ba", "Campbell_etal_2016ab", "Carr_et_Williamson_2004a", "Carr_etal_2017ab", "Frederickson_et_Tumarkin-Deratzian_2014aa", "Frederickson_et_Tumarkin-Deratzian_2014ab", "Frederickson_et_Tumarkin-Deratzian_2014ac", "Frederickson_et_Tumarkin-Deratzian_2014ad", "Garcia_etal_2006a", "Gatesy_etal_2004ab", "Grellet-Tinner_2006a", "Grellet-Tinner_et_Chiappe_2004a", "Grellet-Tinner_et_Makovicky_2006a", "Knoll_2008a", "Kurochkin_1996a", "Lopez-Martinez_et_Vicens_2012a", "Lu_etal_2014aa", "Norden_etal_inpressa", "Pisani_etal_2002a", "Ruiz-Omenaca_etal_1997a", "Ruta_etal_2003ba", "Ruta_etal_2003bb", "Ruta_etal_2007a", "Selles_et_Galobart_2016a", "Sereno_1993a", "Sidor_2001a", "Skutschas_etal_inpressa", "Tanaka_etal_2011a", "Toljagic_et_Butler_2013a", "Tsuihiji_etal_2011aa", "Varricchio_et_Jackson_2004a", "Vila_etal_2017a", "Wilson_2005aa", "Wilson_2005ab", "Zelenitsky_et_Therrien_2008a")
  HigherTaxaToCollapse = c()
  MissingSpecies = "exclude"
  Interval = NULL
  VeilLine = TRUE
  SpeciesToExclude = c()
  IncludeSpecimenLevelOTUs = TRUE
  BackboneConstraint = NULL
  MonophylyConstraint = NULL
  
  # New Options (requires code to actually use them)
  #
  # HigherTaxaToCollapse Vector can be empty.
  # VeilLine TRUE/FALSE (will be in output)
  # SpeciesToExclude Vector of any species to be excluded from the final metatree. E.g., Eshanosaurus, Ricardoestesia.
  # IncludeSpecimenLevelOTUs TRUE/FALSE
  # BackboneConstraint Newick string of backbone constraint (allows taxa not in topology). NULL as default.
  # MonophylyConstraint Newick string of monophyly constraint (excludes taxa not in topology). NULL as default.
  
  # CHECK PARENT IS A DATA SET AND NOT A REFERENCE, E.G., IF ENTER A REFERENCE AS PARENT THEN PARENT TURNS OUT TO HAVE TWO DATA SETS
  # CHECK FOR SPECIES THAT BELONG TO A GENUS DIFFERENT TO THE ONE IN THEIR NAME!
  # NEED TO CATCH ISSUE WHERE GENUS NUMBER IS USED FOR A SPECIES (HARD TO CHECK SO FAR DUE TO INDETERMINATES CONTINGENCY)
  # NEED SOME TEST THAT HELPS CHECK ROOT IS SENSIBLE
  # NEED SOME TEST THAT HELPS DETERMINE IF MULTIPLE OCCURRENCES OF SAME TAXON AFTER RECONCILIATION IS CORRECT OR AN ERROR
  # ADD MORE COMPLEX WEIGHTS BY USING ADDITIONAL CHARACTER STATES! (EACH DATASET TOTAL WEIGHT DETERMINED BY YEAR AND DEPENDENCE THEN SUBDIVIDED ACROSS CHARACTER?) - BUT THIS SEEMS TO SLOW THINGS DRAMATICALLY MAYBE DO BY DUPLICATING CHARACTERS INSTEAD
  # MAKE STR OPTIONAL (SAVES A LITTLE TIME)
  # CHECK THERE ARE MULTIPLE TAXA PRE-RECONCILIATION
  # CHECK INDETS DO NOT GIVE MULTIPLE MATCHES
  # CHECK FOR ABSENT RECON NAMES OR NUMBERS ("")
  
  # HOW TO DELETE DATA SETS THAT STILL CONTRIBUTE TO DEPENDENCE?
  
  # Check MRPDirectory is formatted correctly adn stop and warn user if not:
  if(!all(is.character(MRPDirectory)) || length(MRPDirectory) != 1) stop("MRPDirectory must be a single character string indicating the path to the folder containing the MRP files.")
  
  # Check XMLDirectory is formatted correctly adn stop and warn user if not:
  if(!all(is.character(XMLDirectory)) || length(XMLDirectory) != 1) stop("XMLDirectory must be a single character string indicating the path to the folder containing the XML files.")
  
  # Check TargetClade is formatted correctly adn stop and warn user if not:
  if(!all(is.character(TargetClade)) || length(TargetClade) != 1) stop("TargetClade must be a single character string indicating the desired clade the metatree will represent.")
  
  # Check MissingSpecies respresents a valid option:
  if(length(setdiff(MissingSpecies, c("all", "exclude", "genus"))) > 0) stop("MissingSpecies must be one of \"all\", \"exclude\", or \"genus\".")
  
  # Check VeilLine is a logical and stop and warn user if not:
  if(!is.logical(VeilLine)) stop("VeilLine must be a logical (TRUE or FALSE).")
  
  # Check IncludeSpecimenLevelOTUs is a logical and stop and warn user if not:
  if(!is.logical(IncludeSpecimenLevelOTUs)) stop("IncludeSpecimenLevelOTUs must be a logical (TRUE or FALSE).")
  
  # If not a NULL read backbone constraint (checks it is a valid Newick format):
  if(!is.null(BackboneConstraint)) BackboneConstraintTree <- ape::read.tree(text = BackboneConstraint)
  
  # If not a NULL read monophyly constraint (checks it is a valid Newick format):
  if(!is.null(MonophylyConstraint)) MonophylyConstraintTree <- ape::read.tree(text = MonophylyConstraint)
  
  # List of types of resolution that require finding a senior synonym:
  synonyms <- c("corrected to", "misspelling of", "objective synonym of", "obsolete variant of", "recombined as", "replaced by", "subjective synonym of")
  
  # List of types of resolution that require changing reconciliation to DELETE:
  deletes <- c("nomen dubium", "nomen vanum", "nomen nudum", "nomen oblitum", "invalid subgroup of")
  
  # Print current processing status:
  cat("Reading MRP data...")
  
  # Set working directory as MRP directory:
  setwd(MRPDirectory)
  
  # List MRP files (or just use inslusivedatalist if set):
  MRPFileList <- strsplit(ifelse(is.null(InclusiveDataList), paste(setdiff(gsub("mrp\\.nex", "", list.files()), ExclusiveDataList), "mrp.nex", sep = "", collapse = "%%"), paste(setdiff(InclusiveDataList, ExclusiveDataList), "mrp.nex", sep = "", collapse = "%%")), "%%")[[1]]
  
  # Read in all MRP files and store in a list (include duplicate headers to store parent sibling info later):
  MRPList <- lapply(lapply(as.list(MRPFileList), ReadMorphNexus), function(x) list(x$Matrix_1$Matrix, x$Topper$Header, x$Topper$Header))
  
  # Set names of MRP files:
  names(MRPList) <- MRPFileList
  
  # Print current processing status:
  cat("Done\nReading XML data...")
  
  # Set working directory as XML (i.e., metadata) directory:
  setwd(XMLDirectory)
  
  # List MRP files (or just use inslusivedatalist if set):
  XMLFileList <- strsplit(ifelse(is.null(InclusiveDataList), paste(setdiff(gsub("\\.xml", "", list.files()), ExclusiveDataList), ".xml", sep = "", collapse = "%%"), paste(setdiff(InclusiveDataList, ExclusiveDataList), ".xml", sep = "", collapse = "%%")), "%%")[[1]]
  
  # Check there are no MRPs not listed as XMLs and vice versa (should return empty vector):
  MRPXMLunion <- c(setdiff(gsub("\\.xml", "", XMLFileList), gsub("mrp\\.nex", "", MRPFileList)), setdiff(gsub("mrp\\.nex", "", MRPFileList), gsub("\\.xml", "", XMLFileList)))
  
  # Stop if MRP datasets not listed as XMLs and vice versa:
  if(length(MRPXMLunion) > 0) stop(paste("Datasets do not match (MRP and XML)!:", MRPXMLunion, collapse = " "))
  
  # Update MRP names with stripped down file name:
  names(MRPList) <- gsub("mrp.nex", "", names(MRPList))
  
  # Empty vectors to store error-creating data sets:
  duplicatedtaxonnames <- namematchissues <- vector(mode = "character")
  
  # For each data set:
  for(i in XMLFileList) {
    
    # Get currentfilename:
    currentfilename <- gsub(".xml", "", i)
    
    # Extract XML text:
    XMLString <- readLines(i)
    
    # Extract taxonomic resolution text only:
    TaxonNoNameMatrix <- matrix(unlist(lapply(strsplit(unlist(lapply(strsplit(XMLString[grep("recon_name", XMLString)], "recon_name=\""), '[', 2)), "\" recon_no=\"|\"|</List>|>"), '[', c(1, 2, 4))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("PaleoDBname", "PaleoDBnumber", "OTUName")))[, c(2, 1, 3)]
    
    # Deal with subgenera formatting:
    TaxonNoNameMatrix[, "PaleoDBname"] <- gsub("_\\(|\\)", "", TaxonNoNameMatrix[, "PaleoDBname"])
    
    # Check for spaces in taxon names:
    if(length(grep(" ", TaxonNoNameMatrix[, "OTUName"])) > 0) stop(paste("Found spaces in taxon names in ", i, ".", collapse = ""))
    
    # Check for spaces in taxon numbers:
    if(length(grep(" ", TaxonNoNameMatrix[, "PaleoDBnumber"])) > 0) stop(paste("Found spaces in taxon numbers in ", i, ".", collapse = ""))
    
    # Match OTU names with initial reconciled values:
    taxonmatches <- match(rownames(MRPList[[currentfilename]][[1]]), TaxonNoNameMatrix[, "OTUName"])
    
    # Check everything does match (if not add to match issues vector):
    if(any(is.na(taxonmatches))) namematchissues <- c(namematchissues, currentfilename)
    
    # Check there are not duplicate taxa:
    if(any(duplicated(sort((taxonmatches))))) duplicatedtaxonnames <- c(duplicatedtaxonnames, currentfilename)
    
    # Perform initial reconciliation of matrix names:
    rownames(MRPList[[currentfilename]][[1]]) <- paste(TaxonNoNameMatrix[taxonmatches, "PaleoDBnumber"], TaxonNoNameMatrix[taxonmatches, "PaleoDBname"], sep = "%%%%")
    
    # Extract parent string:
    ParentString <- gsub("\t|<Parent|Parent>|<|>|/| ", "", XMLString[grep("<Parent", XMLString)])
    
    # Extract sibling string:
    SiblingString <- gsub("\t|<Sibling|Sibling>|<|>|/| ", "", XMLString[grep("<Sibling", XMLString)])
    
    # Update names of empty headers to their actual use (parent adn sibling strings):
    names(MRPList[[currentfilename]])[1:3] <- c("matrix", "parent", "sibling")
    
    # Store parent and sibling strings ("" is empty as NULL deletes them from list):
    MRPList[[currentfilename]][c("parent", "sibling")] <- c(ParentString, SiblingString)
    
  }
  
  # Stop and return braces issue (or non-matching taxa issue):
  if(length(namematchissues) > 0) stop(paste(paste("Possible missing braces (<>), rogue period(s) (.), or non-matching taxon names in ", namematchissues, " XML file.", sep = ""), collapse = "\n"))
  
  # Stop and return duplicated OTU name matrices (if any):
  if(length(duplicatedtaxonnames) > 0) stop(paste(paste("Possible duplicated taxa in ", duplicatedtaxonnames, " XML file.", sep = ""), collapse = "\n"))
  
  # Print current processing status:
  cat("Done\nChecking for unsampled parents and siblings...")
  
  # Extract parent and sibling names:
  ParentAndSiblingNames <- sort(unlist(lapply(as.list(unique(unname(unlist(lapply(MRPList, '[', c("parent", "sibling")))))), function(x) x[nchar(x) > 0])))
  
  # Warn user about any unsampled parents and/or siblings:
  if(length(setdiff(ParentAndSiblingNames, names(MRPList))) > 0) print(paste("The following parents and siblings are not in the sample (check they are correct or add them into the sample): ", paste(setdiff(ParentAndSiblingNames, names(MRPList)), collapse = ", "), sep = ""))
  
  # Print current processing status:
  cat("Done\nFinding initial multiple-taxon reconciliations...")
  
  # For each data set:
  for(i in names(MRPList)) {
    
    # Find comma rows (multiple taxa in initial reconciliation):
    commarows <- grep(",", rownames(MRPList[[i]]$matrix))
    
    # If there is at least one multiple-taxon reconciliation:
    if(length(commarows) > 0) {
      
      # For each multiple-taxon reconciliation in reverse order (to avoid later rows not matching):
      for(j in rev(commarows)) {
        
        # Get multiple names of reconciliation:
        multiplenames <- strsplit(rownames(MRPList[[i]]$matrix)[j], "%%%%")[[1]]
        
        # Get multiple-taxon numbers:
        multitaxonnumbers <- strsplit(multiplenames[1], ";")[[1]]
        
        # Get multiple-taxon names:
        multitaxonnames <- strsplit(multiplenames[2], ",")[[1]]
        
        # Check data intergirty with respect to multiple-taxon values:
        if(length(multitaxonnumbers) != length(multitaxonnames)) stop(paste("Problem with multiple-taxon reconciliation(s) in ", i, " (check commas and semi-colons are correct; i.e., of same length).", sep = ""))
        
        # Add new rows at base of matrix:
        MRPList[[i]]$matrix <- rbind(MRPList[[i]]$matrix, matrix(rep(MRPList[[i]]$matrix[j, ], length(multitaxonnumbers)), nrow = length(multitaxonnumbers), byrow = TRUE, dimnames = list(paste(multitaxonnumbers, multitaxonnames, sep = "%%%%"), c())))
        
        # Remove now redundant row from matrix:
        MRPList[[i]]$matrix <- MRPList[[i]]$matrix[-j, , drop = FALSE]
        
      }
      
    }
    
  }
  
  # Print current processing status:
  cat("Done\nRemoving taxa with initial reconciliations of \"DELETE\"...")
  
  # For each data set:
  for(i in names(MRPList)) {
    
    # Find delete rows (initial ones - may remove some later through:
    deleterows <- which(matrix(unlist(strsplit(rownames(MRPList[[i]]$matrix), "%%%%")), ncol = 2, byrow = TRUE)[, 2] == "DELETE")
    
    # If there are deletes:
    if(length(deleterows) > 0) {
      
      # Remove deleted taxon rows:
      MRPList[[i]]$matrix <- MRPList[[i]]$matrix[-deleterows, , drop = FALSE]
      
      # If less than three taxa then just delete all columns (characters) too:
      if(nrow(MRPList[[i]]$matrix) < 3) MRPList[[i]]$matrix <- MRPList[[i]]$matrix[, -c(1:ncol(MRPList[[i]]$matrix)), drop = FALSE]
      
    }
    
  }
  
  # Print current processing status:
  cat("Done\nPruning any data sets where all taxa were marked as \"DELETE\"...")
  
  # If there are any
  if(length(which(unlist(lapply(lapply(MRPList, '[[', "matrix"), ncol)) == 0)) > 0) {
    
    # Get names of data sets to be pruned:
    NamesToRemove <- names(which(unlist(lapply(lapply(MRPList, '[[', "matrix"), ncol)) == 0))
    
    # Prune from MRP List:
    MRPList <- MRPList[-match(NamesToRemove, names(MRPList))]
    
  }
  
  # Print current processing status:
  cat("Done\nSearching for and collapsing pre-reconciliation duplicated taxa...")
  
  # Create empty vector to store data sets with duplicated initially reconciled OTU names:
  duplicatedresolvedOTUs <- vector(mode = "character")
  
  # For each data set:
  for(i in names(MRPList)) {
    
    # Case if any (initially) resolved names are non-unique:
    if(any(duplicated(sort(rownames(MRPList[[i]]$matrix))))) {
      
      # Get duplicated taxon name(s):
      duplicatedtaxa <- unique(sort(rownames(MRPList[[i]]$matrix))[duplicated(sort(rownames(MRPList[[i]]$matrix)))])
      
      # For each duplicated taxon:
      for(j in duplicatedtaxa) {
        
        # Get rows for taxon:
        duplicaterows <- which(rownames(MRPList[[i]]$matrix) == j)
        
        # Build duplicated matrix from other taxa:
        tempMRPmatrix <- matrix(rep(MRPList[[i]]$matrix[-duplicaterows, ], length(duplicaterows)), ncol = ncol(MRPList[[i]]$matrix) * length(duplicaterows), dimnames = list(rownames(MRPList[[i]]$matrix)[-duplicaterows], c()))
        
        # Add duplicated taxon as single row:
        tempMRPmatrix <- rbind(tempMRPmatrix, as.vector(t(MRPList[[i]]$matrix[duplicaterows, ])))
        
        # Add duplicated taxon name:
        rownames(tempMRPmatrix)[nrow(tempMRPmatrix)] <- j
        
        # Update stored MRP matrix:
        MRPList[[i]]$matrix <- tempMRPmatrix
        
      }
      
      # Collapse to just unique characters (has to be outside the above or breaks if more than one duplicate taxon):
      MRPList[[i]]$matrix <- MRPCollapse(MakeMorphMatrix(MRPList[[i]]$matrix, header = "", weights = rep(1, ncol(MRPList[[i]]$matrix)), ordering = rep("unord", ncol(MRPList[[i]]$matrix)), equalise.weights = FALSE))$Matrix_1$Matrix
      
      # Store duplicated OTU names:
      duplicatedresolvedOTUs <- c(duplicatedresolvedOTUs, paste("Duplicated resolved OTUs (", paste(unlist(lapply(strsplit(duplicatedtaxa, "%%%%"), '[', 2)), collapse = ", "), ") in: ", i, ".\n", sep = ""))
      
    }
    
    # Format MRP as NEXUS type matrix:
    MRPmatrix <- MakeMorphMatrix(MRPList[[i]]$matrix, header = "", weights = rep(1, ncol(MRPList[[i]]$matrix)), ordering = rep("unord", ncol(MRPList[[i]]$matrix)), equalise.weights = FALSE)
    
    # Collapse MRP matrix (removes redundant characters created by taxon deletions):
    MRPList[[i]]$matrix <- MRPCollapse(MRPmatrix)$Matrix_1$Matrix
    
  }
  
  # Warn user if there are initially duplicated reconciled taxa (user should check):
  if(length(duplicatedresolvedOTUs) > 0) cat(paste(duplicatedresolvedOTUs, collapse = ""))
  
  # Print current processing status:
  cat("Done\nBuilding initial taxonomy matrix...")
  
  # Create taxonomy matrix to store all taxon resolution data:
  TaxonomyMatrix <- matrix(unlist(strsplit(sort(unique(unlist(lapply(lapply(MRPList, '[[', "matrix"), rownames)))), "%%%%")), ncol = 2, byrow = TRUE, dimnames = list(c(), c("TaxonNo", "TaxonName")))
  
  # Print current processing status:
  cat("Done\nChecking for missing taxon numbers...")
  
  # If any "-1" taxa found stop and tell user:
  if(any(TaxonomyMatrix[, "TaxonNo"] == "-1")) stop(paste("The following taxa have the reconciliation number \"-1\": ", paste(TaxonomyMatrix[TaxonomyMatrix[, "TaxonNo"] == "-1", "TaxonName"], collapse = ", "), sep = ""))
  
  # Print current processing status:
  cat("Done\nChecking for combined taxon numbers...")
  
  # If any "-1" taxa found stop and tell user:
  if(length(grep("&", TaxonomyMatrix[, "TaxonNo"]))) stop(paste("The following taxa have multiple reconciliation numbers: ", paste(TaxonomyMatrix[grep("&", TaxonomyMatrix[, "TaxonNo"]), "TaxonName"], collapse = ", "), sep = ""))
  
  # Print current processing status:
  cat("Done\nBuilding initial Paleobiology Database reconciliation list...")
  
  # Create resolved taxon numbers matrix:
  resolvedtaxonnumbers <- cbind(unique(TaxonomyMatrix[, "TaxonNo"]), PaleobiologyDBTaxaQuerier(taxon_nos = unique(TaxonomyMatrix[, "TaxonNo"]), interval = NULL))
  
  # Deal with subgenera:
  resolvedtaxonnumbers[, "TaxonName"] <- gsub(" \\(|\\)", "", resolvedtaxonnumbers[, "TaxonName"])
  
  # Add column names to first value (input number):
  colnames(resolvedtaxonnumbers)[1] <- "InputNo"
  
  # If specifying an Interval:
  if(!all(is.null(Interval))) {
    
    # Do same query for just taxa in Interval:
    resolvedtaxonnumbersInterval <- cbind(unique(TaxonomyMatrix[, "TaxonNo"]), PaleobiologyDBTaxaQuerier(taxon_nos = unique(TaxonomyMatrix[, "TaxonNo"]), interval = Interval))
    
    # Deal with subgenera:
    resolvedtaxonnumbersInterval[, "TaxonName"] <- gsub(" \\(|\\)", "", resolvedtaxonnumbersInterval[, "TaxonName"])
    
    # Invert variable so only includes taxa outside Interval:
    resolvedtaxonnumbersInterval <- resolvedtaxonnumbers[is.na(resolvedtaxonnumbersInterval[, "TaxonName"]), ]
    
    # Find any nomen dubia etc. to delete:
    deleterows <- which(unlist(lapply(lapply(lapply(as.list(resolvedtaxonnumbersInterval[, "TaxonValidity"]), match, x = deletes), sort), length)) > 0)
    
    # If there are deletes then remove them from the matrix:
    if(length(deleterows) > 0) resolvedtaxonnumbersInterval <- resolvedtaxonnumbersInterval[-deleterows, , drop = FALSE]
    
  }
  
  # Print current processing status:
  cat("Done\nChecking taxon names match with database version...")
  
  # Vector to store any failed matches:
  failedmatches <- vector(mode = "character")
  
  # For each initially reconciled name:
  for(i in 1:nrow(TaxonomyMatrix)) {
    
    # Get resolved name:
    resolvedname <- gsub(" ", "_", resolvedtaxonnumbers[which(resolvedtaxonnumbers[, "InputNo"] == TaxonomyMatrix[i, "TaxonNo"]), "TaxonName"])
    
    # Get input name:
    inputname <- TaxonomyMatrix[i, "TaxonName"]
    
    # Check names truly match (i.e., deals with case of indets where direct match not possible) and store if not:
    if(resolvedname != inputname && any(is.na(match(strsplit(resolvedname, "_")[[1]], strsplit(inputname, "_")[[1]])))) failedmatches <- c(failedmatches, paste("Input name ", inputname, " does not match database name for corresponding number (", resolvedname, ").", sep = ""))
    
  }
  
  # If there are failed matches:
  if(length(failedmatches) > 0) {
    
    # Provide list to user:
    cat(paste(failedmatches, collapse = "\n"))
    
    # Stop function (will break later otherwise):
    stop("")
    
  }
  
  # Print current processing status:
  cat("Done\nChecking taxon validities...")
  
  # Check for any new kind of resolution (should be empty vector):
  newresolutions <- setdiff(sort(unique(resolvedtaxonnumbers[, "TaxonValidity"])), c(deletes, synonyms))
  
  # Stop if new resolutiosn found (need to add these to the resolution types above):
  if(length(newresolutions) > 0) stop(paste("New resolution type found!: ", newresolutions, sep = ""))
  
  # Print current processing status:
  cat("Done\nBuilding synonymy tables...")
  
  # Empty vector to store rows that correspond to some form of junior synonym:
  synonymrows <- c()
  
  # Find all junior synonym rows:
  for(i in synonyms) synonymrows <- sort(c(synonymrows, which(resolvedtaxonnumbers[, "TaxonValidity"] == i)))
  
  # Set junior synonym matrix:
  juniorsynonyms <- resolvedtaxonnumbers[synonymrows, ]
  
  # Create empty matrix to store senior synoyms:
  seniorsynonyms <- matrix(nrow = 0, ncol = 8, dimnames = list(c(), c("OriginalTaxonNo", "ResolvedTaxonNo", "TaxonName", "TaxonRank", "ParentTaxonNo", "TaxonValidity", "AcceptedNumber", "AcceptedName")))
  
  # Reconcile senior synonym with database:
  currenttaxa <- PaleobiologyDBTaxaQuerier(gsub("txn:", "", resolvedtaxonnumbers[synonymrows, "AcceptedNumber"]))
  
  # Deal with subgenera:
  currenttaxa[, "TaxonName"] <- gsub(" \\(|\\)", "", currenttaxa[, "TaxonName"])
  
  # While there are taxa with validity issues:
  while(any(!is.na(currenttaxa[, "TaxonValidity"]))) {
    
    # Isolate rows to check (i.e., just rows where validity isn't NA (i.e., resolved):
    rowstocheck <- which(!is.na(currenttaxa[, "TaxonValidity"]))
    
    # Check just those taxa:
    currenttaxa[rowstocheck, ] <- PaleobiologyDBTaxaQuerier(taxon_nos = gsub("txn:", "", currenttaxa[rowstocheck, "AcceptedNumber"]), taxon_names = currenttaxa[rowstocheck, "AcceptedName"])
    
    # Deal with subgenera:
    currenttaxa[rowstocheck, "TaxonName"] <- gsub(" \\(|\\)", "", currenttaxa[rowstocheck, "TaxonName"])
    
  }
  
  # Make current taxa into senior synonyms list:
  seniorsynonyms <- currenttaxa
  
  # If using an Interval:
  if(!all(is.null(Interval))) {
    
    # Update resolved taxon numbers to valid taxa only:
    resolvedtaxonnumbersInterval[which(!is.na(match(resolvedtaxonnumbersInterval[, "InputNo"], juniorsynonyms[, "InputNo"]))), c("OriginalTaxonNo", "ResolvedTaxonNo", "TaxonName", "TaxonRank", "ParentTaxonNo", "TaxonValidity", "AcceptedNumber", "AcceptedName")] <- seniorsynonyms[match(resolvedtaxonnumbersInterval[, "InputNo"], juniorsynonyms[, "InputNo"])[!is.na(match(resolvedtaxonnumbersInterval[, "InputNo"], juniorsynonyms[, "InputNo"]))], c("OriginalTaxonNo", "ResolvedTaxonNo", "TaxonName", "TaxonRank", "ParentTaxonNo", "TaxonValidity", "AcceptedNumber", "AcceptedName")]
    
  }
  
  # Print current processing status:
  cat("Done\nChecking validity of indeterminate taxon reconciliations...")
  
  # Get list of indeterminates:
  indeterminates <- TaxonomyMatrix[which((unlist(lapply(lapply(strsplit(TaxonomyMatrix[, "TaxonName"], "_"), '==', "aff"), sum)) + unlist(lapply(lapply(strsplit(TaxonomyMatrix[, "TaxonName"], "_"), '==', "cf"), sum)) + unlist(lapply(lapply(strsplit(TaxonomyMatrix[, "TaxonName"], "_"), '==', "indet"), sum)) + unlist(lapply(lapply(strsplit(TaxonomyMatrix[, "TaxonName"], "_"), '==', "sp"), sum))) > 0), "TaxonName"]
  
  # For each indeterminate:
  for(i in indeterminates) {
    
    # Get resolved row number:
    resolvedrownumber <- which(resolvedtaxonnumbers[, "InputNo"] == TaxonomyMatrix[which(TaxonomyMatrix[, "TaxonName"] == i), "TaxonNo"])
    
    # If a possible invalid taxon (validity is not blank):
    if(!is.na(resolvedtaxonnumbers[resolvedrownumber, "TaxonValidity"])) {
      
      # Get accepted number of taxon (may be NA):
      AcceptedNumber <- gsub("txn:", "", resolvedtaxonnumbers[resolvedrownumber, "AcceptedNumber"])
      
      # If accepted number is blank (NA) stop adn warn user taxon is invalid:
      if(is.na(AcceptedNumber)) stop(paste(i, " assigned to a taxon that is invalid, consider renaming.", sep = ""))
      
      # If accepted numebr is different to input number stop and warn user taxon is synonymised:
      if(AcceptedNumber != resolvedtaxonnumbers[resolvedrownumber, "InputNo"]) stop(paste(i, " assigned to a taxon that is invalid, consider renaming.", sep = ""))
      
    }
    
  }
  
  # Print current processing status:
  cat("Done\nDeleting taxa resolved as nomen dubium and the like...")
  
  # Empty vector to store rows that correspond to some form of junior synonym:
  deleterows <- deleterowsInterval <- c()
  
  # Find all junior synonym rows:
  for(i in deletes) deleterows <- sort(c(deleterows, which(resolvedtaxonnumbers[, "TaxonValidity"] == i)))
  
  # Get input numbers that should be deleted:
  numberstodelete <- resolvedtaxonnumbers[deleterows, "InputNo"]
  
  # Create list of taxon numbers to check for data sets with DELETE resolutions:
  taxonnumberslist <- lapply(lapply(lapply(lapply(lapply(lapply(lapply(lapply(MRPList[which(unlist(lapply(lapply(lapply(MRPList, '[[', "matrix"), rownames), length)) > 0)], '[[', 1), rownames), strsplit, split = "%%%%"), unlist), matrix, ncol = 2, byrow = TRUE), '[', ,1), strsplit, split = "&"), unlist)
  
  # Create empty vector to store datasets with deletes (that cna then be processed):
  datasetswithdeletes <- vector(mode = "character")
  
  # For each delete number find data sets that have it and add them to the vector:
  for(i in numberstodelete) datasetswithdeletes <- sort(unique(c(datasetswithdeletes, names(which(unlist(lapply(lapply(taxonnumberslist, '==', i), sum)) > 0)))))
  
  # For each dataset with at least one taxon to delete:
  for(i in datasetswithdeletes) {
    
    # Isolate taxon numbers for current dataset:
    taxonnumbers <- unlist(lapply(strsplit(rownames(MRPList[[i]]$matrix), "%%%%"), '[', 1))
    
    # Convert into list to deal with multi-number higher taxa:
    taxonnumbers <- lapply(lapply(as.list(taxonnumbers), strsplit, split = "&"), unlist)
    
    # Create empty vector to store taxa to delete:
    taxatodelete <- vector(mode = "numeric")
    
    # Find any taxa to delete present in the data set:
    for(j in numberstodelete) taxatodelete <- sort(unique(c(taxatodelete, which(unlist(lapply(lapply(taxonnumbers, '==', j), any))))))
    
    # Find any multi-number higher taxa that are listed as deletes:
    multinumberhighertaxatodelete <- taxatodelete[which(unlist(lapply(taxonnumbers[taxatodelete], length)) > 1)]
    
    # Stop if you find a multi-number higher taxon that is (at least partially) to be deleted):
    if(length(multinumberhighertaxatodelete) > 0) stop(paste("Multi-number higher tax(a) listed to delete: ", paste(unlist(lapply(strsplit(rownames(MRPList[[i]]$matrix), "%%%%"), '[', 2))[multinumberhighertaxatodelete], sep = ","), " (NB: check the status of these).", sep = ""))
    
    # Create MRP matrix:
    MRPmatrix <- MakeMorphMatrix(MRPList[[i]]$matrix[-taxatodelete, , drop = FALSE], header = "", weights = rep(1, ncol(MRPList[[i]]$matrix[-taxatodelete, , drop = FALSE])), ordering = rep("unord", ncol(MRPList[[i]]$matrix[-taxatodelete, , drop = FALSE])), equalise.weights = FALSE)
    
    # Overwrite matrix with collapsed version following taxon deletions:
    MRPList[[i]]$matrix <- MRPCollapse(MRPmatrix)$Matrix_1$Matrix
    
  }
  
  # Print current processing status:
  cat("Done\nReplacing junior synonyms with senior synonyms...")
  
  # Get input numbers that should be replaced:
  numberstosynonymise <- juniorsynonyms[, "InputNo"]
  
  # Create empty vector to store datasets with deletes (that cna then be processed):
  datasetswithjuniorsynonyms <- vector(mode = "character")
  
  # For each delete number find data sets that have it and add them to the vector:
  for(i in numberstosynonymise) datasetswithjuniorsynonyms <- sort(unique(c(datasetswithjuniorsynonyms, names(which(unlist(lapply(lapply(taxonnumberslist, '==', i), sum)) > 0)))))
  
  # For each dataset with at least one taxon to delete:
  for(i in datasetswithjuniorsynonyms) {
    
    # Isolate taxon numbers for current dataset:
    taxonnumbers <- unlist(lapply(strsplit(rownames(MRPList[[i]]$matrix), "%%%%"), '[', 1))
    
    # Convert into list to deal with multi-number higher taxa:
    taxonnumbers <- lapply(lapply(as.list(taxonnumbers), strsplit, split = "&"), unlist)
    
    # Create empty vector to store taxa to synonymise:
    taxatosynonymise <- vector(mode = "numeric")
    
    # Find any taxa to synonymise present in the data set:
    for(j in numberstosynonymise) taxatosynonymise <- sort(unique(c(taxatosynonymise, which(unlist(lapply(lapply(taxonnumbers, '==', j), any))))))
    
    # Find any multi-number higher taxa to synonymise:
    multinumberhighertaxontosynonymise <- unlist(lapply(strsplit(rownames(MRPList[[i]]$matrix), "%%%%"), '[', 2))[taxatosynonymise[which(unlist(lapply(taxonnumbers[taxatosynonymise], length)) > 1)]]
    
    # Stop if you find a multi-number higher taxon that is (at least partially) to be deleted):
    if(length(multinumberhighertaxontosynonymise) > 0) stop(paste("Multi-number higher tax(a) listed to synonymise: ", paste(multinumberhighertaxontosynonymise, sep = ","), " (NB: check the status of these).", sep = ""))
    
    # For each taxon to synonymise:
    for(j in taxatosynonymise) {
      
      # Find corresponding row number in senior synonyms list:
      seniorrownumber <- which(juniorsynonyms[, "InputNo"] == strsplit(rownames(MRPList[[i]]$matrix)[j], "%%%%")[[1]][1])
      
      # Overwrite junior with senior synonym:
      rownames(MRPList[[i]]$matrix)[j] <- paste(gsub("txn:|var:", "", rev(sort(seniorsynonyms[seniorrownumber, c("OriginalTaxonNo", "ResolvedTaxonNo")]))[1]), gsub(" ", "_", seniorsynonyms[seniorrownumber, "TaxonName"]), sep = "%%%%")
      
    }
    
  }
  
  # Print current processing status:
  cat("Done\nBuilding taxonomy...")
  
  # Get a list of the valid OTU names (may be pruned down later to just those in target clade):
  ValidOTUNames <- unique(unlist(lapply(lapply(MRPList, '[[', "matrix"), rownames)))[grep("_", unique(unlist(lapply(lapply(MRPList, '[[', "matrix"), rownames))))]
  
  # Replace junior with senior synonyms in resolved names matrix:
  resolvedtaxonnumbers[synonymrows, c("OriginalTaxonNo", "ResolvedTaxonNo", "TaxonName", "TaxonRank", "ParentTaxonNo", "TaxonValidity", "AcceptedNumber", "AcceptedName")] <- seniorsynonyms
  
  # Overwrite resolved number with resolved taxon number:
  resolvedtaxonnumbers[, "ResolvedTaxonNo"] <- gsub("txn:|var:", "", unlist(lapply(lapply(apply(resolvedtaxonnumbers[, c("OriginalTaxonNo", "ResolvedTaxonNo")], 1, sort), rev), '[', 1)))
  
  # Remove original taxon number column:
  resolvedtaxonnumbers <- resolvedtaxonnumbers[, -which(colnames(resolvedtaxonnumbers) == "OriginalTaxonNo")]
  
  # Remove deleted taxa from resolved names matrix (if there are any):
  if(length(which(!is.na(resolvedtaxonnumbers[, "TaxonValidity"]))) > 0) resolvedtaxonnumbers <- resolvedtaxonnumbers[-which(!is.na(resolvedtaxonnumbers[, "TaxonValidity"])), ]
  
  # Reformat parent taxon numbers into just numbers:
  resolvedtaxonnumbers[, "ParentTaxonNo"] <- gsub("txn:", "", resolvedtaxonnumbers[, "ParentTaxonNo"])
  
  # Collapse resolved matrix to just field with values (i.e., drop valid and senior synonym columns):
  resolvedtaxonnumbers <- resolvedtaxonnumbers[, c("ResolvedTaxonNo", "TaxonName", "TaxonRank", "ParentTaxonNo")]
  
  # If doing something with missing species (i.e., those not currently included as OTUs, but existing in target clade):
  if(MissingSpecies != "exclude") {
    
    # Find all children of target clade:
    AllChildren <- PaleobiologyDBChildFinder(taxon_nos = "1", taxon_names = TargetClade, validonly = TRUE, returnrank = "3", interval = Interval)
    
    # Deal with subgenera:
    AllChildren[, "TaxonName"] <- gsub(" \\(|\\)", "", AllChildren[, "TaxonName"])
    
    # If inserting all missing species get all possible species parent numbers:
    if(MissingSpecies == "all") CurrentSpeciesParentNumbers <- unique(c(gsub("txn:", "", AllChildren[, "ParentTaxonNo"]), resolvedtaxonnumbers[which(resolvedtaxonnumbers[, "TaxonRank"] == 3), "ParentTaxonNo"]))
    
    # If only inserting missing species at genus-level find parent numbers of all current species (i.e., potential genera to add):
    if(MissingSpecies == "genus") CurrentSpeciesParentNumbers <- unique(resolvedtaxonnumbers[which(resolvedtaxonnumbers[, "TaxonRank"] == 3), "ParentTaxonNo"])
    
    # Find any parents not already present in resolved numbers matrix:
    AsYetUnsampledSpeciesParents <- setdiff(CurrentSpeciesParentNumbers, resolvedtaxonnumbers[, "ResolvedTaxonNo"])
    
    # If such parents exist:
    if(length(AsYetUnsampledSpeciesParents) > 0) {
      
      # Find unsampled species parents:
      CurrentSpeciesParents <- PaleobiologyDBTaxaQuerier(taxon_nos = AsYetUnsampledSpeciesParents)
      
      # Deal with subgenera:
      CurrentSpeciesParents[, "TaxonName"] <- gsub(" \\(|\\)", "", CurrentSpeciesParents[, "TaxonName"])
      
      # Find rows corresponding to valid genera:
      ValidGenusRows <- intersect(which(is.na(CurrentSpeciesParents[, "TaxonValidity"])), which(CurrentSpeciesParents[, "TaxonRank"] == "5"))
      
      # If there are valid genera then add these to resolved taxon numbers:
      if(length(ValidGenusRows) > 0) resolvedtaxonnumbers <- rbind(resolvedtaxonnumbers, cbind(unname(gsub("txn:|var:", "", unlist(lapply(lapply(lapply(apply(CurrentSpeciesParents[ValidGenusRows, c("OriginalTaxonNo", "ResolvedTaxonNo"), drop = FALSE], 1, list), unlist), sort, decreasing = TRUE), '[', 1)))), CurrentSpeciesParents[ValidGenusRows, c("TaxonName", "TaxonRank")] , gsub("txn:", "", CurrentSpeciesParents[ValidGenusRows, "ParentTaxonNo"])))
      
    }
    
    # If including all species:
    if(MissingSpecies == "all") {
      
      # Update Valid OTUs accordingly
      ValidOTUNames <- unique(c(ValidOTUNames, gsub(" ", "_", paste(unlist(lapply(strsplit(gsub("NA", "", paste(AllChildren[, "OriginalTaxonNo"], AllChildren[, "ResolvedTaxonNo"], sep = "")), split = "var:|txn:"), '[[', 2)), AllChildren[, "TaxonName"], sep = "%%%%"))))
      
      # Set new children to add to resolved taxon numbers later:
      NewChildren <- matrix(c(gsub("txn:|var:", "", unlist(lapply(apply(AllChildren[, c("OriginalTaxonNo", "ResolvedTaxonNo")], 1, sort, decreasing = TRUE), '[', 1))), AllChildren[, c("TaxonName", "TaxonRank")], gsub("txn:", "", AllChildren[, "ParentTaxonNo"])), ncol = 4)
      
    }
    
    # If including only species assigned to genus-level OTUs:
    if(MissingSpecies == "genus") {
      
      # Get current genus numbers (to check what has already been included):
      CurrentGenusNumbers <- resolvedtaxonnumbers[which(resolvedtaxonnumbers[, "TaxonRank"] == 5), "ResolvedTaxonNo"]
      
      # Get children of sampled genera:
      GeneraChildren <- PaleobiologyDBChildFinder(taxon_nos = CurrentGenusNumbers, validonly = TRUE, returnrank = "3")
      
      # Deal with subgenera:
      GeneraChildren[, "TaxonName"] <- gsub(" \\(|\\)", "", GeneraChildren[, "TaxonName"])
      
      # Update valid OTUs with children of all sampled genera:
      ValidOTUNames <- unique(c(ValidOTUNames, gsub(" ", "_", paste(unlist(lapply(strsplit(gsub("NA", "", paste(GeneraChildren[, "OriginalTaxonNo"], GeneraChildren[, "ResolvedTaxonNo"], sep = "")), split = "var:|txn:"), '[[', 2)), GeneraChildren[, "TaxonName"], sep = "%%%%"))))
      
      # Set new children to add to resolved taxon numbers later:
      NewChildren <- matrix(c(gsub("txn:|var:", "", unlist(lapply(apply(GeneraChildren[, c("OriginalTaxonNo", "ResolvedTaxonNo")], 1, sort, decreasing = TRUE), '[', 1))), GeneraChildren[, c("TaxonName", "TaxonRank")], gsub("txn:", "", GeneraChildren[, "ParentTaxonNo"])), ncol = 4)
      
    }
    
    # Find any new children not already included in resolved taxon numbers list:
    ChildrenToAdd <- setdiff(NewChildren[, 1], resolvedtaxonnumbers[, "ResolvedTaxonNo"])
    
    # If there are children to add then add them to resolved taxon numbers:
    if(length(ChildrenToAdd) > 0) resolvedtaxonnumbers <- rbind(resolvedtaxonnumbers, NewChildren[match(ChildrenToAdd, NewChildren[, 1]), ])
    
  }
  
  # Get initial parent child relationships based on OTUs:
  parentchildrelationships <- paste(unlist(lapply(strsplit(ValidOTUNames, "%%%%"), '[', 1)), resolvedtaxonnumbers[match(unlist(lapply(strsplit(ValidOTUNames, "%%%%"), '[', 1)), resolvedtaxonnumbers[, "ResolvedTaxonNo"]), "ParentTaxonNo"], sep = " belongs to ")
  
  # Find which rows correspond to indeterminate and sp taxa (i.e., those where parent should be initial reconciliation):
  indetsandsps <- sort(c(which(unlist(lapply(lapply(lapply(lapply(lapply(strsplit(ValidOTUNames, "%%%%"), '[', 2), strsplit, split = "_"), unlist), '==', "indet"), any))), which(unlist(lapply(lapply(lapply(lapply(lapply(strsplit(ValidOTUNames, "%%%%"), '[', 2), strsplit, split = "_"), unlist), '==', "sp"), any)))))
  
  # If such taxa exist then update parent child relationships accordingly:
  if(length(indetsandsps) > 0) parentchildrelationships[indetsandsps] <- paste(unlist(lapply(strsplit(ValidOTUNames[indetsandsps], "%%%%"), '[', 1)), unlist(lapply(strsplit(ValidOTUNames[indetsandsps], "%%%%"), '[', 1)), sep = " belongs to ")
  
  # Get list of new children (for which parents are needed) - excludes "Life" which has no parent:
  newchildren <- setdiff(unlist(lapply(strsplit(parentchildrelationships, " belongs to "), '[', 2)), "28595")
  
  # As long as there are still children in need of parents:
  while(length(newchildren) > 0) {
    
    # Find any numbers missing for the taxonomy name resolution matrix:
    missingfromresolutions <- newchildren[which(is.na(match(newchildren, resolvedtaxonnumbers[, "ResolvedTaxonNo"])))]
    
    # If there are such numbers:
    if(length(missingfromresolutions) > 0) {
      
      # Get raw query data for new names
      rawquery <- PaleobiologyDBTaxaQuerier(taxon_nos = missingfromresolutions)
      
      # Deal with subgenera:
      rawquery[, "TaxonName"] <- gsub(" \\(|\\)", "", rawquery[, "TaxonName"])
      
      # Add formatted results of query to resolved names matrix:
      resolvedtaxonnumbers <- rbind(resolvedtaxonnumbers, cbind(gsub("txn:|var:", "", unname(unlist(lapply(lapply(lapply(apply(rawquery[, c("OriginalTaxonNo", "ResolvedTaxonNo"), drop = FALSE], 1, list), unlist), sort, decreasing = TRUE), '[', 1)))), rawquery[, c("TaxonName", "TaxonRank"), drop = FALSE], gsub("txn:", "", rawquery[, "ParentTaxonNo"])))
      
    }
    
    # Add new parent child relationships to list:
    parentchildrelationships <- c(parentchildrelationships, paste(resolvedtaxonnumbers[match(newchildren, resolvedtaxonnumbers[, "ResolvedTaxonNo"]), "ResolvedTaxonNo"], resolvedtaxonnumbers[match(newchildren, resolvedtaxonnumbers[, "ResolvedTaxonNo"]), "ParentTaxonNo"], sep = " belongs to "))
    
    # Update new children:
    newchildren <- setdiff(resolvedtaxonnumbers[match(newchildren, resolvedtaxonnumbers[, "ResolvedTaxonNo"]), "ParentTaxonNo"], "28595")
    
  }
  
  # If Life is missing then add it at bottom:
  if(all(!resolvedtaxonnumbers[, "ResolvedTaxonNo"] == "28595")) resolvedtaxonnumbers <- rbind(resolvedtaxonnumbers, c("28595", "Life", "25", NA))
  
  # Convert parent-child relationships into a matrix (columns for child and parent):
  parentchildmatrix <- matrix(unlist(strsplit(parentchildrelationships, split = " belongs to ")), ncol = 2, byrow = TRUE, dimnames = list(c(), c("Child", "Parent")))
  
  # Update parent-child matrix with child names:
  parentchildmatrix[, "Child"] <- resolvedtaxonnumbers[match(parentchildmatrix[, "Child"], resolvedtaxonnumbers[, "ResolvedTaxonNo"]), "TaxonName"]
  
  # Update parent-child matrix with parent names:
  parentchildmatrix[, "Parent"] <- resolvedtaxonnumbers[match(parentchildmatrix[, "Parent"], resolvedtaxonnumbers[, "ResolvedTaxonNo"]), "TaxonName"]
  
  # Add valid OTU names into parent-child matrix:
  parentchildmatrix[c(1:length(ValidOTUNames)), "Child"] <- unlist(lapply(strsplit(ValidOTUNames, "%%%%"), '[', 2))
  
  # Add missing taxon ("Life") to parent-child matrix:
  parentchildmatrix[which(is.na(parentchildmatrix[, "Parent"])), "Parent"] <- "Life"
  
  # Create empty taxonomy MRP matrix:
  TaxonomyMRP <- matrix(0, nrow = length(ValidOTUNames), ncol = length(sort(unique(parentchildmatrix[, "Parent"]))), dimnames = list(sort(unlist(lapply(strsplit(ValidOTUNames, "%%%%"), '[', 2))), sort(unique(parentchildmatrix[, "Parent"]))))
  
  # Remove duplicates:
  parentchildmatrix <- matrix(unlist(strsplit(unique(paste(parentchildmatrix[, "Child"], parentchildmatrix[, "Parent"], sep = "%%%%")), "%%%%")), ncol = 2, byrow = TRUE, dimnames = list(c(), c("Child", "Parent")))
  
  # Check for duplicate names:
  if(any(duplicated(sort(unlist(lapply(strsplit(ValidOTUNames, "%%%%"), '[[', 2)))))) stop(paste("The following OTU names are duplicated in the database (check and correct): ", paste(sort(unlist(lapply(strsplit(ValidOTUNames, "%%%%"), '[[', 2)))[duplicated(sort(unlist(lapply(strsplit(ValidOTUNames, "%%%%"), '[[', 2))))], sep = ", "), sep = ""))
  
  # For each OTU (traces up through hierarchy until its presence in every higher taxon to which it belongs is assigned)::
  for(i in 1:length(ValidOTUNames)) {
    
    # Set starting current child taxon:
    currentchild <- parentchildmatrix[i, "Child"]
    
    # Set starting current parent taxon:
    currentparent <- setdiff(parentchildmatrix[which(parentchildmatrix[, "Child"] == currentchild), "Parent"], currentchild)
    
    # Record presence of child in parent in taxonomy matrix:
    TaxonomyMRP[parentchildmatrix[i, "Child"], currentparent] <- 1
    
    # As long as the parent is not "Life" (top of taxonomic hierarchy not reached):
    while(currentparent != "Life") {
      
      # Update child with previous parent:
      currentchild <- currentparent
      
      # Update parent with new parent:
      currentparent <- setdiff(parentchildmatrix[which(parentchildmatrix[, "Child"] == currentchild), "Parent"], currentchild)
      
      # Record presence of child in parent:
      TaxonomyMRP[parentchildmatrix[i, "Child"], currentparent] <- 1
      
    }
    
  }
  
  # Print current processing status:
  cat("Done\nDealing with subgenera...")
  
  # Find any subgenera names (as supraspecific column names only):
  subgenerarows <- grep("\\(", colnames(TaxonomyMRP))
  
  # If subgenera found make these single names (i.e., removes parentheses that will screw up Newick trees later):
  if(length(subgenerarows) > 0) colnames(TaxonomyMRP)[subgenerarows] <- gsub("\\(|\\)| ", "", colnames(TaxonomyMRP)[subgenerarows])
  
  # Correct subgenera in MRP list too:
  MRPList <- lapply(MRPList, function(x) {NamesToCheck <- rownames(x$matrix); NamesToCheck <- gsub("_\\(|\\)", "", NamesToCheck); rownames(x$matrix) <- NamesToCheck; x})
  
  # Print current processing status:
  cat("Done\nTidying up taxonomy...")
  
  # Check target clade is actually found:
  if(length(which(colnames(TaxonomyMRP) == TargetClade)) == 0) stop("Target clade not found in taxonomy. Check spelling/Paleobiology database validity.")
  
  # Work out which taxa are actually valid OTUs (belong to target clade):
  NewValidOTUs <- names(which(TaxonomyMRP[, TargetClade] == 1))
  
  # Find any subspecies names:
  SubspeciesNames <- NewValidOTUs[unlist(lapply(strsplit(NewValidOTUs, split = ""), function(x) all(c(sum(x == "_") == 2, length(grep("[:A-Z:]", x)) == 1))))]
  
  # Only continue if subspecies were found:
  if(length(SubspeciesNames) > 0) {
    
    # Find any subspecues where species is not in sample:
    SubspeciesWhereSpeciesIsNotFound <- SubspeciesNames[!unlist(lapply(as.list(SubspeciesNames), function(x) any(NewValidOTUs == paste(strsplit(x, split = "_")[[1]][1:2], collapse = "_"))))]
    
    # If subspecies without sampled species :
    if(length(SubspeciesWhereSpeciesIsNotFound) > 0) {
      
      # Rename these with species names:
      NewValidOTUs[match(SubspeciesWhereSpeciesIsNotFound, NewValidOTUs)] <- unlist(lapply(strsplit(SubspeciesWhereSpeciesIsNotFound, split = "_"), function(x) paste(x[1:2], collapse = "_")))
      
      # Rename taxonomy MRP rownames with species names too:
      rownames(TaxonomyMRP)[match(SubspeciesWhereSpeciesIsNotFound, rownames(TaxonomyMRP))] <- unlist(lapply(strsplit(SubspeciesWhereSpeciesIsNotFound, split = "_"), function(x) paste(x[1:2], collapse = "_")))
      
    }
    
    # Collapse subspecies back to just the names where the species is already sampled:
    SubspeciesNames <- setdiff(SubspeciesNames, SubspeciesWhereSpeciesIsNotFound)
    
    # If there are subspecies where species already exists then prune these from the taxonomy MRP:
    if(length(SubspeciesNames) > 0) NewValidOTUs <- NewValidOTUs[-match(SubspeciesNames, NewValidOTUs)]
    
  }
  
  # Modify this if using intervals:
  if(!all(is.null(Interval))) NewValidOTUs <- setdiff(NewValidOTUs, gsub(" ", "_", resolvedtaxonnumbersInterval[, "TaxonName"]))
  
  # Can now strip out numbers from taxon names:
  for(i in 1:length(MRPList)) if(!is.null(rownames(MRPList[[i]]$matrix))) rownames(MRPList[[i]]$matrix) <- unlist(lapply(strsplit(rownames(MRPList[[i]]$matrix), "%%%%"), '[', 2))
  
  # Collapse taxonomy MRP to just new valid taxa:
  TaxonomyMRP <- TaxonomyMRP[NewValidOTUs, ]
  
  # Make taxonomy MRP into list:
  TaxonomyMRPlist <- split(TaxonomyMRP, rep(1:ncol(TaxonomyMRP), each = nrow(TaxonomyMRP)))
  
  # Add column names to list:
  names(TaxonomyMRPlist) <- colnames(TaxonomyMRP)
  
  # Find higher taxa for which every taxon is present:
  redundanthighertaxa <- colnames(TaxonomyMRP)[intersect(which(unlist(lapply(lapply(TaxonomyMRPlist, unique), length)) == 1), which(unlist(lapply(lapply(TaxonomyMRPlist, unique), '[', 1)) == 1))]
  
  # Empty higher taxa:
  emptyhighertaxa <- colnames(TaxonomyMRP)[intersect(which(unlist(lapply(lapply(TaxonomyMRPlist, unique), length)) == 1), which(unlist(lapply(lapply(TaxonomyMRPlist, unique), '[', 1)) == 0))]
  
  # Find taxonomic autapomorphies (those with just one OTU and hence redundant):
  taxonomicautapomorphies <- names(which(unlist(lapply(TaxonomyMRPlist, sum)) == 1))
  
  # Collapse taxonomy MRP by removing constant characters (i.e., most of the subgroups just established - not autapomorphies as they can be substitutes later!):
  TaxonomyMRP <- TaxonomyMRP[, -match(c(redundanthighertaxa, emptyhighertaxa), colnames(TaxonomyMRP)), drop = FALSE]
  
  # Print current processing status:
  cat("Done\nSubstituting valid OTUs for supraspecific taxa...")
  
  # Find datasets with (valid) surpaspecific OTUs:
  datasetswithsupraspecificOTUs <- which(unlist(lapply(lapply(lapply(lapply(MRPList, '[[', "matrix"), rownames), intersect, y = colnames(TaxonomyMRP)), length)) > 0)
  
  # If such data sets exist:
  if(length(datasetswithsupraspecificOTUs)) {
    
    # For each such data set:
    for(i in datasetswithsupraspecificOTUs) {
      
      # Find higher taxa that will need to be replaced:
      highertaxatoreplace <- intersect(rownames(MRPList[[i]]$matrix), colnames(TaxonomyMRP))
      
      # For each higher taxon to replace:
      for(j in highertaxatoreplace) {
        
        # Find substitue names from taxonomy:
        substitutenames <- names(which(TaxonomyMRP[, j] == 1))
        
        # Add these to end of matrix using coding for higher taxon:
        MRPList[[i]]$matrix <- rbind(MRPList[[i]]$matrix, matrix(rep(MRPList[[i]]$matrix[j, ], length(substitutenames)), nrow = length(substitutenames), byrow = TRUE, dimnames = list(substitutenames, c())))
        
        # Remove now replaced higher taxon from matrix:
        MRPList[[i]]$matrix <- MRPList[[i]]$matrix[-which(rownames(MRPList[[i]]$matrix) == j), , drop = FALSE]
        
      }
      
    }
    
  }
  
  # Print current processing status:
  cat("Done\nRetracting subspecies into species...")
  
  # Replace all subspecies with their species name:
  MRPList <- lapply(MRPList, function(x) {UnderscoreAndCapitalCounts <- matrix(unlist(lapply(strsplit(rownames(x$matrix), split = ""), function(y) c(sum(y == "_"), length(grep("[:A-Z:]", y))))), ncol = 2, byrow = TRUE, dimnames = list(c(), c("Underscores", "Capitals"))); SubspeciesRows <- intersect(which(UnderscoreAndCapitalCounts[, "Underscores"] == 2), which(UnderscoreAndCapitalCounts[, "Capitals"] == 1)); if(length(SubspeciesRows) > 0) rownames(x$matrix)[SubspeciesRows] <- unlist(lapply(strsplit(rownames(x$matrix)[SubspeciesRows], split = "_"), function(z) paste(z[1:2], collapse = "_"))); x})
  
  # Find datatsets with duplicated OTU names:
  DuplicatedOTUDatasets <- which(unlist(lapply(MRPList, function(x) any(duplicated(rownames(x$matrix))))))
  
  # If duplicated OTU names then collapse these matrices to just unique taxa:
  if(length(DuplicatedOTUDatasets) > 0) MRPList[DuplicatedOTUDatasets] <- lapply(MRPList[DuplicatedOTUDatasets], function(x) {MorphMatrix <- list(list("", NULL), list(NA, "STANDARD", x$matrix, rep("unord", ncol(x$matrix)), rep(1, ncol(x$matrix)), rep(0, ncol(x$matrix)), rep(1, ncol(x$matrix)), list(c("0", "1"), "?", "-"))); names(MorphMatrix) <- c("Topper", "Matrix_1"); names(MorphMatrix[[1]]) <- c("Header", "StepMatrices"); names(MorphMatrix[[2]]) <- c("BlockName", "Datatype", "Matrix", "Ordering", "Weights", "MinVals", "MaxVals", "Characters"); names(MorphMatrix[[2]][[8]]) <- c("Symbols", "Missing", "Gap"); x$matrix <- CompactifyMatrix(MorphMatrix)$Matrix_1$Matrix; x})
  
  # Print current processing status:
  cat("Done\nFurther tidying of taxonomy...")
  
  # If applying an Interval:
  if(!all(is.null(Interval))) {
    
    # Find any rows to delete (because they represent taxa outside the Interval):
    RowsToDelete <- sort(match(gsub(" ", "_", resolvedtaxonnumbersInterval[resolvedtaxonnumbersInterval[, "TaxonRank"] == "3", "TaxonName"]), rownames(TaxonomyMRP)))
    
    # If found then remove these from the taxonomy MRP:
    if(length(RowsToDelete) > 0) TaxonomyMRP <- TaxonomyMRP[-RowsToDelete, , drop = FALSE]
    
  }
  
  # Collapse taxonomy MRP by removing autapomorphic characters (if any):
  if(length(taxonomicautapomorphies) > 0) TaxonomyMRP <- TaxonomyMRP[, -match(taxonomicautapomorphies, colnames(TaxonomyMRP)), drop = FALSE]
  
  # Overwrite taxonomy MRP list with strings for each column:
  TaxonomyMRPlist <- apply(TaxonomyMRP, 2, paste, collapse = "")
  
  # Find any duplicated MRP strings:
  duplicatedMRPstrings <- rle(sort(TaxonomyMRPlist))$values[which(rle(sort(TaxonomyMRPlist))$lengths > 1)]
  
  # If there are duplicated columns (i.e., redundant MRP characters in the taxonomy):
  if(length(duplicatedMRPstrings) > 0) {
    
    # For each duplicated character:
    for(i in duplicatedMRPstrings) {
      
      # Get duplicated columns for current MRP string:
      duplicatedcolumns <- which(TaxonomyMRPlist == i)
      
      # Form new column name by collapsing higher taxa that are duplicated to a single string:
      newcolumnname <- paste(names(duplicatedcolumns), collapse = "_et_")
      
      # Overwrite first duplicated column name with new collapsed name:
      colnames(TaxonomyMRP)[duplicatedcolumns[1]] <- newcolumnname
      
      # Remove redundant columns from matrix:
      TaxonomyMRP <- TaxonomyMRP[, -duplicatedcolumns[2:length(duplicatedcolumns)], drop = FALSE]
      
      # Remove redundant columns from list:
      TaxonomyMRPlist <- TaxonomyMRPlist[-duplicatedcolumns[2:length(duplicatedcolumns)]]
      
    }
    
  }
  
  # Print current processing status:
  cat("Done\nRemoving outgroups, empty supraspecific taxa and those outside the sampling interval...")
  
  # Find any remaining taxa that now need to be deleted (outgroups to target clade and empty higher taxa):
  TaxaToDelete <- setdiff(unlist(lapply(lapply(MRPList, '[[', "matrix"), rownames)), NewValidOTUs)
  
  # If applying an Interval then add taxa outside of it to the deletes list:
  if(!all(is.null(Interval))) TaxaToDelete <- unique(c(TaxaToDelete, gsub(" ", "_", resolvedtaxonnumbersInterval[resolvedtaxonnumbersInterval[, "TaxonRank"] == "3", "TaxonName"])))
  
  # Find data sets with taxa to delete:
  datasetswithtaxatodelete <- which(unlist(lapply(lapply(lapply(lapply(MRPList, '[[', "matrix"), rownames), intersect, y = TaxaToDelete), length)) > 0)
  
  # If there are such data sets:
  if(length(datasetswithtaxatodelete) > 0) {
    
    # For each data set with taxa to delete (done in reverse order in case deletions are found):
    for(i in rev(datasetswithtaxatodelete)) {
      
      # Collapse matrix to just valid OTUs:
      MRPList[[i]]$matrix <- MRPList[[i]]$matrix[intersect(rownames(MRPList[[i]]$matrix), NewValidOTUs), , drop = FALSE]
      
      # Get logical for how many characters are variable:
      VariableCharacters <- unlist(lapply(lapply(lapply(apply(matrix(apply(MRPList[[i]]$matrix, 2, paste, collapse = "")), 1, strsplit, split = ""), unlist), unique), length)) > 1
      
      # As long as there are variable characters:
      if(any(VariableCharacters)) {
        
        # If there are less than three taxa left collapse to a zero-column matrix:
        if(nrow(MRPList[[i]]$matrix) < 3) MRPList[[i]]$matrix <- MRPList[[i]]$matrix[, -c(1:ncol(MRPList[[i]]$matrix)), drop = FALSE]
        
        # If there are at least three taxa left (minimum for meaning) collapse to just meaningful:
        if(nrow(MRPList[[i]]$matrix) >= 3) MRPList[[i]]$matrix <- MRPCollapse(MakeMorphMatrix(MRPList[[i]]$matrix, header = "", weights = rep(1, ncol(MRPList[[i]]$matrix)), ordering = rep("unord", ncol(MRPList[[i]]$matrix)), equalise.weights = FALSE))$Matrix_1$Matrix
        
        # If there are no variable characters:
      } else {
        
        # Remove ith data set from MRP list (as has no meaningful characters):
        MRPList <- MRPList[-i]
        
      }
      
    }
    
  }
  
  # Print current processing status:
  cat("Done\nProducing taxonomy tree...")
  
  # Duplicated Txonomy MRP to create a collapsable version for generating taxonomy Newick string:
  TaxonomyMRPNewick <- TaxonomyMRP
  
  # Get order of columns to collapse to form MRP
  columncollapseorder <- order(apply(TaxonomyMRPNewick, 2, sum))
  
  # For each column ("clade") in order from samllest to largest:
  for(i in columncollapseorder) {
    
    # Get taxa present in current clade:
    taxonrows <- which(TaxonomyMRPNewick[, i] == 1)
    
    # Create new partial Newick string for current clade (node):
    newNewickstring <- paste("(", paste(names(taxonrows), collapse = ","), ")", colnames(TaxonomyMRPNewick)[i], sep = "")
    
    # Replace row name with new Newick string:
    rownames(TaxonomyMRPNewick)[taxonrows[1]] <- newNewickstring
    
    # Remove now redundant taxa from Newick matrix:
    TaxonomyMRPNewick <- TaxonomyMRPNewick[-taxonrows[2:length(taxonrows)], , drop = FALSE]
    
  }
  
  # Complete Newick string and over write Taxonomy Newick matrix:
  TaxonomyMRPNewick <- paste("(", paste(rownames(TaxonomyMRPNewick), collapse = ","), ")", TargetClade, ";", sep = "")
  
  # Ladderize taxonomy tree for neatness!:
  TaxonomyMRPTree <- ladderize(read.tree(text = TaxonomyMRPNewick))
  
  # Print current processing status:
  cat("Done\nAdding NAs for indet. and sp. subclades to taxonomy MRP...")
  
  # Find indeterminates and sps:
  indetsandsps <- NewValidOTUs[which((unlist(lapply(lapply(strsplit(NewValidOTUs, "_"), '==', "indet"), sum)) + unlist(lapply(lapply(strsplit(NewValidOTUs, "_"), '==', "sp"), sum))) > 0)]
  
  # If there are indeterminates and/or sps:
  if(length(indetsandsps) > 0) {
    
    # For each such taxon:
    for(i in indetsandsps) {
      
      # Find higher taxon to which it belongs:
      highertaxon <- colnames(TaxonomyMRP)[which(unlist(lapply(lapply(strsplit(colnames(TaxonomyMRP), "_et_"), '==', strsplit(i, "_")[[1]][1]), sum)) == 1)]
      
      # Find any sub (suprapseicifc taxa) for that higher taxon:
      subtaxa <- colnames(TaxonomyMRP)[which(unlist(lapply(lapply(lapply(split(TaxonomyMRP[names(which(TaxonomyMRP[, highertaxon] == 1)), , drop = FALSE], rep(1:ncol(TaxonomyMRP[names(which(TaxonomyMRP[, highertaxon] == 1)), , drop = FALSE]), each = nrow(TaxonomyMRP[names(which(TaxonomyMRP[, highertaxon] == 1)), , drop = FALSE]))), sort), unique), length)) == 2)]
      
      # If these exist then set ith taxon as being NA with respect to belonging to the subtax(a):
      if(length(subtaxa) > 1) TaxonomyMRP[i, subtaxa] <- NA
      
    }
    
  }
  
  # Print current processing status:
  cat("Done\nGetting weighting data (publication year and dependencies)...")
  
  # Get publication years for each data set:
  publicationyears <- as.numeric(gsub("[:a-z:]", "", gsub("inpress", strsplit(as.character(Sys.Date()), "-")[[1]][1], unlist(lapply(lapply(strsplit(names(MRPList), "_"), rev), '[', 1)))))
  
  # Get parent (data set) data:
  parentsdata <- unlist(lapply(MRPList, '[[', "parent"))
  
  # Get sibling (data set) data:
  siblingsdata <- unlist(lapply(MRPList, '[[', "sibling"))
  
  # Get multiple-offspring data sets:
  multioffspringdatasets <- unique(sort(parentsdata[which(nchar(parentsdata) > 0)])[which(duplicated(sort(parentsdata[which(nchar(parentsdata) > 0)])))])
  
  # If there are such data sets:
  if(length(multioffspringdatasets) > 0) {
    
    # Update siblings data accordingly:
    for(i in multioffspringdatasets) siblingsdata[which(parentsdata == i)] <- paste(names(which(parentsdata == i)), collapse = "%%%%")
    
  }
  
  # Create pool of parents present amongst data sets (children of them must also be present by definition):
  parentspresentinpool <- intersect(parentsdata, names(MRPList))
  
  # Create empty vector to store data sets to remove from pool:
  datasetstoremove <- vector(mode = "character")
  
  # For each parent in pool:
  for(i in parentspresentinpool) {
    
    # Case if parent made redundant:
    if(any(unlist(lapply(lapply(lapply(lapply(MRPList[which(parentsdata == i)], '[[', 1), rownames), setdiff, x = rownames(MRPList[[i]]$matrix)), length)) == 0)) {
      
      # Add data set to remove list:
      datasetstoremove <- c(datasetstoremove, i)
      
      # Case if parent not redundant (becomes sibling):
    } else {
      
      # Update siblings data with parent:
      siblingsdata[sort(c(i, names(MRPList[which(parentsdata == i)])))] <- paste(sort(c(names(MRPList[which(parentsdata == i)]), unlist(strsplit(siblingsdata[sort(c(i, names(MRPList[which(parentsdata == i)])))], "%%%%")))), collapse = "%%%%")
      
    }
    
  }
  
  # Expunge data sets to remove from siblings lists:
  siblingsdata <- unlist(lapply(lapply(lapply(lapply(siblingsdata, strsplit, split = "%%%%"), unlist), setdiff, y = datasetstoremove), paste, collapse = "%%%%"))
  
  # Remove now redundant data sets from siblings:
  siblingsdata <- siblingsdata[-match(datasetstoremove, names(MRPList))]
  
  # Remove sibling relationships of size one (effectively pointless!):
  if(length(intersect(which(nchar(siblingsdata) > 0), setdiff(c(1:length(siblingsdata)), grep("%%%%", siblingsdata)))) > 0) siblingsdata[intersect(which(nchar(siblingsdata) > 0), setdiff(c(1:length(siblingsdata)), grep("%%%%", siblingsdata)))] <- ""
  
  # Remove now redundant data sets from publication years (if there are any):
  if(length(datasetstoremove) > 0) publicationyears <- publicationyears[-match(datasetstoremove, names(MRPList))]
  
  # Remove now redundant data sets from MRP list (if there are any):
  if(length(datasetstoremove) > 0) MRPList <- MRPList[-match(datasetstoremove, names(MRPList))]
  
  # Print current processing status:
  cat("Done\nCalculating weights...")
  
  # Equation 1 in Supplementary Information of Lloyd et al. (2016):
  publicationyearweights <- 10 * 2^(0.5 * (publicationyears - min(publicationyears)))
  
  # Get indepenece weights (using siblings data):
  independenceweights <- 1 / apply(rbind(rep(1, length(MRPList)), unlist(lapply(lapply(lapply(siblingsdata, strsplit, split = "%%%%"), unlist), length))), 2, max)
  
  # Get N characters weight (to avoid ucnertain data sets swamping the signal):
  Ncharacterweights <- 1 / apply(rbind(rep(1, length(MRPList)), unlist(lapply(lapply(MRPList, '[[', "matrix"), ncol))), 2, max)
  
  # Combine weights (just a product for now, but other options should be explored!:
  publicationweights <- publicationyearweights * independenceweights * Ncharacterweights
  
  # Make range of weights span 990 (difference between 10 and 1000):
  publicationweights <- (1 / (diff(range(publicationweights)) / 990)) * publicationweights
  
  # Make minimum weight 10 (and by extension maximum weight 1000):
  publicationweights <- (10 - min(publicationweights)) + publicationweights
  
  # Restrict weights to 2dp:
  publicationweights <- round(publicationweights, 2)
  
  # Print current processing status:
  cat("Done\nBuilding MRP matrix...")
  
  # Create empty full MRP matrix:
  FullMRPMatrix <- matrix(nrow = length(NewValidOTUs), ncol = 0, dimnames = list(NewValidOTUs, c()))
  
  # Create empty character weights vector:
  characterweights <- vector(mode = "numeric")
  
  # For each data set:
  for(i in 1:length(MRPList)) {
    
    # Find taxa not present in full matrix:
    taxanotinmatrix <- setdiff(NewValidOTUs, rownames(MRPList[[i]]$matrix))
    
    # Add data set block (including missing taxa as NAs) into full matrix:
    FullMRPMatrix <- cbind(FullMRPMatrix, rbind(MRPList[[i]]$matrix, matrix(NA, nrow = length(taxanotinmatrix), ncol = ncol(MRPList[[i]]$matrix), dimnames = list(taxanotinmatrix, rep(names(MRPList)[i], ncol(MRPList[[i]]$matrix)))))[rownames(FullMRPMatrix), ])
    
    # Add character weights using publication weights:
    characterweights <- c(characterweights, rep(publicationweights[i], ncol(MRPList[[i]]$matrix)))
    
  }
  
  # Add taxonomy into matrix:
  FullMRPMatrix <- cbind(FullMRPMatrix, TaxonomyMRP[rownames(FullMRPMatrix), ])
  
  # Add taxonomy weights:
  characterweights <- c(characterweights, rep(1, ncol(TaxonomyMRP)))
  
  # Add all zero outgroup:
  FullMRPMatrix <- rbind(rep("0", ncol(FullMRPMatrix)), FullMRPMatrix)
  
  # Add all zero outgroup taxon name:
  rownames(FullMRPMatrix)[1] <- "allzero"
  
  # Convert into Claddis formatted matrix:
  FullMRPMatrix <- MakeMorphMatrix(FullMRPMatrix, header = "", weights = characterweights, ordering = rep("ord", ncol(FullMRPMatrix)), equalise.weights = FALSE)
  
  # Print current processing status:
  cat("Done\nPerforming Safe Taxonomic Reduction...")
  
  # Perform STR on full matrix:
  STRdata <- SafeTaxonomicReduction(FullMRPMatrix)
  
  # Create additional STR matrix from full matrix:
  STRMRPMatrix <- FullMRPMatrix
  
  # If STR removed taxa then create reduced matrix ready for output:
  if(nrow(STRdata$str.list) > 0) STRMRPMatrix$Matrix_1$Matrix <- STRdata$reduced.matrix$Matrix_1$Matrix
  
  # Print current processing status:
  cat("Done\nCompiling and returning output...")
  
  # Compile output:
  output <- list(FullMRPMatrix, STRMRPMatrix, TaxonomyMRPTree, STRdata$str.list)
  
  # Add names:
  names(output) <- c("FullMRPMatrix", "STRMRPMatrix", "TaxonomyTree", "SafelyRemovedTaxa")
  
  # Print current processing status:
  cat("Done")
  
  # Return output:
  return(output)
  
}
