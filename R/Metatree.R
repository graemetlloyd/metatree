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
#' @param BackboneConstraint Either a file name of one of the source data sets or a Newick string to represent a backbone constraint (will enforce topology in final metatree but allows taxa not in topology to fall out inside the constraint). This is not required and the default (NULL) will mean no constraint is applied.
#' @param MonophylyConstraint Either a file name of one of the source data sets or a Newick string to represent a monophyly constraint (will enforce topology in final metatree and force taxa not in topology to fall outside the constraint). This is not required and the default (NULL) will mean no constraint is applied.
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
  InclusiveDataList <- sort(c(GetFilesForClade("matricht.html"), "Bickelmann_etal_2009a", "Caldwell_1996a", "Chen_etal_2014ba", "Chen_etal_2014bb", "deBraga_et_Rieppel_1997a", "Gauthier_etal_1988b", "Laurin_et_Reisz_1995a", "Muller_2004a", "Reisz_etal_2011a", "Rieppel_et_Reisz_1999a", "Rieppel_et_deBraga_1996a", "Young_2003a"))
  ExclusiveDataList <- c("Averianov_inpressa", "Bravo_et_Gaete_2015a", "Brocklehurst_etal_2013a", "Brocklehurst_etal_2015aa", "Brocklehurst_etal_2015ab", "Brocklehurst_etal_2015ac", "Brocklehurst_etal_2015ad", "Brocklehurst_etal_2015ae", "Brocklehurst_etal_2015af", "Bronzati_etal_2012a", "Bronzati_etal_2015ab", "Brusatte_etal_2009ba", "Campbell_etal_2016ab", "Carr_et_Williamson_2004a", "Carr_etal_2017ab", "Frederickson_et_Tumarkin-Deratzian_2014aa", "Frederickson_et_Tumarkin-Deratzian_2014ab", "Frederickson_et_Tumarkin-Deratzian_2014ac", "Frederickson_et_Tumarkin-Deratzian_2014ad", "Garcia_etal_2006a", "Gatesy_etal_2004ab", "Grellet-Tinner_2006a", "Grellet-Tinner_et_Chiappe_2004a", "Grellet-Tinner_et_Makovicky_2006a", "Knoll_2008a", "Kurochkin_1996a", "Lopez-Martinez_et_Vicens_2012a", "Lu_etal_2014aa", "Norden_etal_inpressa", "Pisani_etal_2002a", "Ruiz-Omenaca_etal_1997a", "Ruta_etal_2003ba", "Ruta_etal_2003bb", "Ruta_etal_2007a", "Selles_et_Galobart_2016a", "Sereno_1993a", "Sidor_2001a", "Skutschas_etal_inpressa", "Tanaka_etal_2011a", "Toljagic_et_Butler_2013a", "Tsuihiji_etal_2011aa", "Varricchio_et_Jackson_2004a", "Vila_etal_2017a", "Wilson_2005aa", "Wilson_2005ab", "Zelenitsky_et_Therrien_2008a")
  HigherTaxaToCollapse = c("Cymbospondylidae", "Grippiidae")
  MissingSpecies = "exclude"
  Interval = c("Triassic", "Jurassic")
  VeilLine = TRUE
  SpeciesToExclude = c("Californosaurus_perrini", "Toretocnemus_californicus", "Toretocnemus_zitteli", "Hudsonelpidia_brevirostris")
  IncludeSpecimenLevelOTUs = TRUE
  BackboneConstraint = "Moon_inpressa"
  MonophylyConstraint = NULL
  
  # New Options (requires code to actually use them)
  #
  # BackboneConstraint Newick string of backbone constraint (allows taxa not in topology). NULL as default. Allow to be an input file too.
  # MonophylyConstraint Newick string of monophyly constraint (excludes taxa not in topology). NULL as default. Allow to be an input file too.
  
  # FOR HIGHER TAXA TO COLLAPSE HAVE TO ALSO EDIT CONSTRAINT TREES TOO (AND CHECK THEY CAN EVEN MESH!)
  # CHECK FOR SPECIES THAT BELONG TO A GENUS DIFFERENT TO THE ONE IN THEIR NAME!
  # NEED TO CATCH ISSUE WHERE GENUS NUMBER IS USED FOR A SPECIES (HARD TO CHECK SO FAR DUE TO INDETERMINATES CONTINGENCY)
  # NEED SOME TEST THAT HELPS CHECK ROOT IS SENSIBLE (CAN DO WHEN MATCH TAXONOMY TO DATA SETS FOR CHUNKING)
  # NEED SOME TEST THAT HELPS DETERMINE IF MULTIPLE OCCURRENCES OF SAME TAXON AFTER RECONCILIATION IS CORRECT OR AN ERROR
  # ADD MORE COMPLEX WEIGHTS BY USING ADDITIONAL CHARACTER STATES! (EACH DATASET TOTAL WEIGHT DETERMINED BY YEAR AND DEPENDENCE THEN SUBDIVIDED ACROSS CHARACTER?) - BUT THIS SEEMS TO SLOW THINGS DRAMATICALLY MAYBE DO BY DUPLICATING CHARACTERS INSTEAD
  # MAKE STR OPTIONAL (SAVES A LITTLE TIME)
  # CHECK THERE ARE MULTIPLE TAXA PRE-RECONCILIATION
  # CHECK INDETS DO NOT GIVE MULTIPLE MATCHES
  # ADD INPUT WEIGHT OPTION (WILL WANT TO WEIGHT FOR DATA SET NOT CHARACTERS AS THESE CAN CHANGE IN PROCESSING, E.G. HAVE EVERY CHARACTER BE SAME WEIGHT IN DATA SET)
  # ADD INVERSE OPTION OF SPECIES TO EXCLUDE (SPECIES TO INCLUDE)
  # CHUNK AND WORK OUT HOW TO CALL TNT AND PARALLELISE
  # NEED GOOD WAY TO DEAL WITH WEIGHTING OF DATA SETS WITH LOTS OF CHARACTERS VERSSUS THOSE WITH FEW (NORMALISING BY N TIPS)
  
  # HOW TO DELETE DATA SETS THAT STILL CONTRIBUTE TO DEPENDENCE? (DELETE MATRIX BUT DO NOT REMOVE FROM MRP LIST)
  
  # Subfunction that gives just MRPs where matrix is still intact (has rows and columns):
  ActiveMRP <- function(MRPList) unname(which(unlist(lapply(MRPList, function(x) prod(dim(x$Matrix)))) > 0))

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
  
  # If backbone constraint is set:
  if(!is.null(BackboneConstraint)) {
    
    # Check is only a single value and stop and warn user if not:
    if(length(BackboneConstraint) > 1) stop("BackboneConstraint must be a single value. (Cannot apply two backbone constraints simultaneously.)")
    
    # Check is a string and stop and warn user if not:
    if(!is.character(BackboneConstraint)) stop("BackboneConstraint must be a text string. Reformat and try again.")
    
    # Try and work out if backbone is a Newick:
    BackboneIsNewick <- length(grep("\\(|\\)", strsplit(BackboneConstraint, ",")[[1]])) > 0
    
    # If backbone is a Newick try and read it in to check for errors:
    if(BackboneIsNewick) BackboneConstraintTree <- ape::read.tree(text = BackboneConstraint)
    
  }
  
  # If monophyly onstraint is set:
  if(!is.null(MonophylyConstraint)) {
    
    # Check is only a single value and stop and warn user if not:
    if(length(MonophylyConstraint) > 1) stop("MonophylyConstraint must be a single value. (Cannot apply two monophyly constraints simultaneously.)")
    
    # Check is a string and stop and warn user if not:
    if(!is.character(MonophylyConstraint)) stop("BackboneConstraint must be a text string. Reformat and try again.")
    
    # Try and work out if monophyly is a Newick:
    MonophylyIsNewick <- length(grep("\\(|\\)", strsplit(MonophylyConstraint, ",")[[1]])) > 0
    
    # If monophyly is a Newick trya nd read it in to check for errors:
    if(MonophylyIsNewick) MonophylyConstraintTree <- ape::read.tree(text = MonophylyConstraint)
    
  }
 
  # List of types of resolution that require finding a senior synonym:
  synonyms <- c("corrected to", "misspelling of", "objective synonym of", "obsolete variant of", "recombined as", "replaced by", "subjective synonym of")
  
  # List of types of resolution that require changing reconciliation to DELETE:
  DeletionCategories <- c("nomen dubium", "nomen vanum", "nomen nudum", "nomen oblitum", "invalid subgroup of")
  
  # Print current processing status:
  cat("Reading MRP data...")
  
  # Set working directory as MRP directory:
  setwd(MRPDirectory)
  
  # List MRP files (or just use inclusivedatalist if set):
  MRPFileList <- strsplit(ifelse(length(InclusiveDataList) > 0, paste(setdiff(sort(unique(InclusiveDataList)), sort(unique(ExclusiveDataList))), "mrp.nex", sep = "", collapse = "%%"), paste(setdiff(gsub("mrp\\.nex", "", list.files()), sort(unique(ExclusiveDataList))), "mrp.nex", sep = "", collapse = "%%")), "%%")[[1]]
  
  # If a backbone constraint is used and is a file check it is present in the source data and stop and warn user if not:
  if(!is.null(BackboneConstraint) && !BackboneIsNewick) if(is.na(match(paste(BackboneConstraint, "mrp.nex", sep = ""), MRPFileList))) stop("Backbone constraint file not found in data. Check name and try again.")
  
  # If a monophyly constraint is used and is a file check it is present in the source data and stop and warn user if not:
  if(!is.null(MonophylyConstraint) && !MonophylyIsNewick) if(is.na(match(paste(MonophylyConstraint, "mrp.nex", sep = ""), MRPFileList))) stop("Monophyly constraint file not found in data. Check name and try again.")
  
  # Read in all MRP files and store in a list (include duplicate headers to store parent sibling info later):
  MRPList <- lapply(lapply(as.list(MRPFileList), Claddis::ReadMorphNexus), function(x) {y <- list(x$Matrix_1$Matrix, x$Matrix_1$Weights, "", "", ""); names(y) <- c("Matrix", "Weights", "FileName", "Parent", "Sibling"); y})
  
  # Set names of MRP files:
  names(MRPList) <- gsub("mrp.nex", "", MRPFileList)
  
  # Print current processing status:
  cat("Done\nReading XML data...")
  
  # Set working directory as XML (i.e., metadata) directory:
  setwd(XMLDirectory)
  
  # List MRP files (or just use inslusivedatalist if set):
  XMLFileList <- strsplit(ifelse(length(InclusiveDataList) > 0, paste(setdiff(sort(unique(InclusiveDataList)), sort(unique(ExclusiveDataList))), ".xml", sep = "", collapse = "%%"), paste(setdiff(gsub("\\.xml", "", list.files()), sort(unique(ExclusiveDataList))), ".xml", sep = "", collapse = "%%")), "%%")[[1]]
  
  # Check there are no MRPs not listed as XMLs and vice versa (should return empty vector):
  MRPXMLunion <- c(setdiff(gsub("\\.xml", "", XMLFileList), gsub("mrp\\.nex", "", MRPFileList)), setdiff(gsub("mrp\\.nex", "", MRPFileList), gsub("\\.xml", "", XMLFileList)))
  
  # Stop if MRP datasets not listed as XMLs and vice versa:
  if(length(MRPXMLunion) > 0) stop(paste("Datasets do not match (MRP and XML)!:", MRPXMLunion, collapse = " "))
  
  # Read in all XML files and store in a list:
  XMLList <- lapply(as.list(XMLFileList), function(x) ReadMetatreeXML(x))
  
  # Add names to XML list:
  names(XMLList) <- gsub(".xml", "", XMLFileList)
  
  # Collapse to just pertinent information:
  XMLList <- lapply(XMLList, function(x) {y <- list(); y[["TaxonMatrix"]] <- x$SourceTree$Taxa$TagContents; y[["FileName"]] <- unname(unlist(x$SourceTree$Filename)); y[["Parent"]] <- unname(unlist(x$SourceTree$Parent)); y[["Sibling"]] <- unname(unlist(x$SourceTree$Sibling)); y})
  
  # Find any files that contain duplicated taxon names:
  FilesWithDuplicatedTaxonNames <- names(XMLList)[which(unlist(lapply(XMLList, function(x) any(duplicated(x$TaxonMatrix[, "ListValue"])))))]
  
  # If duplicate names were found stop and warn user:
  if(length(FilesWithDuplicatedTaxonNames) > 0) stop(paste("The following files contain duplicate taxon names: ", paste(FilesWithDuplicatedTaxonNames, collapse = ", "), ". Ensure all taxon names are unique and try again.", sep = ""))
  
  # Find any taxon names that do not match between MRP and XML:
  TaxonMismatches <- mapply(function(x, y) {MRPNames <- rownames(x$Matrix); XMLNames <- y$TaxonMatrix[, "ListValue"]; c(setdiff(MRPNames, XMLNames), setdiff(XMLNames, MRPNames))}, x = MRPList[names(MRPList)], y = XMLList[names(MRPList)])
  
  # Find any files with mismatching taxon names between MRP and XML:
  FilesWithTaxonMismatches <- names(TaxonMismatches)[which(unlist(lapply(TaxonMismatches, function(x) length(x))) > 0)]
  
  # If such files are found then stop and warn user:
  if(length(FilesWithTaxonMismatches) > 0) stop(paste("The following files contain mismatching taxon names between the MRP and XML versions: ", paste(FilesWithTaxonMismatches, collapse = ", "), ". Ensure all taxon names match and try again.", sep = ""))
  
  # Compile any name issues:
  NameIssues <- lapply(XMLList, function(x) {TaxonMatrix <- x$TaxonMatrix; SpacesFound <- c(grep(" ", TaxonMatrix[, "recon_name"]), grep(" ", TaxonMatrix[, "recon_no"]), grep(" ", TaxonMatrix[, "ListValue"])); EmptyValuesFound <- c(which(TaxonMatrix[, "recon_name"] == ""), which(TaxonMatrix[, "recon_no"] == ""), which(TaxonMatrix[, "ListValue"] == "")); RogueNumberCharacters <- setdiff(unique(unlist(strsplit(TaxonMatrix[, "recon_no"], ""))), c(0:9, ";", "-")); RogueNameCharacters <- setdiff(unique(c(unlist(strsplit(TaxonMatrix[, "recon_name"], "")), unlist(strsplit(TaxonMatrix[, "ListValue"], "")))), c(LETTERS, letters, 0:9, "_", ",")); y <- list(SpacesFound, EmptyValuesFound, RogueNumberCharacters, RogueNameCharacters); names(y) <- c("SpacesFound", "EmptyValuesFound", "RogueNumberCharacters", "RogueNameCharacters"); y})
  
  # Find any files with spaces in taxon names:
  FilesWithSpaces <- names(NameIssues)[unlist(lapply(NameIssues, function(x) length(x$SpacesFound))) > 0]
  
  # Find any values with empty values for taxon names:
  FilesWithEmptyValues <- names(NameIssues)[unlist(lapply(NameIssues, function(x) length(x$EmptyValuesFound))) > 0]
  
  # Files with rogue values in the recon number field:
  FilesWithRogueTaxonNumbers <- names(NameIssues)[unlist(lapply(NameIssues, function(x) length(x$RogueNumberCharacters))) > 0]
  
  # Files with rogue values in the name fields:
  FilesWithRogueTaxonNames <- names(NameIssues)[unlist(lapply(NameIssues, function(x) length(x$RogueNameCharacters))) > 0]
  
  # If issues with spaces in names stop and warn user:
  if(length(FilesWithSpaces) > 0) stop(paste("The following files contain spaces in the taxonomic reconciliation (names or numbers): ", paste(FilesWithSpaces, collapse = ", "), ". Remove spaces and try again.", sep = ""))
  
  # If issues with empty names stop and warn user:
  if(length(FilesWithEmptyValues) > 0) stop(paste("The following files contain empty values in the taxonomic reconciliation (names or numbers): ", paste(FilesWithEmptyValues, collapse = ", "), ". Ensure all values are filled and try again.", sep = ""))
  
  # If issues with rogue characters in number field stop and warn user:
  if(length(FilesWithRogueTaxonNumbers) > 0) stop(paste("The following files contain rogue values in the taxonomic reconciliation (numbers): ", paste(FilesWithRogueTaxonNumbers, collapse = ", "), ". Ensure all taxon numbers only include semicolon(s) (the separating character) or dashes (for negative values) and try again.", sep = ""))
  
  # If issues with rogue characters in name field stop and warn user:
  if(length(FilesWithRogueTaxonNames) > 0) stop(paste("The following files contain rogue values in the taxonomic reconciliation (names): ", paste(FilesWithRogueTaxonNames, collapse = ", "), ". Ensure all taxon names are formed from alphanumerics, commas (the separating character) or underscores and try again.", sep = ""))
  
  # Reconcile OTU names with XML version:
  MRPList <- mapply(function(x, y) {rownames(x$Matrix)[unlist(lapply(as.list(rownames(x$Matrix)), function(z) which(y$TaxonMatrix[, "ListValue"] == z)))] <- paste(y$TaxonMatrix[, "recon_no"], y$TaxonMatrix[, "recon_name"], sep = "%%%%"); x$FileName <- y$FileName; if(!is.null(y$Parent)) x$Parent <- y$Parent; if(!is.null(y$Sibling)) x$Sibling <- y$Sibling; x}, x = MRPList[names(MRPList)], y = XMLList[names(MRPList)], SIMPLIFY = FALSE)

  # Print current processing status:
  cat("Done\nChecking for unsampled parents and siblings...")
  
  # Extract parent and sibling names:
  ParentAndSiblingNames <- sort(unlist(lapply(as.list(unique(unname(unlist(lapply(MRPList, '[', c("Parent", "Sibling")))))), function(x) x[nchar(x) > 0])))
  
  # Warn user about any unsampled parents and/or siblings:
  if(length(setdiff(ParentAndSiblingNames, names(MRPList))) > 0) print(paste("The following parents and siblings are not in the sample (check they are correct or add them into the sample): ", paste(setdiff(ParentAndSiblingNames, names(MRPList)), collapse = ", "), sep = ""))
  
  # Print current processing status:
  cat("Done\nFinding initial multiple-taxon reconciliations...")
  
  # Subfunction to make multi-taxon reconciliations unique OTUs:
  SeparateMultiTaxonReconciliations <- function(ListBlock) {
  
    # Find comma rows (multiple taxa in initial reconciliation):
    commarows <- grep(",", rownames(ListBlock$Matrix))
    
    # If there is at least one multiple-taxon reconciliation:
    if(length(commarows) > 0) {
      
      # For each multiple-taxon reconciliation in reverse order (to avoid later rows not matching):
      for(j in rev(commarows)) {
        
        # Get multiple names of reconciliation:
        multiplenames <- strsplit(rownames(ListBlock$Matrix)[j], "%%%%")[[1]]
        
        # Get multiple-taxon numbers:
        multitaxonnumbers <- strsplit(multiplenames[1], ";")[[1]]
        
        # Get multiple-taxon names:
        multitaxonnames <- strsplit(multiplenames[2], ",")[[1]]
        
        # Check data integrity with respect to multiple-taxon values:
        if(length(multitaxonnumbers) != length(multitaxonnames)) stop(paste("Problem with multiple-taxon reconciliation(s) in ", ListBlock$FileName, " (check commas and semi-colons are correct; i.e., of same length).", sep = ""))
        
        # Add new rows at base of matrix:
        ListBlock$Matrix <- rbind(ListBlock$Matrix, matrix(rep(ListBlock$Matrix[j, ], length(multitaxonnumbers)), nrow = length(multitaxonnumbers), byrow = TRUE, dimnames = list(paste(multitaxonnumbers, multitaxonnames, sep = "%%%%"), c())))
        
        # Remove now redundant row from matrix:
        ListBlock$Matrix <- ListBlock$Matrix[-j, , drop = FALSE]
        
      }
      
    }
    
    # Return updated list block:
    return(ListBlock)
    
  }
  
  # Separate out multi-taxon reconcilations:
  MRPList <- lapply(MRPList, SeparateMultiTaxonReconciliations)
  
  # If excluding specimen-level OTUs:
  if(!IncludeSpecimenLevelOTUs) {
    
    # Print current processing status:
    cat("Done\nRemoving specimen-level OTUs...")
    
    # Convert specimen-level OTUs to taxa to DELETE:
    MRPList <- lapply(MRPList, function(x) {RowNamesToDelete <- which(unlist(lapply(strsplit(rownames(x$Matrix), split = ""), function(y) sum(y == "_") > 2))); if(length(RowNamesToDelete) > 0) rownames(x$Matrix)[RowNamesToDelete] <- "0%%%%DELETE"; x})
    
  }
  
  # Print current processing status:
  cat("Done\nRemoving taxa with initial reconciliations of \"DELETE\"...")
  
  # Remove any taxa reconciled as DELETE:
  MRPList <- lapply(MRPList, function(x) {DeleteRows <- which(unlist(lapply(strsplit(rownames(x$Matrix), split = "%%%%"), function(y) y[2])) == "DELETE"); if(length(DeleteRows) > 0) x$Matrix <- x$Matrix[-DeleteRows, , drop = FALSE]; x})
  
  # Prune matrices following deletion:
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- PisaniMRPPrune(Claddis::MakeMorphMatrix(x$Matrix, weights = x$Weights, ignore.duplicate.taxa = TRUE)); x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})

  # Print current processing status:
  cat("Done\nSearching for and collapsing pre-reconciliation duplicated taxa...")
  
  # Collapse any duplicate taxon names:
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- Claddis::MakeMorphMatrix(x$Matrix, weights = x$Weights, ignore.duplicate.taxa = TRUE); if(any(duplicated(rownames(y$Matrix_1$Matrix)))) {DuplicateNames <- setdiff(unlist(lapply(strsplit(rownames(y$Matrix_1$Matrix)[duplicated(rownames(y$Matrix_1$Matrix))], split = "%%%%"), '[', 2)), "DELETE"); if(length(DuplicateNames) > 0) cat(paste("\nDuplicate resolved OTU name(s) found in ", x$FileName, ": ", paste(DuplicateNames, collapse = ", "), ". Check this is correct.", sep = "")); y <- CollapseDuplicateTaxonMRP(y)}; x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
  
  # Print current processing status:
  cat("Done\nBuilding initial taxonomy matrix...")
  
  # Create taxonomy matrix to store all taxon resolution data:
  TaxonomyMatrix <- do.call(rbind, strsplit(unique(unname(unlist(lapply(MRPList[ActiveMRP(MRPList)], function(x) rownames(x$Matrix))))), split ="%%%%"))
  
  # Add column names:
  colnames(TaxonomyMatrix) <- c("TaxonNo", "TaxonName")
  
  # Print current processing status:
  cat("Done\nChecking for missing taxon numbers...")
  
  # If any "-1" taxa found stop and tell user:
  if(any(TaxonomyMatrix[, "TaxonNo"] == "-1")) stop(paste("The following taxa have the reconciliation number \"-1\": ", paste(TaxonomyMatrix[TaxonomyMatrix[, "TaxonNo"] == "-1", "TaxonName"], collapse = ", "), sep = ""))
  
  # Print current processing status:
  cat("Done\nBuilding initial Paleobiology Database reconciliation list...")
  
  # Create resolved taxon numbers matrix:
  ResolvedTaxonNumbers <- cbind(unique(TaxonomyMatrix[, "TaxonNo"]), PaleobiologyDBTaxaQuerier(taxon_nos = unique(TaxonomyMatrix[, "TaxonNo"]), interval = NULL))
  
  # Deal with subgenera:
  ResolvedTaxonNumbers[, "TaxonName"] <- gsub(" \\(|\\)", "", ResolvedTaxonNumbers[, "TaxonName"])
  
  # Add column names to first value (input number):
  colnames(ResolvedTaxonNumbers)[1] <- "InputNo"
  
  # If specifying an Interval:
  if(!all(is.null(Interval))) {
    
    # Do same query for just taxa in Interval:
    ResolvedTaxonNumbersInterval <- cbind(unique(TaxonomyMatrix[, "TaxonNo"]), PaleobiologyDBTaxaQuerier(taxon_nos = unique(TaxonomyMatrix[, "TaxonNo"]), interval = Interval))
    
    # Deal with subgenera:
    ResolvedTaxonNumbersInterval[, "TaxonName"] <- gsub(" \\(|\\)", "", ResolvedTaxonNumbersInterval[, "TaxonName"])
    
    # Invert variable so only includes taxa outside Interval:
    ResolvedTaxonNumbersInterval <- ResolvedTaxonNumbers[is.na(ResolvedTaxonNumbersInterval[, "TaxonName"]), ]
    
    # Find any nomen dubia etc. to delete:
    DeleteRows <- which(unlist(lapply(lapply(lapply(as.list(ResolvedTaxonNumbersInterval[, "TaxonValidity"]), match, x = DeletionCategories), sort), length)) > 0)
    
    # If there are deletes then remove them from the matrix:
    if(length(DeleteRows) > 0) ResolvedTaxonNumbersInterval <- ResolvedTaxonNumbersInterval[-DeleteRows, , drop = FALSE]
    
  }
  
  # Print current processing status:
  cat("Done\nChecking taxon names match with database version...")
  
  # Vector to store any failed matches:
  failedmatches <- vector(mode = "character")
  
  # For each initially reconciled name:
  for(i in 1:nrow(TaxonomyMatrix)) {
    
    # Get resolved name:
    resolvedname <- gsub(" ", "_", ResolvedTaxonNumbers[which(ResolvedTaxonNumbers[, "InputNo"] == TaxonomyMatrix[i, "TaxonNo"]), "TaxonName"])
    
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
  newresolutions <- setdiff(sort(unique(ResolvedTaxonNumbers[, "TaxonValidity"])), c(DeletionCategories, synonyms))
  
  # Stop if new resolutiosn found (need to add these to the resolution types above):
  if(length(newresolutions) > 0) stop(paste("New resolution type found!: ", newresolutions, sep = ""))
  
  # Print current processing status:
  cat("Done\nBuilding synonymy tables...")
  
  # Empty vector to store rows that correspond to some form of junior synonym:
  synonymrows <- c()
  
  # Find all junior synonym rows:
  for(i in synonyms) synonymrows <- sort(c(synonymrows, which(ResolvedTaxonNumbers[, "TaxonValidity"] == i)))
  
  # Set junior synonym matrix:
  JuniorSynonyms <- ResolvedTaxonNumbers[synonymrows, ]
  
  # Create empty matrix to store senior synoyms:
  SeniorSynonyms <- matrix(nrow = 0, ncol = 8, dimnames = list(c(), c("OriginalTaxonNo", "ResolvedTaxonNo", "TaxonName", "TaxonRank", "ParentTaxonNo", "TaxonValidity", "AcceptedNumber", "AcceptedName")))
  
  # Reconcile senior synonym with database:
  currenttaxa <- PaleobiologyDBTaxaQuerier(gsub("txn:", "", ResolvedTaxonNumbers[synonymrows, "AcceptedNumber"]))
  
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
  SeniorSynonyms <- currenttaxa
  
  # If using an Interval:
  if(!all(is.null(Interval))) {
    
    # Update resolved taxon numbers to valid taxa only:
    ResolvedTaxonNumbersInterval[which(!is.na(match(ResolvedTaxonNumbersInterval[, "InputNo"], JuniorSynonyms[, "InputNo"]))), c("OriginalTaxonNo", "ResolvedTaxonNo", "TaxonName", "TaxonRank", "ParentTaxonNo", "TaxonValidity", "AcceptedNumber", "AcceptedName")] <- SeniorSynonyms[match(ResolvedTaxonNumbersInterval[, "InputNo"], JuniorSynonyms[, "InputNo"])[!is.na(match(ResolvedTaxonNumbersInterval[, "InputNo"], JuniorSynonyms[, "InputNo"]))], c("OriginalTaxonNo", "ResolvedTaxonNo", "TaxonName", "TaxonRank", "ParentTaxonNo", "TaxonValidity", "AcceptedNumber", "AcceptedName")]
    
  }
  
  # Only complete this step if including specimen-level OTUs (there will not be any at this stage anyway if set as FALSE):
  if(IncludeSpecimenLevelOTUs) {
    
    # Print current processing status:
    cat("Done\nChecking validity of indeterminate taxon reconciliations...")
    
    # Get list of indeterminates:
    indeterminates <- TaxonomyMatrix[which((unlist(lapply(lapply(strsplit(TaxonomyMatrix[, "TaxonName"], "_"), '==', "aff"), sum)) + unlist(lapply(lapply(strsplit(TaxonomyMatrix[, "TaxonName"], "_"), '==', "cf"), sum)) + unlist(lapply(lapply(strsplit(TaxonomyMatrix[, "TaxonName"], "_"), '==', "indet"), sum)) + unlist(lapply(lapply(strsplit(TaxonomyMatrix[, "TaxonName"], "_"), '==', "sp"), sum))) > 0), "TaxonName"]
    
    # For each indeterminate:
    for(i in indeterminates) {
      
      # Get resolved row number:
      resolvedrownumber <- which(ResolvedTaxonNumbers[, "InputNo"] == TaxonomyMatrix[which(TaxonomyMatrix[, "TaxonName"] == i), "TaxonNo"])
      
      # If a possible invalid taxon (validity is not blank):
      if(!is.na(ResolvedTaxonNumbers[resolvedrownumber, "TaxonValidity"])) {
        
        # Get accepted number of taxon (may be NA):
        AcceptedNumber <- gsub("txn:", "", ResolvedTaxonNumbers[resolvedrownumber, "AcceptedNumber"])
        
        # If accepted number is blank (NA) stop adn warn user taxon is invalid:
        if(is.na(AcceptedNumber)) stop(paste(i, " assigned to a taxon that is invalid, consider renaming.", sep = ""))
        
        # If accepted numebr is different to input number stop and warn user taxon is synonymised:
        if(AcceptedNumber != ResolvedTaxonNumbers[resolvedrownumber, "InputNo"]) stop(paste(i, " assigned to a taxon that is invalid, consider renaming.", sep = ""))
        
      }
      
    }
    
  }

  # Print current processing status:
  cat("Done\nDeleting taxa resolved as nomen dubium and the like...")
  
  # Get input numbers that should be deleted:
  NumbersToDelete <- ResolvedTaxonNumbers[unlist(lapply(as.list(DeletionCategories), function(x) which(ResolvedTaxonNumbers[, "TaxonValidity"] == x))), "InputNo"]
  
  # As long as there are numbers to delete:
  if(length(NumbersToDelete) > 0) {
    
    # Remove any taxa to delete:
    MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {TaxonNumbers <- do.call(rbind, strsplit(rownames(x$Matrix), split = "%%%%"))[, 1]; DeleteRows <- sort(match(NumbersToDelete, TaxonNumbers)); if(length(DeleteRows) > 0) x$Matrix <- x$Matrix[-DeleteRows, , drop = FALSE]; x})
    
    # Prune matrices following deletion:
    MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- PisaniMRPPrune(Claddis::MakeMorphMatrix(x$Matrix, weights = x$Weights, ignore.duplicate.taxa = TRUE)); x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
    
  }
  
  # Print current processing status:
  cat("Done\nReplacing junior synonyms with senior synonyms...")
  
  # Rebuild junior synonyms into a vector of names:
  JuniorSynonymsVector <- paste(JuniorSynonyms[, "InputNo"], gsub(" ", "_", JuniorSynonyms[, "TaxonName"]), sep = "%%%%")
  
  # Rebuild senior synonyms into a vector of names:
  SeniorSynonymsVector <- paste(unlist(lapply(apply(SeniorSynonyms[, c("OriginalTaxonNo", "ResolvedTaxonNo")], 1, as.list), function(x) unname(gsub("txn:|var:", "", unlist(x)[!is.na(unlist(x))][1])))), gsub(" ", "_", SeniorSynonyms[, "TaxonName"]), sep = "%%%%")
  
  # Build synonym matrix:
  SynonymyMatrix <- cbind(JuniorSynonymsVector, SeniorSynonymsVector)
  
  # Replace junior synonyms with senior synonyms:
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {SynonymyRows <- sort(match(SynonymyMatrix[, 1], rownames(x$Matrix))); if(length(SynonymyRows) > 0) rownames(x$Matrix)[SynonymyRows] <- SynonymyMatrix[match(rownames(x$Matrix)[SynonymyRows], SynonymyMatrix[, 1]), 2]; x})

  # Collapse any duplicate taxa created by this substitution:
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- Claddis::MakeMorphMatrix(x$Matrix, weights = x$Weights, ignore.duplicate.taxa = TRUE); if(any(duplicated(rownames(y$Matrix_1$Matrix)))) {DuplicateNames <- rownames(y$Matrix_1$Matrix)[duplicated(rownames(y$Matrix_1$Matrix))]; y <- CollapseDuplicateTaxonMRP(y)}; x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
  
  # Prune characters made redundant by these collapses:
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- PisaniMRPPrune(Claddis::MakeMorphMatrix(x$Matrix, weights = x$Weights, ignore.duplicate.taxa = TRUE)); x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
  
  # GOT TO HERE WITH REFACTOR (BUT HAVE JUMPED AROUND A BUNCH, SO...)

  # Print current processing status:
  cat("Done\nBuilding taxonomy...")
  
  # Get a list of the valid OTU names (may be pruned down later to just those in target clade):
  ValidOTUNames <- unique(unlist(lapply(lapply(MRPList, '[[', "Matrix"), rownames)))[grep("_", unique(unlist(lapply(lapply(MRPList, '[[', "Matrix"), rownames))))]
  
  # Replace junior with senior synonyms in resolved names matrix:
  ResolvedTaxonNumbers[synonymrows, c("OriginalTaxonNo", "ResolvedTaxonNo", "TaxonName", "TaxonRank", "ParentTaxonNo", "TaxonValidity", "AcceptedNumber", "AcceptedName")] <- SeniorSynonyms
  
  # Overwrite resolved number with resolved taxon number:
  ResolvedTaxonNumbers[, "ResolvedTaxonNo"] <- gsub("txn:|var:", "", unlist(lapply(lapply(apply(ResolvedTaxonNumbers[, c("OriginalTaxonNo", "ResolvedTaxonNo")], 1, sort), rev), '[', 1)))
  
  # Remove original taxon number column:
  ResolvedTaxonNumbers <- ResolvedTaxonNumbers[, -which(colnames(ResolvedTaxonNumbers) == "OriginalTaxonNo")]
  
  # Remove deleted taxa from resolved names matrix (if there are any):
  if(length(which(!is.na(ResolvedTaxonNumbers[, "TaxonValidity"]))) > 0) ResolvedTaxonNumbers <- ResolvedTaxonNumbers[-which(!is.na(ResolvedTaxonNumbers[, "TaxonValidity"])), ]
  
  # Reformat parent taxon numbers into just numbers:
  ResolvedTaxonNumbers[, "ParentTaxonNo"] <- gsub("txn:", "", ResolvedTaxonNumbers[, "ParentTaxonNo"])
  
  # Collapse resolved matrix to just field with values (i.e., drop valid and senior synonym columns):
  ResolvedTaxonNumbers <- ResolvedTaxonNumbers[, c("ResolvedTaxonNo", "TaxonName", "TaxonRank", "ParentTaxonNo")]
  
  # If doing something with missing species (i.e., those not currently included as OTUs, but existing in target clade):
  if(MissingSpecies != "exclude") {
    
    # Find all children of target clade:
    AllChildren <- PaleobiologyDBChildFinder(taxon_nos = "1", taxon_names = TargetClade, validonly = TRUE, returnrank = "3", interval = Interval)
    
    # Deal with subgenera:
    AllChildren[, "TaxonName"] <- gsub(" \\(|\\)", "", AllChildren[, "TaxonName"])
    
    # If inserting all missing species get all possible species parent numbers:
    if(MissingSpecies == "all") CurrentSpeciesParentNumbers <- unique(c(gsub("txn:", "", AllChildren[, "ParentTaxonNo"]), ResolvedTaxonNumbers[which(ResolvedTaxonNumbers[, "TaxonRank"] == 3), "ParentTaxonNo"]))
    
    # If only inserting missing species at genus-level find parent numbers of all current species (i.e., potential genera to add):
    if(MissingSpecies == "genus") CurrentSpeciesParentNumbers <- unique(ResolvedTaxonNumbers[which(ResolvedTaxonNumbers[, "TaxonRank"] == 3), "ParentTaxonNo"])
    
    # Find any parents not already present in resolved numbers matrix:
    AsYetUnsampledSpeciesParents <- setdiff(CurrentSpeciesParentNumbers, ResolvedTaxonNumbers[, "ResolvedTaxonNo"])
    
    # If such parents exist:
    if(length(AsYetUnsampledSpeciesParents) > 0) {
      
      # Find unsampled species parents:
      CurrentSpeciesParents <- PaleobiologyDBTaxaQuerier(taxon_nos = AsYetUnsampledSpeciesParents)
      
      # Deal with subgenera:
      CurrentSpeciesParents[, "TaxonName"] <- gsub(" \\(|\\)", "", CurrentSpeciesParents[, "TaxonName"])
      
      # Find rows corresponding to valid genera:
      ValidGenusRows <- intersect(which(is.na(CurrentSpeciesParents[, "TaxonValidity"])), which(CurrentSpeciesParents[, "TaxonRank"] == "5"))
      
      # If there are valid genera then add these to resolved taxon numbers:
      if(length(ValidGenusRows) > 0) ResolvedTaxonNumbers <- rbind(ResolvedTaxonNumbers, cbind(unname(gsub("txn:|var:", "", unlist(lapply(lapply(lapply(apply(CurrentSpeciesParents[ValidGenusRows, c("OriginalTaxonNo", "ResolvedTaxonNo"), drop = FALSE], 1, list), unlist), sort, decreasing = TRUE), '[', 1)))), CurrentSpeciesParents[ValidGenusRows, c("TaxonName", "TaxonRank")] , gsub("txn:", "", CurrentSpeciesParents[ValidGenusRows, "ParentTaxonNo"])))
      
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
      CurrentGenusNumbers <- ResolvedTaxonNumbers[which(ResolvedTaxonNumbers[, "TaxonRank"] == 5), "ResolvedTaxonNo"]
      
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
    ChildrenToAdd <- setdiff(NewChildren[, 1], ResolvedTaxonNumbers[, "ResolvedTaxonNo"])
    
    # If there are children to add then add them to resolved taxon numbers:
    if(length(ChildrenToAdd) > 0) ResolvedTaxonNumbers <- rbind(ResolvedTaxonNumbers, NewChildren[match(ChildrenToAdd, NewChildren[, 1]), ])
    
  }
  
  # Get initial parent child relationships based on OTUs:
  parentchildrelationships <- paste(unlist(lapply(strsplit(ValidOTUNames, "%%%%"), '[', 1)), ResolvedTaxonNumbers[match(unlist(lapply(strsplit(ValidOTUNames, "%%%%"), '[', 1)), ResolvedTaxonNumbers[, "ResolvedTaxonNo"]), "ParentTaxonNo"], sep = " belongs to ")
  
  # If including specimen level OTUs:
  if(IncludeSpecimenLevelOTUs) {
    
    # Find which rows correspond to indeterminate and sp taxa (i.e., those where parent should be initial reconciliation):
    indetsandsps <- sort(c(which(unlist(lapply(lapply(lapply(lapply(lapply(strsplit(ValidOTUNames, "%%%%"), '[', 2), strsplit, split = "_"), unlist), '==', "indet"), any))), which(unlist(lapply(lapply(lapply(lapply(lapply(strsplit(ValidOTUNames, "%%%%"), '[', 2), strsplit, split = "_"), unlist), '==', "sp"), any)))))
    
    # If such taxa exist then update parent child relationships accordingly:
    if(length(indetsandsps) > 0) parentchildrelationships[indetsandsps] <- paste(unlist(lapply(strsplit(ValidOTUNames[indetsandsps], "%%%%"), '[', 1)), unlist(lapply(strsplit(ValidOTUNames[indetsandsps], "%%%%"), '[', 1)), sep = " belongs to ")
    
  }
  
  # Get list of new children (for which parents are needed) - excludes "Life" which has no parent:
  newchildren <- setdiff(unlist(lapply(strsplit(parentchildrelationships, " belongs to "), '[', 2)), "28595")
  
  # As long as there are still children in need of parents:
  while(length(newchildren) > 0) {
    
    # Find any numbers missing for the taxonomy name resolution matrix:
    missingfromresolutions <- newchildren[which(is.na(match(newchildren, ResolvedTaxonNumbers[, "ResolvedTaxonNo"])))]
    
    # If there are such numbers:
    if(length(missingfromresolutions) > 0) {
      
      # Get raw query data for new names
      rawquery <- PaleobiologyDBTaxaQuerier(taxon_nos = missingfromresolutions)
      
      # Deal with subgenera:
      rawquery[, "TaxonName"] <- gsub(" \\(|\\)", "", rawquery[, "TaxonName"])
      
      # Add formatted results of query to resolved names matrix:
      ResolvedTaxonNumbers <- rbind(ResolvedTaxonNumbers, cbind(gsub("txn:|var:", "", unname(unlist(lapply(lapply(lapply(apply(rawquery[, c("OriginalTaxonNo", "ResolvedTaxonNo"), drop = FALSE], 1, list), unlist), sort, decreasing = TRUE), '[', 1)))), rawquery[, c("TaxonName", "TaxonRank"), drop = FALSE], gsub("txn:", "", rawquery[, "ParentTaxonNo"])))
      
    }
    
    # Add new parent child relationships to list:
    parentchildrelationships <- c(parentchildrelationships, paste(ResolvedTaxonNumbers[match(newchildren, ResolvedTaxonNumbers[, "ResolvedTaxonNo"]), "ResolvedTaxonNo"], ResolvedTaxonNumbers[match(newchildren, ResolvedTaxonNumbers[, "ResolvedTaxonNo"]), "ParentTaxonNo"], sep = " belongs to "))
    
    # Update new children:
    newchildren <- setdiff(ResolvedTaxonNumbers[match(newchildren, ResolvedTaxonNumbers[, "ResolvedTaxonNo"]), "ParentTaxonNo"], "28595")
    
  }
  
  # If Life is missing then add it at bottom:
  if(all(!ResolvedTaxonNumbers[, "ResolvedTaxonNo"] == "28595")) ResolvedTaxonNumbers <- rbind(ResolvedTaxonNumbers, c("28595", "Life", "25", NA))
  
  # Convert parent-child relationships into a matrix (columns for child and parent):
  parentchildmatrix <- matrix(unlist(strsplit(parentchildrelationships, split = " belongs to ")), ncol = 2, byrow = TRUE, dimnames = list(c(), c("Child", "Parent")))
  
  # Update parent-child matrix with child names:
  parentchildmatrix[, "Child"] <- ResolvedTaxonNumbers[match(parentchildmatrix[, "Child"], ResolvedTaxonNumbers[, "ResolvedTaxonNo"]), "TaxonName"]
  
  # Update parent-child matrix with parent names:
  parentchildmatrix[, "Parent"] <- ResolvedTaxonNumbers[match(parentchildmatrix[, "Parent"], ResolvedTaxonNumbers[, "ResolvedTaxonNo"]), "TaxonName"]
  
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
  MRPList <- lapply(MRPList, function(x) {NamesToCheck <- rownames(x$Matrix); NamesToCheck <- gsub("_\\(|\\)", "", NamesToCheck); rownames(x$Matrix) <- NamesToCheck; x})
  
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
  if(!all(is.null(Interval))) NewValidOTUs <- setdiff(NewValidOTUs, gsub(" ", "_", ResolvedTaxonNumbersInterval[, "TaxonName"]))
  
  # Can now strip out numbers from taxon names:
  for(i in 1:length(MRPList)) if(!is.null(rownames(MRPList[[i]]$Matrix))) rownames(MRPList[[i]]$Matrix) <- unlist(lapply(strsplit(rownames(MRPList[[i]]$Matrix), "%%%%"), '[', 2))
  
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
  datasetswithsupraspecificOTUs <- which(unlist(lapply(lapply(lapply(lapply(MRPList, '[[', "Matrix"), rownames), intersect, y = colnames(TaxonomyMRP)), length)) > 0)
  
  # If such data sets exist:
  if(length(datasetswithsupraspecificOTUs) > 0) {
    
    # For each such data set:
    for(i in datasetswithsupraspecificOTUs) {
      
      # Find higher taxa that will need to be replaced:
      highertaxatoreplace <- intersect(rownames(MRPList[[i]]$Matrix), colnames(TaxonomyMRP))
      
      # For each higher taxon to replace:
      for(j in highertaxatoreplace) {
        
        # Find substitue names from taxonomy:
        substitutenames <- names(which(TaxonomyMRP[, j] == 1))
        
        # Add these to end of matrix using coding for higher taxon:
        MRPList[[i]]$Matrix <- rbind(MRPList[[i]]$Matrix, matrix(rep(MRPList[[i]]$Matrix[j, ], length(substitutenames)), nrow = length(substitutenames), byrow = TRUE, dimnames = list(substitutenames, c())))
        
        # Remove now replaced higher taxon from matrix:
        MRPList[[i]]$Matrix <- MRPList[[i]]$Matrix[-which(rownames(MRPList[[i]]$Matrix) == j), , drop = FALSE]
        
      }
      
    }
    
  }
  
  # Print current processing status:
  cat("Done\nRetracting subspecies into species...")
  
  # Replace all subspecies with their species name:
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {UnderscoreAndCapitalCounts <- matrix(unlist(lapply(strsplit(rownames(x$Matrix), split = ""), function(y) c(sum(y == "_"), length(grep("[:A-Z:]", y))))), ncol = 2, byrow = TRUE, dimnames = list(c(), c("Underscores", "Capitals"))); SubspeciesRows <- intersect(which(UnderscoreAndCapitalCounts[, "Underscores"] == 2), which(UnderscoreAndCapitalCounts[, "Capitals"] == 1)); if(length(SubspeciesRows) > 0) rownames(x$Matrix)[SubspeciesRows] <- unlist(lapply(strsplit(rownames(x$Matrix)[SubspeciesRows], split = "_"), function(z) paste(z[1:2], collapse = "_"))); x})
  
  # Collapse any duplicate taxon names:
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- Claddis::MakeMorphMatrix(x$Matrix, weights = x$Weights, ignore.duplicate.taxa = TRUE); if(any(duplicated(rownames(y$Matrix_1$Matrix)))) {DuplicateNames <- rownames(y$Matrix_1$Matrix)[duplicated(rownames(y$Matrix_1$Matrix))]; if(length(DuplicateNames) > 0) cat(paste("\nDuplicate resolved OTU name(s) found post higher-taxon substitution in ", x$FileName, ": ", paste(DuplicateNames, collapse = ", "), ". Check this is correct.", sep = "")); y <- CollapseDuplicateTaxonMRP(y)}; x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
  
  # Print current processing status:
  cat("Done\nFurther tidying of taxonomy...")
  
  # If applying an Interval:
  if(!all(is.null(Interval))) {
    
    # Find any rows to delete (because they represent taxa outside the Interval):
    RowsToDelete <- sort(match(gsub(" ", "_", ResolvedTaxonNumbersInterval[ResolvedTaxonNumbersInterval[, "TaxonRank"] == "3", "TaxonName"]), rownames(TaxonomyMRP)))
    
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
  TaxaToDelete <- setdiff(unlist(lapply(lapply(MRPList, '[[', "Matrix"), rownames)), NewValidOTUs)
  
  # If applying an Interval then add taxa outside of it to the deletes list:
  if(!all(is.null(Interval))) TaxaToDelete <- unique(c(TaxaToDelete, gsub(" ", "_", ResolvedTaxonNumbersInterval[ResolvedTaxonNumbersInterval[, "TaxonRank"] == "3", "TaxonName"])))
  
  # If there are species to exclude:
  if(length(SpeciesToExclude) > 0) {
    
    # Build vector of all current OTU names:
    OTUNames <- unique(unlist(lapply(MRPList, function(x) rownames(x$Matrix))))
    
    # Find any missing names (in exclude list but not in tree):
    MissingNames <- setdiff(SpeciesToExclude, OTUNames)
    
    # If any are found stop and warn user:
    if(length(MissingNames) > 0) stop(paste("The following SpeciesToExclude were not actually found in the data: ", paste(MissingNames, collapse = ", "), ". Check they are spelled correctly and try again.", sep = ""))
    
    # Add species to exclude to taxa to delete:
    TaxaToDelete <- unique(c(TaxaToDelete, SpeciesToExclude))
    
    # Remove species to exclude from Taxonomy MRP (as lomg as they are still there):
    if(length(intersect(SpeciesToExclude, rownames(TaxonomyMRP))) > 0) TaxonomyMRP <- TaxonomyMRP[-match(intersect(SpeciesToExclude, rownames(TaxonomyMRP)), rownames(TaxonomyMRP)), , drop = FALSE]
    
    # Find any columns to delete (duplicated, autapomorphic or constant):
    ColumnsToDalete <- unique(c(which(duplicated(apply(TaxonomyMRP, 2, paste, collapse = ""))), unname(which(apply(TaxonomyMRP, 2, sum) < 2))))
    
    # If columns are to be deleted then delete them:
    if(length(ColumnsToDalete) > 0) TaxonomyMRP <- TaxonomyMRP[, -ColumnsToDalete]
  
  }
  
  # Delete taxa from every matrix they occur in:
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {DeleteRows <- match(intersect(TaxaToDelete, rownames(x$Matrix)), rownames(x$Matrix)); if(length(DeleteRows) > 0) x$Matrix <- x$Matrix[-DeleteRows, , drop = FALSE]; x})
  
  # Prune redundant characters from matrices following taxon deletion(s):
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- PisaniMRPPrune(Claddis::MakeMorphMatrix(x$Matrix, weights = x$Weights, ignore.duplicate.taxa = TRUE)); x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
  
  # Print current processing status:
  cat("Done\nProducing taxonomy tree...")
  
  # Duplicated Taxonomy MRP to create a collapsable version for generating taxonomy Newick string:
  TaxonomyMRPNewick <- TaxonomyMRP
  
  # Get order of columns to collapse to form MRP
  columncollapseorder <- order(apply(TaxonomyMRPNewick, 2, sum))
  
  # For each column ("clade") in order from smallest to largest:
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

  # If there are higher taxa to collapse:
  if(length(HigherTaxaToCollapse) > 0) {
    
    # Print current processing status:
    cat("Done\nCollapsing higher taxa...")
    
    # Find any higher taxa actually present in target clade:
    HigherTaxaInTargetClade <- colnames(TaxonomyMRP)
    
    # Find any missing names (in collapse list but not in taxonomy):
    MissingHigherTaxa <- setdiff(HigherTaxaToCollapse, HigherTaxaInTargetClade)
    
    # If any are found stop and warn user:
    if(length(MissingHigherTaxa) > 0) stop(paste("The following HigherTaxaToCollapse were not actually found in the data: ", paste(MissingHigherTaxa, collapse = ", "), ". Check they are spelled correctly and are valid (according to the Paleobiology Database) and try again.", sep = ""))
    
    # Check the clades are all unique (not internested) and if not then stop and warn user:
    if(any(duplicated(unlist(lapply(as.list(HigherTaxaToCollapse), function(x) rownames(TaxonomyMRP[TaxonomyMRP[, colnames(TaxonomyMRP) == x] == 1, ])))))) stop("HigherTaxaToCollapse contains clades that are internested. Remove the internested clades and try again.")
    
    # Build taxa to rename matrix:
    TaxaToRenameMatrix <- do.call(rbind, lapply(as.list(HigherTaxaToCollapse), function(x) unname(cbind(x, rownames(TaxonomyMRP[TaxonomyMRP[, colnames(TaxonomyMRP) == x] == 1, ])))))
    
    # Build list of each clade's species compliment:
    CladeContentsList <- lapply(as.list(TaxonomyMRPTree$node.label), function(x) {NodeNumber <- which(TaxonomyMRPTree$node.label == x) + Ntip(TaxonomyMRPTree); TaxonomyMRPTree$tip.label[FindDescendants(n = NodeNumber, tree = TaxonomyMRPTree)]})
    
    # Add names of clades:
    names(CladeContentsList) <- TaxonomyMRPTree$node.label
    
    # Find any subsumed clades (to be removed from taxonomy MRP):
    SubsumedClades <- unlist(lapply(as.list(HigherTaxaToCollapse), function(x) {CurrentClade <- which(names(CladeContentsList) == x); TempCladeContents <- CladeContentsList[-CurrentClade]; names(which(unlist(lapply(TempCladeContents, function(x) length(setdiff(x, CladeContentsList[[CurrentClade]])))) == 0))}))
    
    # Build block to add to taxonomy MRP out of first taxon inside each clade to collapse:
    BlockToAddToTaxonomyMRP <- do.call(rbind, lapply(as.list(HigherTaxaToCollapse), function(x) TaxonomyMRP[TaxaToRenameMatrix[which(TaxaToRenameMatrix[, 1] == x)[1], 2], ]))
    
    # Add uppercase rownames to
    rownames(BlockToAddToTaxonomyMRP) <- toupper(HigherTaxaToCollapse)
    
    # Add to taxonomy MRP:
    TaxonomyMRP <- rbind(TaxonomyMRP, BlockToAddToTaxonomyMRP)
    
    # Collapse taxonomy MRP down by removing clades and the species from the collapsed clades:
    TaxonomyMRP <- TaxonomyMRP[-unlist(lapply(as.list(TaxaToRenameMatrix[, 2]), function(x) which(rownames(TaxonomyMRP) == x))), -unlist(lapply(as.list(c(HigherTaxaToCollapse, SubsumedClades)), function(x) which(colnames(TaxonomyMRP) == x))), drop = FALSE]
    
    # Remove all but one collapsed taxa from each clade from the tree:
    TaxonomyMRPTree <- drop.tip(TaxonomyMRPTree, TaxaToRenameMatrix[-unlist(lapply(as.list(HigherTaxaToCollapse), function(x) which(TaxaToRenameMatrix[, 1] == x)[1])), 2])
    
    # Ladderize taxonomy tree for neatness!:
    TaxonomyMRPTree <- ladderize(TaxonomyMRPTree)
    
    # Build tips to replace matrix:
    TipsToReplaceMatrix <- do.call(rbind, lapply(as.list(HigherTaxaToCollapse), function(x) TaxaToRenameMatrix[which(TaxaToRenameMatrix[, 1] == x)[1], ]))
    
    # Replace tip names in tree:
    TaxonomyMRPTree$tip.label[unlist(lapply(apply(TipsToReplaceMatrix, 1, as.list), function(x) {x <- unlist(x); which(TaxonomyMRPTree$tip.label == x[2])}))] <- toupper(TipsToReplaceMatrix[, 1])
    
    # Replace names in MRP matrices:
    MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {CurrentRownames <- rownames(x$Matrix); NamesToReplace <- intersect(CurrentRownames, TaxaToRenameMatrix[, 2]); if(length(NamesToReplace) > 0) rownames(x$Matrix)[match(NamesToReplace, rownames(x$Matrix))] <- toupper(TaxaToRenameMatrix[match(NamesToReplace, TaxaToRenameMatrix[, 2]), 1]); x})
    
    # Collapse any duplicate taxa created by this substitution (very likely!):
    MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- Claddis::MakeMorphMatrix(x$Matrix, weights = x$Weights, ignore.duplicate.taxa = TRUE); if(any(duplicated(rownames(y$Matrix_1$Matrix)))) {DuplicateNames <- rownames(y$Matrix_1$Matrix)[duplicated(rownames(y$Matrix_1$Matrix))]; y <- CollapseDuplicateTaxonMRP(y)}; x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
    
    # Update new valid OTUs:
    NewValidOTUs <- sort(rownames(TaxonomyMRP))
    
    # HAVE TO EDIT CONSTRAINT TREES TOO, INCLUDING DELETING TAXA!
    
  }
  
  # If specimen-level OTUs are included:
  if(IncludeSpecimenLevelOTUs) {
    
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
    
  }
  
  # For empty data sets make sure matrix is zero-by-zero and weights have no lengths:
  MRPList <- lapply(MRPList, function(x) {MatrixSize <- nrow(x$Matrix) * ncol(x$Matrix); if(MatrixSize == 0) {x$Matrix <- matrix(nrow = 0, ncol = 0); x$Weights <- vector(mode = "numeric")}; x})
  
  # Print current processing status:
  cat("Done\nGetting weighting data (publication year and dependencies)...")
  
  # Add publication year to each data set (in presses are ascribed the current year):
  MRPList <- lapply(MRPList, function(x) {x$PublicationYear <- gsub("[:A-Z:a-z:]|_", "", gsub("inpress", strsplit(as.character(Sys.Date()), "-")[[1]][1], x$FileName)); x})
  
  # Find any missing parents:
  MissingParents <- setdiff(unique(unname(unlist(lapply(MRPList, function(x) x$Parent[nchar(x$Parent) > 0])))), names(MRPList))
  
  # Find all parent data set names:
  ParentDataSets <- sort(unique(unname(unlist(lapply(MRPList, function(x) x$Parent[nchar(x$Parent) > 0])))))
  
  # Get child data sets for each parent:
  ChildDataSets <- lapply(as.list(ParentDataSets), function(x) unname(unlist(mapply(function(x, y) y$FileName[y$Parent == x], x = x, y = MRPList))))
  
  # Add names to child data sets:
  names(ChildDataSets) <- ParentDataSets
  
  # Now include any grandchildren, greatgrandchildren etc.:
  ChildDataSets <- lapply(ChildDataSets, function(x) sort(unique(c(x, unname(unlist(ChildDataSets[intersect(x, names(ChildDataSets))]))))))
  
  # Add sibling relationships to data sets with shared parents and update parents field with grandparents, greatgrandparents etc.:
  MRPList <- lapply(MRPList, function(x) {SiblingVector <- c(x$Sibling, setdiff(ChildDataSets[[match(x$Parent, ParentDataSets)]], x$FileName)); if(any(nchar(SiblingVector)) > 0) SiblingVector <- SiblingVector[nchar(SiblingVector) > 0]; x$Sibling <- unique(SiblingVector); x$Parent <- names(which(unlist(lapply(ChildDataSets, function(y) length(intersect(y, x$FileName)))) > 0)); x})
  
  # Get any redundnat parent data sets (all taxa included in at least one child data set):
  RedundantParents <- unique(unname(unlist(lapply(MRPList[ActiveMRP(MRPList)], function(x) {ActiveParent <- setdiff(x$Parent[nchar(x$Parent) > 0], MissingParents); if(length(ActiveParent) > 0) if(length(setdiff(rownames(MRPList[[ActiveParent]]$Matrix), rownames(x$Matrix))) == 0) x$Parent}))))
  
  # If redundant parents were found then collapsee these data sets back to empty matrix and weights::
  if(length(RedundantParents) > 0) MRPList[RedundantParents] <- lapply(MRPList[RedundantParents], function(x) {x$Matrix <- matrix(ncol = 0, nrow = 0); x$Weights <- vector(mode = "numeric"); x})
  
  # Find any remaining active parent data sets:
  ActiveParents <- names(which(unlist(lapply(MRPList[setdiff(ParentDataSets, MissingParents)], function(x) nrow(x$Matrix) * ncol(x$Matrix))) > 0))
  
  # If active parents remain:
  if(length(ActiveParents) > 0) {
    
    # Add children of active parent to siblings:
    MRPList[ActiveParents] <- lapply(MRPList[ActiveParents], function(x) {SiblingVector <- unique(c(ChildDataSets[[x$FileName]], x$Sibling)); x$Sibling <- unique(SiblingVector[nchar(SiblingVector) > 0]); x})
    
    # Add parent as sibling of offspring:
    MRPList <- lapply(MRPList, function(x) {SiblingVector <- c(x$Sibling, intersect(x$Parent, ActiveParents)); if(any(nchar(SiblingVector) > 0)) SiblingVector <- SiblingVector[nchar(SiblingVector) > 0]; x$Sibling <- sort(unique(SiblingVector)); x})
    
  }
  
  # Find any empty data sets to remove:
  RemovedSourceData <- sort(names(which(unlist(lapply(MRPList, function(x) nrow(x$Matrix) * ncol(x$Matrix))) == 0)))
  
  # Remove data sets from MRPList:
  MRPList[RemovedSourceData] <- NULL
  
  # Remove any dead siblings:
  MRPList <- lapply(MRPList, function(x) {x$Sibling <- intersect(x$Sibling, names(MRPList)); x})
  
  # If usinga  veil line:
  if(VeilLine) {
    
    # Print current processing status:
    cat("Done\nApplying veil line...")
    
    # Start with currnet year as veil year:
    CurrentVeilYear <- as.numeric(strsplit(as.character(Sys.Date()), "-")[[1]][1])
    
    # Set current taxa included as being from current veil year to present:
    CurrentTaxaIncluded <- unique(unlist(lapply(MRPList[as.numeric(unlist(lapply(MRPList, function(x) x$PublicationYear))) >= CurrentVeilYear], function(x) rownames(x$Matrix))))
    
    # Make a stop point (where all taxa are sampled):
    StopPoint <- length(unique(unlist(lapply(MRPList, function(x) rownames(x$Matrix)))))
    
    # While not all taxa are included in current sample:
    while(length(CurrentTaxaIncluded) < StopPoint) {
      
      # Increment one year back in time:
      CurrentVeilYear <- CurrentVeilYear - 1

      # Update current taxa included:
      CurrentTaxaIncluded <- unique(unlist(lapply(MRPList[as.numeric(unlist(lapply(MRPList, function(x) x$PublicationYear))) >= CurrentVeilYear], function(x) rownames(x$Matrix))))
      
    }
    
    # Find any data sets to remove (from veil year or older):
    DataSetsToRemove <- which(as.numeric(unlist(lapply(MRPList, function(x) x$PublicationYear))) <= CurrentVeilYear)
    
    # If data sets to remove:
    if(length(DataSetsToRemove) > 0) {
      
      # Add to removed surce data vector:
      RemovedSourceData <- sort(c(RemovedSourceData, names(MRPList)[DataSetsToRemove]))
      
      # Remove from MRP list:
      MRPList <- MRPList[-DataSetsToRemove]
      
      # Remove any new dead siblings:
      MRPList <- lapply(MRPList, function(x) {x$Sibling <- intersect(x$Sibling, names(MRPList)); x})
    
    }

  # If not using veil line:
  } else {
    
    # Set current veil year as oldest data set:
    CurrentVeilYear <- min(as.numeric(unlist(lapply(MRPList, function(x) x$PublicationYear))))
    
  }
  
  
  
  
  
  
  
  
  ######
  
  # Print current processing status:
  cat("Done\nCalculating weights...")
  
  # Equation 1 in Supplementary Information of Lloyd et al. (2016):
  publicationyearweights <- 10 * 2^(0.5 * (publicationyears - CurrentVeilYear))
  
  # Get independece weights (using siblings data):
  independenceweights <- 1 / apply(rbind(rep(1, length(MRPList)), unlist(lapply(lapply(lapply(siblingsdata, strsplit, split = "%%%%"), unlist), length))), 2, max)
  
  # Get N characters weight (to avoid uncertain data sets swamping the signal):
  Ncharacterweights <- 1 / apply(rbind(rep(1, length(MRPList)), unlist(lapply(lapply(MRPList, '[[', "Matrix"), ncol))), 2, max)
  
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
    taxanotinmatrix <- setdiff(NewValidOTUs, rownames(MRPList[[i]]$Matrix))
    
    # Add data set block (including missing taxa as NAs) into full matrix:
    FullMRPMatrix <- cbind(FullMRPMatrix, rbind(MRPList[[i]]$Matrix, matrix(NA, nrow = length(taxanotinmatrix), ncol = ncol(MRPList[[i]]$Matrix), dimnames = list(taxanotinmatrix, rep(names(MRPList)[i], ncol(MRPList[[i]]$Matrix)))))[rownames(FullMRPMatrix), ])
    
    # Add character weights using publication weights:
    characterweights <- c(characterweights, rep(publicationweights[i], ncol(MRPList[[i]]$Matrix)))
    
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
  output <- list(FullMRPMatrix, STRMRPMatrix, TaxonomyMRPTree, STRdata$str.list, RemovedSourceData, CurrentVeilYear)
  
  # Add names:
  names(output) <- c("FullMRPMatrix", "STRMRPMatrix", "TaxonomyTree", "SafelyRemovedTaxa", "RemovedSourceData", "CurrentVeilYear")
  
  # Print current processing status:
  cat("Done")
  
  # Return output:
  return(output)
  
}
