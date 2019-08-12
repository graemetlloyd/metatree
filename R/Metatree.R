#' Builds metatree from source data
#'
#' @descrioption
#'
#' Builds a metatree data set from a set of source data.
#'
#' @param MRPDirectory The directory in which the MRP files are to be read from.
#' @param XMLDirectory The directory in which the XML files are to be read from.
#' @param TargetClade The name of the target clade of the metatree (e.g., "Dinosauria").
#' @param InclusiveDataList A vector of the data sest to include in the metatree. Can be left empty to just read all files in \code{MRPDirectory} abd \code{XMLDirectory}.
#' @param ExclusiveDataList A vector of any data sets to exclude from the metatree. Can be left empty if all data sets in \code{MRPDirectory} abd \code{XMLDirectory} are valid. (Intended to exclude things like oogenera or footprint analyses, other supertree data sets etc.)
#' @param HigherTaxaToCollapse Vector of any higher taxa to collapse (e.g., if you are focused on relationships in a stem-group). NB: It is very important that these are safely monophyletic or the results can be compromised if not.
#' @param SpeciesToExclude Vector of any species to be excluded from the final metatree. E.g., Eshanosaurus, Ricardoestesia.
#' @param MissingSpecies What to do with species assigned to the target clade, but not present in the source data. Options are: "exclude" (excludes these missing species), "genus" (include those species in a genus-level polytomy if the genus is sampled in the source data), and "all"
#' @param Interval If restricting the sample to a specific interval of geologic time then use this option (passed to \link{PaleobiologyDBChildFinder} which should be consulted for formatting). Default is NULL (no restriction on ages of tips to be included).
#' @param VeilLine A logical indiicating whther to remove older data sets that do not increase taxonomic coverage (TRUE; the default) or not (FALSE). See Lloyd et al. (2016).
#' @param IncludeSpecimenLevelOTUs A logical indicating whether specimen-level OTUs should (TRUE; the default) or should not (FALSE) be included in the metatree.
#' @param BackboneConstraint The file name of one of the source data sets to represent a backbone constraint (will enforce topology in final metatree but allows taxa not in topology to fall out inside the constraint). This is not required and the default (NULL) will mean no constraint is applied.
#' @param MonophylyConstraint The file name of one of the source data sets to represent a monophyly constraint (will enforce topology in final metatree and forces taxa not in topology to fall outside the constraint). This is not required and the default (NULL) will mean no constraint is applied.
#' @param RelativeWeights A numeric vector of four values (default \code{c(1, 1, 1, 1)}) giving the respective weights to use for: 1) the input weights (the weights read in from the source MRP files), 2) the publication year weights (from equation 1 in the supplement of Lloyd et al. 2016), 3) the data set dependency weights (1 / the number of "sibling" data sets; see Lloyd et al. 2016), and 4) the within-matrix weights of individual clades (1 / number of conflicitng clades). Zeroes exclude particular weighting types. E.g., to only use input weights use \code{c(1, 0, 0, 0)}.
#' @param WeightCombination How to combine the weights above. Must be one of either "product" or "sum". Note product will exclude zero weight values to avoid zero weight output. E.g., if only using input weights the result of combining weights will not be all zeroes simply because the other types of weight are set at zero.
#' @param ReportContradictionsToScreen Logical indicating whether or not to print any taxonomy-phylogeny contradictions to the screen.
#'
#' @details
#'
#' \bold{Introduction}
#'
#' Broadly speaking this function is an implementation and extension of the approach to generating composite phylogeneic trees laid out in Lloyd et al. (2016), which itself builds on Lloyd et al. (2008), namely the "metatree" approach. Metatrees primarily differ from formal supertrees in that, instead of published trees (figures in source publications), the input data are the original character-taxon matrices and/or sequence alignments. This is superior to supertrees if you want to: 1) standardise the way input data are analysed, 2) choose the optimality criterion instead of being forced to use whatever the original study used, 3) include non-focal species that may have been removed from published figures (improving taxon overlap), and 4) more properly incorporate phylogenetic uncertainty rather than being restricted to the use of consensus topologies.
#'
#' \bold{Formatting input data}
#'
#' The implementation here assumes this data set reanalysis has already been performed, and the results have been encoded in NEXUS (Maddison et al. 1997) format using Matrix Representation with Parsimony (MRP; Baum 1992; Ragan 1992).
#.
#. XML foramtting
#'
#' \bold{Input directories}
#'
#' What these are and default of using all files, matching names etc.
#'
#' \bold{Data set inclusion and exclusion}
#'
#' Inclusive and exclusive lists.
#'
#' \bold{Taxonomic options}
#'
#' Collapsing higher taxa and species to remove. Levels of inclusivity (all, exclude, genus).
#'
#' \bold{Specifying a sampling interval}
#'
#' Time and why desirable (stem-groups). See also collapsing higher taxa. Issue with specimen-level OTUs.
#'
#' \bold{Use of veil lines}
#'
#' Alongside a priori exclusion of data sets can also use a veil line. Advantages: reduces final data set size, avoids older signal completely (although can also be dealt with with weights (see below).
#'
#' \bold{Specimen-level OTUs}
#'
#' As this is aimed at fossil data a common issue is use pf specimen-level OTUs in source data sets that (currently) lack a formal species designation. Issues with time and harder to check validity (may have since been named).
#'
#' \bold{Applying constraint(s)}
#'
#' Two ways to do this. MRP is better as original data can be passed through pipeline plus MRP allows more complex forms of constraint than regular topology-level ones.
#'
#' \bold{Weighting source data}
#'
#' Two ways to do this. MRP is better as original data can be passed through pipeline plus MRP allows more complex forms of constraint than regular topology-level ones.
#'
#' \bold{Taxonomy-phylogeny contradictions}
#'
#' Partly a test that correct outgroups were used, but also way to potentially inform required updates to Paleobiology Database taxonomy.
#'
#' \bold{Outputs}
#'
#' Text.
#'
#' @return TBC.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Baum, B. R., 1992. Combining trees as a way of combining data sets for phylogenetic inference, and the desirability of combining gene trees. Taxon, 41, 3-10.
#'
#' Lloyd, G. T., Davis, K. E., Pisani, D., Tarver, J. E., Ruta, M., Sakamoto, M., Hone, D. W. E., Jennings, R. & Benton, M. J., 2008. Dinosaurs and the Cretaceous Terrestrial Revolution. Proceedings of the Royal Society B, 275, 2483-2490.
#'
#' Lloyd, G. T., Bapst, D. W., Friedman, M. and Davis, K. E., 2016. Probabilistic divergence time estimation without branch lengths: dating the origins of dinosaurs, avian flight, and crown birds. Biology Letters, 12, 20160609.
#'
#' Maddison, D. R., Swofford, D. L. and Maddison, W. P., 1997. NEXUS: an extensible file format for systematic information. Systematic Biology, 46, 590-621.
#'
#' Ragan, M., 1992. Phylogenetic inference based on matrix representation of trees. Molecular Phylogenetics and Evolution, 1, 113-126.
#'
#' @examples
#'
#' # Local test for Ichthyopterygia:
#' #Metatree(MRPDirectory = "/Users/eargtl/Documents/Homepage/www.graemetlloyd.com/mrp", XMLDirectory = "/Users/eargtl/Documents/Homepage/www.graemetlloyd.com/xml", TargetClade = "Ichthyopterygia", InclusiveDataList = sort(c(GetFilesForClade("matricht.html"), "Bickelmann_etal_2009a", "Caldwell_1996a", "Chen_etal_2014ba", "Chen_etal_2014bb", "deBraga_et_Rieppel_1997a", "Gauthier_etal_1988b", "Laurin_et_Reisz_1995a", "Muller_2004a", "Reisz_etal_2011a", "Rieppel_et_Reisz_1999a", "Rieppel_et_deBraga_1996a", "Young_2003a")), ExclusiveDataList = c("Averianov_inpressa", "Bravo_et_Gaete_2015a", "Brocklehurst_etal_2013a", "Brocklehurst_etal_2015aa", "Brocklehurst_etal_2015ab", "Brocklehurst_etal_2015ac", "Brocklehurst_etal_2015ad", "Brocklehurst_etal_2015ae", "Brocklehurst_etal_2015af", "Bronzati_etal_2012a", "Bronzati_etal_2015ab", "Brusatte_etal_2009ba", "Campbell_etal_2016ab", "Carr_et_Williamson_2004a", "Carr_etal_2017ab", "Frederickson_et_Tumarkin-Deratzian_2014aa", "Frederickson_et_Tumarkin-Deratzian_2014ab", "Frederickson_et_Tumarkin-Deratzian_2014ac", "Frederickson_et_Tumarkin-Deratzian_2014ad", "Garcia_etal_2006a", "Gatesy_etal_2004ab", "Grellet-Tinner_2006a", "Grellet-Tinner_et_Chiappe_2004a", "Grellet-Tinner_et_Makovicky_2006a", "Knoll_2008a", "Kurochkin_1996a", "Lopez-Martinez_et_Vicens_2012a", "Lu_etal_2014aa", "Norden_etal_inpressa", "Pisani_etal_2002a", "Ruiz-Omenaca_etal_1997a", "Ruta_etal_2003ba", "Ruta_etal_2003bb", "Ruta_etal_2007a", "Selles_et_Galobart_2016a", "Sereno_1993a", "Sidor_2001a", "Skutschas_etal_inpressa", "Tanaka_etal_2011a", "Toljagic_et_Butler_2013a", "Tsuihiji_etal_2011aa", "Varricchio_et_Jackson_2004a", "Vila_etal_2017a", "Wilson_2005aa", "Wilson_2005ab", "Zelenitsky_et_Therrien_2008a"), MissingSpecies = "exclude", BackboneConstraint = "Moon_inpressa", RelativeWeights = c(0, 100, 10, 1))
#'
#' @export Metatree
Metatree <- function(MRPDirectory, XMLDirectory, TargetClade = "", InclusiveDataList = c(), ExclusiveDataList = c(), HigherTaxaToCollapse = c(), SpeciesToExclude = c(), MissingSpecies = "exclude", Interval = NULL, VeilLine = TRUE, IncludeSpecimenLevelOTUs = TRUE, BackboneConstraint = NULL, MonophylyConstraint = NULL, RelativeWeights = c(1, 1, 1, 1), WeightCombination = "sum", ReportContradictionsToScreen = FALSE) {
  
  # DOUBLE CHECK PARENT REPLACEMENT LINE NOW MULTIPLE PARENTS EXIST - SEEMS TO WORK BUT MIGHT NOT.
  # Weights are also super slow (IntraMatrixWeights really?). Can this be sped up somehow? E.g., way STR is.
  
  # DEFO NEED TO IMPROVE NAME CHECKER AFTER RECON COS THAT IS WHERE MOST OF THE LATER ISSUES ARE (E.G. GENUS NAME NUMBER USED FOR SPECIES)
  # FOR HIGHER TAXA TO COLLAPSE HAVE TO ALSO EDIT CONSTRAINT TREES (AND CHECK THEY CAN EVEN MESH!)
  # CHECK FOR SPECIES THAT BELONG TO A GENUS DIFFERENT TO THE ONE IN THEIR NAME!
  # NEED TO CATCH ISSUE WHERE GENUS NUMBER IS USED FOR A SPECIES (HARD TO CHECK SO FAR DUE TO INDETERMINATES CONTINGENCY); I.E., A SAFER CHECK THAT RECON_NAME MATCHES DATABASE NAME FOR RECON_NO
  # NEED SOME TEST THAT HELPS DETERMINE IF MULTIPLE OCCURRENCES OF SAME TAXON AFTER RECONCILIATION IS CORRECT OR AN ERROR
  # CHECK THERE ARE MULTIPLE TAXA PRE-RECONCILIATION
  # CHECK INDETS DO NOT GIVE MULTIPLE MATCHES
  # ADD BLOCK NAMES TO OUTPUT (I.E. SOURCE DATA); COULD BUILD MRP MATRIX THIS WAY? (ADDITIONAL OUTPUT TYPE PROBABLY BEST WAY TO DO IT)
  # THINK ABOUT BETTER WASY TO HANDLE WEIGHTS AS CHARACTERS GET PRUNED/MERGED/DUPLICATED ETC.
  # CHECK INTERACTION BETWEEN COLLAPSED HIGHER TAXA AND EMPTY HIGHER TAXA THAT CURENTLY GET DELETED
  # ADD INCLUDED DATA SETS TO OUTPUT ALONGSIDE REMOVED!
  
  # OPTIONS TO ADD IN FUTURE:
  #
  # 1. Way to do historical metatrees (should be easy as just prune out starting data sets then run as normal - except setting current year to historical year).
  # 2. Allow Purvis coding instead of Baum and Ragan. (Involves dealing with NA issues that would be caused currently.)
  # 3. Proper chunking and maybe even terminal/TNT calls.
  # 4. Species to include option as alternative to species to exclude.
  # 5. Make Safe Taxonomic Reduction optional.
  # 6. Make operable as web site callsas well as/instead of local directory calls. (Needs proper XMLs to be available.)
  
  # Subfunction that gives just MRPs where matrix is still intact (has rows and columns):
  ActiveMRP <- function(MRPList) unname(which(unlist(lapply(MRPList, function(x) prod(dim(x$Matrix)))) > 0))
  
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
  
  # Subfunction to find contradicting MRP characters between a string (single character) and a matrix (multiple characters):
  MRPCharacterContradiction <- function(MRPCharacterString, MRPCharacterMatrix) {
    
    # Check MRP string has names and stop and warn user if not:
    if(is.null(names(MRPCharacterString))) stop("MRPCharacterString must have names. Add and try again.")
    
    # Check MRP matrix has names and stop and warn user if not:
    if(is.null(rownames(MRPCharacterMatrix))) stop("MRPCharacterMatrix must have row names. Add and try again.")
    
    # Check MRP string and matrix match in size and stop and warn user if not:
    if(length(MRPCharacterString) != nrow(MRPCharacterMatrix)) stop("MRPCharacterString must have the same length as the number of rows of MRPCharacterMatrix. Check data and try again.")
    
    # Check names of MRP string and matrix match and stop and warn user if not:
    if(!all(sort(names(MRPCharacterString)) == sort(rownames(MRPCharacterMatrix)))) stop("MRPCharacterString names must match row names of MRPCharacterMatrix. Check names and try again.")
    
    # Check only zeroes and ones are coded and stop and warn user if not:
    if(length(setdiff(unique(c(MRPCharacterString, MRPCharacterMatrix)), c("0", "1"))) > 0) stop("Both MRPCharacterString and MRPCharacterMatrix must consist exclusively of the characters \"0\" and \"1\". Check data and try again.")
    
    # Check MRP string contains both zeroe and ones and stop and warn user if not:
    if(length(unique(MRPCharacterString)) < 2) stop("MRPCharacterString must contain both zeroes and ones.")
    
    # Check every MRP matrix column contains both zeroe and ones and stop and warn user if not:
    if(length(unique(as.vector(MRPCharacterMatrix))) < 2) stop("MRPCharacterString must contain both zeroes and ones.")
    
    # Find names corresponding to scores of "0" in the MRP string:
    StringZeroNames <- names(which(MRPCharacterString == "0"))
    
    # Find names corresponding to scores of "1" in the MRP string:
    StringOneNames <- names(which(MRPCharacterString == "1"))
    
    # Find any matrix columns that contradict the string character (both "0" and "1" coded for zero and one matches in MRP string):
    ContradictionColumns <- which(apply(MRPCharacterMatrix, 2, function(x) length(unique(x[StringZeroNames])) == 2 && length(unique(x[StringOneNames])) == 2))
    
    # Return contradicting columns:
    return(ContradictionColumns)
    
  }
  
  # Subfunction to produce intra matrix weights:
  MRPIntraMatrixWeights <- function(MRPMatrix) {
    
    # Get a list of characters that contradict with each character in turn:
    ContradictionList <- lapply(apply(MRPMatrix, 2, as.list), function(x) MRPCharacterContradiction(unlist(x), MRPMatrix))
    
    # Add every character to its' own list:
    ContradictionList <- mapply(function(x, y) c(x, y), x = as.list(1:ncol(MRPMatrix)), y = ContradictionList)
    
    # Build vector of intra matrix weights:
    IntraMatrixWeights <- 1 / unlist(lapply(ContradictionList, function(x) length(unique(unlist(ContradictionList[x])))))
    
    # Return intra matrix weights:
    return(IntraMatrixWeights)
    
  }
  
  # Subfunction to collapse vector of string to single formatted string:
  WriteListAsString <- function(ListOfItems, OxfordComma = TRUE) {
    
    # If not using Oxford comma set final bridge as not having one:
    if(!OxfordComma) FinalBridge <- " and "
    
    # If using Oxford comma set final bridge as having one:
    if(OxfordComma) FinalBridge <- ", and "
    
    # If list is a single item make that the output:
    if(length(ListOfItems) == 1) Output <- ListOfItems
    
    # If list is two items join with a simple and:
    if(length(ListOfItems) == 2) Output <- paste(ListOfItems, collapse = " and ")
    
    # If three or more items format as list with Oxford Comma option included:
    if(length(ListOfItems) > 2) Output <- paste(paste(ListOfItems[1:(length(ListOfItems) - 2)], collapse = ", "), paste(ListOfItems[(length(ListOfItems) - 1):length(ListOfItems)], collapse = FinalBridge), sep = ", ")
    
    # Return output:
    return(Output)
    
  }
  
  # Subfunction to find taxonomy-phylogeny contradictions and turn them into warning messages:
  ListContradictions <- function(TaxonomyMRP, MRPMatrix, ContradictionTaxa, DataSetName) {
    
    # Subfunction to build warning messages:
    BuildWarningMessages <- function(BoundMatrix, ContradictionTaxon, DataSetName) {
      
      # Create four blocks for matching:
      ZeroZero <- ZeroOne <- OneZero <- OneOne <- BoundMatrix
      
      # Fill blocks with zeroes:
      ZeroOne[, 1] <- OneZero[, 2] <- ZeroZero[1:length(ZeroZero)] <- "0"
      
      # Fill blocks with ones:
      ZeroOne[, 2] <- OneZero[, 1] <- OneOne[1:length(OneOne)] <- "1"
      
      # Build empty list of names:
      NamesList <- list()
      
      # Store names where scores are zero and zero:
      NamesList[["ZeroZeroNames"]] <- names(which(apply(BoundMatrix == ZeroZero, 1, all)))
      
      # Store names where scores are zero and one:
      NamesList[["ZeroOneNames"]] <- names(which(apply(BoundMatrix == ZeroOne, 1, all)))
      
      # Store names where scores are one and zero:
      NamesList[["OneZeroNames"]] <- names(which(apply(BoundMatrix == OneZero, 1, all)))
      
      # Store names where scores are one and one:
      NamesList[["OneOneNames"]] <- names(which(apply(BoundMatrix == OneOne, 1, all)))
      
      # Prune down to just the minimum list lengths (the likley problem candidate(s)):
      NamesList <- NamesList[unlist(lapply(NamesList, length)) == min(unlist(lapply(NamesList, length)))]
      
      # Collapse each item of the list to a singel string:
      NamesList <- lapply(NamesList, WriteListAsString)
      
      # Find datasets with second case (apparent taxa outside clade that should be inside):
      CaseOne <- match(c("ZeroZeroNames", "OneZeroNames"), names(NamesList))
      
      # Find datasets with second case (apparent taxa inside clade that should be outside):
      CaseTwo <- match(c("ZeroOneNames", "OneOneNames"), names(NamesList))
      
      # Remove NAs from first case matches:
      CaseOne <- CaseOne[!is.na(CaseOne)]
      
      # Remove NAs from second case matches:
      CaseTwo <- CaseTwo[!is.na(CaseTwo)]
      
      # If first case exists then reformat string with message to user:
      if(length(CaseOne) > 0) NamesList[CaseOne] <- lapply(NamesList[CaseOne], function(x) paste("In ", DataSetName, " the following taxa were found outside ", ContradictionTaxon, " when taxonomy suggests they should be inside: ", x, ". Check data set and/or taxonomy that this is correct.\n", sep = ""))
      
      # If second case exists then reformat string with message to user:
      if(length(CaseTwo) > 0) NamesList[CaseTwo] <- lapply(NamesList[CaseTwo], function(x) paste("In ", DataSetName, " the following taxa were found inside ", ContradictionTaxon, " when taxonomy suggests they should be outside: ", x, ". Check data set and/or taxonomy that this is correct.\n", sep = ""))
      
      # Reformat names list as a vector for output:
      Output <- unname(unlist(NamesList))
      
      # If multiple values then reformat into a single value:
      if(length(Output) > 1) Output <- paste("In ", DataSetName, " one of the following is true:\n", paste(paste(paste(1:length(Output), ". T", sep = ""), gsub(paste("In ", DataSetName, " t| Check data set and/or taxonomy that this is correct.\n", sep = ""), "", Output), rep("\n", length(Output)), sep = ""), collapse = ""), "Check data set and/or taxonomy that this is correct.\n", sep = "")
      
      # Return output:
      return(Output)
      
    }
    
    # Build output into unique values (as can get duplicates):
    Output <- unique(unlist(lapply(as.list(ContradictionTaxa), function(x) lapply(as.list(MRPCharacterContradiction(TaxonomyMRP[, x], MRPMatrix)), function(y) {BoundMatrix <- cbind(TaxonomyMRP[rownames(MRPMatrix), x], MRPMatrix[, y]); BuildWarningMessages(BoundMatrix, x, DataSetName)}))))
    
    # Return output:
    return(Output)
    
  }

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
  
  # Set defualt of constraint in use to FALSE:
  ConstraintInUse <- FALSE
  
  # Check that there is a maximum of one constraint tree being used and stop and warn user if so.
  # Technically it ought to be possible to do this, but it leaves open some potentially horrendous disasters best avoided for now:
  if(!is.null(BackboneConstraint) && !is.null(MonophylyConstraint)) stop("Cannot currently apply a backbone constraint and a monophyly constraint simultaneously.")
  
  # If backbone constraint is set:
  if(!is.null(BackboneConstraint)) {
    
    # Set constraint in use to TRUE:
    ConstraintInUse <- TRUE
    
    # Check is only a single value and stop and warn user if not:
    if(length(BackboneConstraint) > 1) stop("BackboneConstraint must be a single value. (Cannot apply two backbone constraints simultaneously.)")
    
    # Check is a string and stop and warn user if not:
    if(!is.character(BackboneConstraint)) stop("BackboneConstraint must be a text string. Reformat and try again.")
    
    # Set constraint data set to backbone constraint:
    ConstraintDataSet <- BackboneConstraint
    
    # Set constraint type to backbone:
    ConstraintType <- "backbone"

  }
  
  # If monophyly onstraint is set:
  if(!is.null(MonophylyConstraint)) {
    
    # Set constraint in use to TRUE:
    ConstraintInUse <- TRUE

    # Check is only a single value and stop and warn user if not:
    if(length(MonophylyConstraint) > 1) stop("MonophylyConstraint must be a single value. (Cannot apply two monophyly constraints simultaneously.)")
    
    # Check is a string and stop and warn user if not:
    if(!is.character(MonophylyConstraint)) stop("MonophylyConstraint must be a text string. Reformat and try again.")
    
    # Set constraint data set to monophyly constraint:
    ConstraintDataSet <- MonophylyConstraint
    
    # Set constraint type to monophyly:
    ConstraintType <- "monophyly"
    
  }
  
  # Check VeilLine is a logical and stop and warn user if not:
  if(!is.logical(ReportContradictionsToScreen)) stop("ReportContradictionsToScreen must be a logical (TRUE or FALSE).")
  
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
  
  # Check there are four relative weights values and stop and warn user if not:
  if(length(RelativeWeights) != 4) stop("RelativeWeights must consist of exactly four values. Fix and try again.")
  
  # Check relative weights are numeric and stop and warn user if not:
  if(!is.numeric(RelativeWeights)) stop("RelativeWeights must consist of numeric values. Fix and try again.")
  
  # Check there are no negative weights and at least one positive weight is being used and stop and wanr user if not:
  if(sum(RelativeWeights) <= 0 || any(RelativeWeights < 0)) stop("RelativeWeights must include at least one positive value and negative values are not allowed. Fix and try again.")
  
  # Check only a single weight combination value is being used and stop and warn user if not:
  if(length(WeightCombination) != 1) stop("WeightCombination must consist of a single value.")
  
  # Check weight combintion is a valid option and stop and warn if not:
  if(length(setdiff(WeightCombination, c("product", "sum"))) > 0) stop("WeightCombination must be one of \"product\" or \"sum\" only. Check spelling and try again.")
  
  # Read in all MRP files and store in a list (include duplicate headers to store parent sibling info later):
  MRPList <- lapply(lapply(as.list(MRPFileList), Claddis::ReadMorphNexus), function(x) {y <- list(x$Matrix_1$Matrix, x$Matrix_1$Weights, "", "", ""); names(y) <- c("Matrix", "Weights", "FileName", "Parent", "Sibling"); y})
  
  # Set names of MRP files:
  names(MRPList) <- gsub("mrp.nex", "", MRPFileList)
  
  # Find maximum input weight:
  MaximumInputWeight <- max(unlist(lapply(MRPList, function(x) x$Weights)))
  
  # Resacle all input weights zero to one by dividing through by maximum input weight:
  MRPList <- lapply(MRPList, function(x) {x$Weights <- x$Weights / MaximumInputWeight; x})
  
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
  
  # Read in all XML files and store in a list (reformatting subgenera as GenusSubgenus along the way):
  XMLList <- lapply(as.list(XMLFileList), function(x) {y <- ReadMetatreeXML(x); y$SourceTree$Taxa$TagContents[, "recon_name"] <- gsub("_\\(|\\)", "", y$SourceTree$Taxa$TagContents[, "recon_name"]); y})
  
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
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- PisaniMRPPrune(Claddis::MakeMorphMatrix(x$Matrix, Weights = x$Weights, ignore.duplicate.taxa = TRUE)); x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})

  # Print current processing status:
  cat("Done\nSearching for and collapsing pre-reconciliation duplicated taxa...")
  
  # Collapse any duplicate taxon names:
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- Claddis::MakeMorphMatrix(x$Matrix, Weights = x$Weights, ignore.duplicate.taxa = TRUE); if(any(duplicated(rownames(y$Matrix_1$Matrix)))) {DuplicateNames <- setdiff(unlist(lapply(strsplit(rownames(y$Matrix_1$Matrix)[duplicated(rownames(y$Matrix_1$Matrix))], split = "%%%%"), '[', 2)), "DELETE"); if(length(DuplicateNames) > 0) cat(paste("\nDuplicate resolved OTU name(s) found in ", x$FileName, ": ", paste(DuplicateNames, collapse = ", "), ". Check this is correct.", sep = "")); y <- CollapseDuplicateTaxonMRP(y)}; x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
  
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
    MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- PisaniMRPPrune(Claddis::MakeMorphMatrix(x$Matrix, Weights = x$Weights, ignore.duplicate.taxa = TRUE)); x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
    
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
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- Claddis::MakeMorphMatrix(x$Matrix, Weights = x$Weights, ignore.duplicate.taxa = TRUE); if(any(duplicated(rownames(y$Matrix_1$Matrix)))) {DuplicateNames <- rownames(y$Matrix_1$Matrix)[duplicated(rownames(y$Matrix_1$Matrix))]; y <- CollapseDuplicateTaxonMRP(y)}; x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
  
  # Prune characters made redundant by these collapses:
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- PisaniMRPPrune(Claddis::MakeMorphMatrix(x$Matrix, Weights = x$Weights, ignore.duplicate.taxa = TRUE)); x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
  
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
  ParentChildRelationships <- paste(unlist(lapply(strsplit(ValidOTUNames, "%%%%"), '[', 1)), ResolvedTaxonNumbers[match(unlist(lapply(strsplit(ValidOTUNames, "%%%%"), '[', 1)), ResolvedTaxonNumbers[, "ResolvedTaxonNo"]), "ParentTaxonNo"], sep = " belongs to ")
  
  # If including specimen level OTUs:
  if(IncludeSpecimenLevelOTUs) {
    
    # Find which rows correspond to indeterminate and sp taxa (i.e., those where parent should be initial reconciliation):
    indetsandsps <- sort(c(which(unlist(lapply(lapply(lapply(lapply(lapply(strsplit(ValidOTUNames, "%%%%"), '[', 2), strsplit, split = "_"), unlist), '==', "indet"), any))), which(unlist(lapply(lapply(lapply(lapply(lapply(strsplit(ValidOTUNames, "%%%%"), '[', 2), strsplit, split = "_"), unlist), '==', "sp"), any)))))
    
    # If such taxa exist then update parent child relationships accordingly:
    if(length(indetsandsps) > 0) ParentChildRelationships[indetsandsps] <- paste(unlist(lapply(strsplit(ValidOTUNames[indetsandsps], "%%%%"), '[', 1)), unlist(lapply(strsplit(ValidOTUNames[indetsandsps], "%%%%"), '[', 1)), sep = " belongs to ")
    
  }
  
  # Get list of new children (for which parents are needed) - excludes "Life" which has no parent:
  newchildren <- setdiff(unlist(lapply(strsplit(ParentChildRelationships, " belongs to "), '[', 2)), "28595")
  
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
    ParentChildRelationships <- c(ParentChildRelationships, paste(ResolvedTaxonNumbers[match(newchildren, ResolvedTaxonNumbers[, "ResolvedTaxonNo"]), "ResolvedTaxonNo"], ResolvedTaxonNumbers[match(newchildren, ResolvedTaxonNumbers[, "ResolvedTaxonNo"]), "ParentTaxonNo"], sep = " belongs to "))
    
    # Update new children:
    newchildren <- setdiff(ResolvedTaxonNumbers[match(newchildren, ResolvedTaxonNumbers[, "ResolvedTaxonNo"]), "ParentTaxonNo"], "28595")
    
  }
  
  # If Life is missing then add it at bottom:
  if(all(!ResolvedTaxonNumbers[, "ResolvedTaxonNo"] == "28595")) ResolvedTaxonNumbers <- rbind(ResolvedTaxonNumbers, c("28595", "Life", "25", NA))
  
  # Convert parent-child relationships into a matrix (columns for child and parent):
  parentchildmatrix <- matrix(unlist(strsplit(ParentChildRelationships, split = " belongs to ")), ncol = 2, byrow = TRUE, dimnames = list(c(), c("Child", "Parent")))
  
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
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- Claddis::MakeMorphMatrix(x$Matrix, Weights = x$Weights, ignore.duplicate.taxa = TRUE); if(any(duplicated(rownames(y$Matrix_1$Matrix)))) {DuplicateNames <- rownames(y$Matrix_1$Matrix)[duplicated(rownames(y$Matrix_1$Matrix))]; if(length(DuplicateNames) > 0) cat(paste("\nDuplicate resolved OTU name(s) found post higher-taxon substitution in ", x$FileName, ": ", paste(DuplicateNames, collapse = ", "), ". Check this is correct.", sep = "")); y <- CollapseDuplicateTaxonMRP(y)}; x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
  
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
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- PisaniMRPPrune(Claddis::MakeMorphMatrix(x$Matrix, Weights = x$Weights, ignore.duplicate.taxa = TRUE)); x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
  
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
    
    # Add uppercase rownames to new block:
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
    MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- Claddis::MakeMorphMatrix(x$Matrix, Weights = x$Weights, ignore.duplicate.taxa = TRUE); if(any(duplicated(rownames(y$Matrix_1$Matrix)))) {DuplicateNames <- rownames(y$Matrix_1$Matrix)[duplicated(rownames(y$Matrix_1$Matrix))]; y <- CollapseDuplicateTaxonMRP(y)}; x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
    
    # Prune constant characters and collapse duplicated characters:
    MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- PisaniMRPPrune(Claddis::MakeMorphMatrix(x$Matrix, Weights = x$Weights, ignore.duplicate.taxa = TRUE)); x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
    
    # Update new valid OTUs:
    NewValidOTUs <- sort(rownames(TaxonomyMRP))
    
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
  
  # Set current year (used multiple times later):
  CurrentYear <- strsplit(as.character(Sys.Date()), "-")[[1]][1]
  
  # Add publication year to each data set (in presses are ascribed the current year):
  MRPList <- lapply(MRPList, function(x) {x$PublicationYear <- gsub("[:A-Z:a-z:]|_|-", "", gsub("inpress", CurrentYear, x$FileName)); x})
  
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
    
    # Start with current year as veil year:
    CurrentVeilYear <- as.numeric(CurrentYear)
    
    # Set current taxa included as being from current veil year (to present):
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
    
    # Find any data sets to remove (older than veil year):
    DataSetsToRemove <- which(as.numeric(unlist(lapply(MRPList, function(x) x$PublicationYear))) < CurrentVeilYear)
    
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
  
  # Print current processing status:
  cat("Done\nCalculating weights...")
  
  # Reformat weights to be input weights, publication year weights (equation 1 in Supplementary Information of Lloyd et al. 2016),
  # data set dependence weights (1 / (N siblings + 1)), and clade contradiction weights (1 / frequency of contradictory clades).
  # All weights are set on a zero to one scale initially and then multiplied by RelativeWeights:
  MRPList <- lapply(MRPList, function(x) {InputWeights <- x$Weights; x$Weights <- NULL; x$InputWeights <- RelativeWeights[1] * InputWeights; x$PublicationYearWeights <- RelativeWeights[2] * (rep(((2 ^ (0.5 * (as.numeric(x$PublicationYear) - CurrentVeilYear + 1))) - 1) / ((2 ^ (0.5 * (as.numeric(CurrentYear) - CurrentVeilYear + 1))) - 1), length(InputWeights))); x$DataSetDependenceWeights <- RelativeWeights[3] * rep(1 / (length(x$Sibling) + 1), length(InputWeights)); x$CladeContradictionWeights <- RelativeWeights[4] * MRPIntraMatrixWeights(x$Matrix); x})
  
  # If using sum to combine weights collapse weights to just their sum:
  if(WeightCombination == "sum") MRPList <- lapply(MRPList, function(x) {x$Weights <- apply(rbind(x$InputWeights, x$PublicationYearWeights, x$DataSetDependenceWeights, x$CladeContradictionWeights), 2, sum); x$InputWeights <- NULL; x$PublicationYearWeights <- NULL; x$DataSetDependenceWeights <- NULL; x$CladeContradictionWeights <- NULL; x})
  
  # If using product to combine weights collapse weights to just their product (excluding zeroes) and remove other weights from list:
  if(WeightCombination == "product") MRPList <- lapply(MRPList, function(x) {x$Weights <- apply(rbind(x$InputWeights, x$PublicationYearWeights, x$DataSetDependenceWeights, x$CladeContradictionWeights), 2, function(y) {z <- as.character(y); z[z == "0"] <- NA; prod(as.numeric(z), na.rm = TRUE)}); x$InputWeights <- NULL; x$PublicationYearWeights <- NULL; x$DataSetDependenceWeights <- NULL; x$CladeContradictionWeights <- NULL; x})
  
  # Get current maximum weight:
  MaximumWeight <- max(unlist(lapply(MRPList, function(x) x$Weights)))
  
  # Get current minimum weight:
  MinimumWeight <- min(unlist(lapply(MRPList, function(x) x$Weights)))
  
  # Calculate the multiplication factor for weight rescaling (10 to 1000):
  MultiplicationFactor <- 1 / ((MaximumWeight - MinimumWeight) / 990)
  
  # Calculate the addition factor for weight rescaling (10 to 1000):
  AdditionFactor <- 10 - (MinimumWeight * MultiplicationFactor)
  
  # Rescale weights (10 to 1000) and round results to two decimal places (best TNT can cope with):
  MRPList <- lapply(MRPList, function(x) {x$Weights <- round((x$Weights * MultiplicationFactor) + AdditionFactor, 2); x})
  
  # Print current processing status:
  cat("Done\nChecking for phylogeny-taxonomy contradictions...")
  
  # Do first pass to find any contradictions between taxonomy MRP and the fully reconciled source matrix:
  MRPList <- lapply(MRPList, function(x) {TaxonomyMRPStrings <- TaxonomyMRP[rownames(x$Matrix), ]; TaxonomyMRPStrings[is.na(TaxonomyMRPStrings)] <- "0"; TaxonomyMRPStrings <- TaxonomyMRPStrings[, apply(TaxonomyMRPStrings, 2, function(y) length(unique(y))) == 2]; x$TaxonomyContradictions <- names(which(unlist(lapply(apply(TaxonomyMRPStrings, 2, list), function(z) length(MRPCharacterContradiction(unlist(z), x$Matrix))) > 0))); x$TaxonomyContradictionProportion <- length(x$TaxonomyContradictions) / ncol(TaxonomyMRPStrings); x})
  
  # Store monophyletic taxa (those not contradicted by any phylogenetic characters - useful for chunking larger data if found)
  MonophyleticTaxa <- setdiff(colnames(TaxonomyMRP), unique(unlist(lapply(MRPList, function(x) x$TaxonomyContradictions))))
  
  # If reporting contradiction issues to the screen:
  if(ReportContradictionsToScreen) {
    
    # Find any data sets with taxonomy-phylogeny contradictions:
    ContradictionIssueDataSets <- names(unlist(lapply(MRPList, function(x) length(x$TaxonomyContradictions))) > 0)
    
    # If there are data sets with contradictions:
    if(length(ContradictionIssueDataSets) > 0) {
      
      # Build vector of contradiction warnings:
      ContradictionWarnings <- unname(unlist(lapply(MRPList[ContradictionIssueDataSets], function(x) {TaxonomyMRPSubset <- TaxonomyMRP[rownames(x$Matrix), x$TaxonomyContradictions, drop = FALSE]; if(any(is.na(TaxonomyMRPSubset))) TaxonomyMRPSubset[is.na(TaxonomyMRPSubset)] <- "0"; ListContradictions(TaxonomyMRP = TaxonomyMRPSubset, MRPMatrix = x$Matrix, ContradictionTaxa = x$TaxonomyContradictions, DataSetName = x$FileName)})))
      
      # Print warnings to screen:
      cat(ContradictionWarnings)
      
      # NEED TO BREAK THIS DOWN FURTHER AS CLEARLY SOME REDUNDANCY! (E.G. GROUPING HIGHER TAXA WITH SAME ISSUE, OR DATA SETS WITH SAME ISSUE)
      
    }
    
  }
  
  # If using a constraint:
  if(ConstraintInUse) {
    
    # Print current processing status:
    cat("Done\nApply constraint tree...")
    
    # If a monophyly constraint add all other taxa outside the constraint:
    if(ConstraintType == "monophyly") MRPList[[grep(ConstraintDataSet, names(MRPList))]]$Matrix <- rbind(MRPList[[grep(ConstraintDataSet, names(MRPList))]]$Matrix, matrix("0", ncol = ncol(MRPList[[grep("Constraint", names(MRPList))]]$Matrix), nrow = length(setdiff(NewValidOTUs, rownames(MRPList[[grep(ConstraintDataSet, names(MRPList))]]$Matrix))), dimnames = list(setdiff(NewValidOTUs, rownames(MRPList[[grep(ConstraintDataSet, names(MRPList))]]$Matrix)), c())))
    
    # Get combined weight of all non-constraint data (need to know to correctly weight the constraint data):
    NonConstraintWeightsTotal <- sum(unname(unlist(lapply(MRPList[-grep(ConstraintDataSet, names(MRPList))], function(x) x$Weights))))
    
    # Embiggen MRP matrix so that weights are high enough to ensure contstraint gets implemented:
    MRPList[[ConstraintDataSet]]$Matrix <- metatree::EmbiggenMatrix(Claddis::MakeMorphMatrix(MRPList[[ConstraintDataSet]]$Matrix, Weights = MRPList[[ConstraintDataSet]]$Weights), N = ceiling(NonConstraintWeightsTotal / 1000))$Matrix_1$Matrix
    
    # Ensure current constraint weights are set at maximum possible value (1000):
    MRPList[[ConstraintDataSet]]$Weights <- rep(1000, ncol(MRPList[[ConstraintDataSet]]$Matrix))
    
  }

  # Print current processing status:
  cat("Done\nBuilding MRP matrix...")
  
  # Add in missing taxa as NAs to every taxon:
  MRPList <- lapply(MRPList, function(x) {MissingTaxa <- setdiff(rownames(TaxonomyMRP), rownames(x$Matrix)); if(length(MissingTaxa) > 0) x$Matrix <- rbind(x$Matrix, matrix(nrow = length(MissingTaxa), ncol = ncol(x$Matrix), dimnames = list(MissingTaxa, c()))); x$Matrix <- x$Matrix[rownames(TaxonomyMRP), , drop = FALSE]; x})
  
  # Build full MRP matrix:
  FullMRPMatrix <- MakeMorphMatrix(CharacterTaxonMatrix = cbind(do.call(cbind, lapply(MRPList, function(x) x$Matrix)), TaxonomyMRP), Weights = c(unname(unlist(lapply(MRPList, function(x) x$Weights))), rep(1, ncol(TaxonomyMRP))))
  
  # Add all zero outgroup to matrix:
  FullMRPMatrix$Matrix_1$Matrix <- rbind(matrix("0", nrow = 1, ncol = ncol(FullMRPMatrix$Matrix_1$Matrix), dimnames = list("allzero", c())), FullMRPMatrix$Matrix_1$Matrix)

  # Print current processing status:
  cat("Done\nPerforming Safe Taxonomic Reduction...")
  
  # Perform STR on full matrix:
  STRData <- SafeTaxonomicReduction(FullMRPMatrix)
  
  # Create additional STR matrix from full matrix:
  STRMRPMatrix <- STRData$reduced.matrix
  
  # Print current processing status:
  cat("Done\nCompiling and returning output...")
  
  # Compile output:
  Output <- list(FullMRPMatrix, STRMRPMatrix, TaxonomyMRPTree, STRData$str.list, RemovedSourceData, CurrentVeilYear)
  
  # Add names:
  names(Output) <- c("FullMRPMatrix", "STRMRPMatrix", "TaxonomyTree", "SafelyRemovedTaxa", "RemovedSourceData", "CurrentVeilYear")
  
  # Print current processing status:
  cat("Done")
  
  # Return output:
  return(Output)
  
}

#MRPDirectory <- "/Users/eargtl/Documents/Homepage/www.graemetlloyd.com/mrp" # MRP file directory
#XMLDirectory <- "/Users/eargtl/Documents/Homepage/www.graemetlloyd.com/xml" # XML file directory
#TargetClade <- "Cetacea"
#InclusiveDataList <- c("Aguirre-Fernandez_etal_2009a", "Aguirre-Fernandez_et_Fordyce_2014a", "Albright_etal_inpressa", "Arnold_etal_2005a", "Benoit_etal_2011a", "Bianucci_2013a", "Bianucci_et_Gingerich_2011a", "Bianucci_et_Landini_2006a", "Bianucci_etal_2007a", "Bianucci_etal_2010a", "Bianucci_etal_2013a", "Bianucci_etal_2016a", "Bianucci_etal_2018aa", "Bianucci_etal_2018ab", "Bisconti_2008a", "Bisconti_2010a", "Bisconti_et_Bosselaers_2016a", "Bisconti_etal_2013a", "Bisconti_etal_2019a", "Bisconti_inpressa", "Boersma_et_Pyenson_2015a", "Boersma_et_Pyenson_2016a", "Boersma_etal_2017a", "Boessenecker_et_Fordyce_2015a", "Boessenecker_et_Fordyce_2017a", "Boessenecker_et_Fordyce_inpressa", "Boessenecker_etal_2017a", "Bosselaers_et_Post_2010a", "Bouetel_et_de_Muizon_2006a", "Buono_et_Cozzuol_2013a", "Buono_etal_2017a", "Churchill_etal_2012a", "Churchill_etal_2016a", "Collareta_etal_2017a", "Colpaert_etal_inpressa", "Demere_etal_2008a", "Dooley_etal_2004a", "Ekdale_etal_2011a", "El_Adli_etal_2014a", "Fajardo-Mellor_etal_2006a", "Fitzgerald_2010a", "Fordyce_1994a", "Fordyce_et_Marx_2013a", "Fordyce_et_Marx_2016a", "Fordyce_et_Marx_2018a", "Geisler_2001a", "Geisler_et_Luo_1996a", "Geisler_et_Sanders_2003a", "Geisler_et_Theodor_2009a", "Geisler_et_Uhen_2003a", "Geisler_et_Uhen_2005a", "Geisler_etal_2005a", "Geisler_etal_2011a", "Geisler_etal_2012a", "Geisler_etal_2014a", "Geisler_etal_2017a", "Gibson_etal_inpressa", "Godfrey_etal_2016a", "Godfrey_etal_2017a", "Goldin_2018a", "Goldin_et_Startsev_2014a", "Goldin_et_Startsev_2017a", "Goldin_et_Steeman_2015a", "Goldin_et_Zvonok_2013a", "Goldin_etal_2014a", "Goldin_etal_inpressa", "Hernandez_Cisneros_2018a", "Heyning_1997a", "Kimura_et_Hasegawa_2010a", "Lambert_2008a", "Lambert_2008b", "Lambert_et_Louwye_2016a", "Lambert_etal_2008a", "Lambert_etal_2009a", "Lambert_etal_2010a", "Lambert_etal_2013a", "Lambert_etal_2014a", "Lambert_etal_2015a", "Lambert_etal_2017a", "Lambert_etal_2017b", "Lambert_etal_inpressa", "Lambert_etal_inpressb", "Lambert_etal_inpressc", "Lambert_etal_inpressda", "Lambert_etal_inpressdb", "Luo_et_Marsh_1996a", "Martinez-Caceres_etal_2017a", "Marx_2011a", "Marx_et_Fordyce_2015a", "Marx_et_Fordyce_2016a", "Marx_etal_2015a", "Marx_etal_2016a", "Marx_etal_2017a", "McGowen_etal_2009a", "Messenger_et_McGuire_1998a", "Mijan_etal_inpressa", "Muizon_etal_2019a", "Murakami_etal_2012a", "Murakami_etal_2012b", "Murakami_etal_2014a", "Murakami_etal_2014b", "Murakami_etal_inpressa", "Nelson_et_Uhen_inpressa", "OLeary_et_Gatesy_2008a", "Paolucci_etal_inpressa", "Peredo_et_Pyenson_2018a", "Peredo_et_Uhen_2016a", "Peredo_etal_inpressa", "Post_etal_2017a", "Pyenson_etal_2015a", "Racicot_etal_2019aa", "Racicot_etal_2019ab", "Ramassamy_2016a", "Sanders_et_Geisler_inpressa", "Solis-Anorve_etal_2019a", "Spaulding_etal_2009a", "Steeman_2007a", "Tanaka_et_Fordyce_2014a", "Tanaka_et_Fordyce_2015a", "Tanaka_et_Fordyce_2016a", "Tanaka_et_Fordyce_inpressa", "Tanaka_et_Watanabe_inpressa", "Tanaka_etal_2018a", "Theodor_et_Foss_2005a", "Thewissen_etal_2001a", "Thewissen_etal_2007a", "Tsai_et_Fordyce_2018a", "Tsai_et_Fordyce_inpressa", "Tsai_et_Fordyce_inpressb", "Uhen_1999a", "Uhen_2004a", "Uhen_2014a", "Uhen_et_Gingerich_2001a", "Velez-Juarbe_etal_2015a", "Velez-Juarbe_inpressa", "Viglino_etal_2019a", "Viglino_etal_inpressa", "Wichura_etal_inpressa")
#ExclusiveDataList <- c("Averianov_inpressa", "Bravo_et_Gaete_2015a", "Brocklehurst_etal_2013a", "Brocklehurst_etal_2015aa", "Brocklehurst_etal_2015ab", "Brocklehurst_etal_2015ac", "Brocklehurst_etal_2015ad", "Brocklehurst_etal_2015ae", "Brocklehurst_etal_2015af", "Bronzati_etal_2012a", "Bronzati_etal_2015ab", "Brusatte_etal_2009ba", "Campbell_etal_2016ab", "Carr_et_Williamson_2004a", "Carr_etal_2017ab", "Frederickson_et_Tumarkin-Deratzian_2014aa", "Frederickson_et_Tumarkin-Deratzian_2014ab", "Frederickson_et_Tumarkin-Deratzian_2014ac", "Frederickson_et_Tumarkin-Deratzian_2014ad", "Garcia_etal_2006a", "Gatesy_etal_2004ab", "Grellet-Tinner_2006a", "Grellet-Tinner_et_Chiappe_2004a", "Grellet-Tinner_et_Makovicky_2006a", "Knoll_2008a", "Kurochkin_1996a", "Lopez-Martinez_et_Vicens_2012a", "Lu_etal_2014aa", "Norden_etal_inpressa", "Pisani_etal_2002a", "Ruiz-Omenaca_etal_1997a", "Ruta_etal_2003ba", "Ruta_etal_2003bb", "Ruta_etal_2007a", "Selles_et_Galobart_2016a", "Sereno_1993a", "Sidor_2001a", "Skutschas_etal_inpressa", "Tanaka_etal_2011a", "Toljagic_et_Butler_2013a", "Tsuihiji_etal_2011aa", "Varricchio_et_Jackson_2004a", "Vila_etal_2017a", "Wilson_2005aa", "Wilson_2005ab", "Zelenitsky_et_Therrien_2008a")
#HigherTaxaToCollapse = c()
#SpeciesToExclude = c()
#MissingSpecies = "exclude"
#Interval = NULL
#VeilLine = TRUE
#IncludeSpecimenLevelOTUs = TRUE
#BackboneConstraint = "McGowen_etal_2009a"
#MonophylyConstraint = NULL
#RelativeWeights = c(0, 100, 10, 1)
#WeightCombination = "sum"
#ReportContradictionsToScreen = FALSE

#MRPDirectory <- "~/Dropbox/Mammal_Supertree/Primates/GraemeVersion/MRP" # MRP file directory
#XMLDirectory <- "~/Dropbox/Mammal_Supertree/Primates/GraemeVersion/XML" # XML file directory
#TargetClade <- "Primates"
#InclusiveDataList <- c()
#ExclusiveDataList <- c("Averianov_inpressa", "Bravo_et_Gaete_2015a", "Brocklehurst_etal_2013a", "Brocklehurst_etal_2015aa", "Brocklehurst_etal_2015ab", "Brocklehurst_etal_2015ac", "Brocklehurst_etal_2015ad", "Brocklehurst_etal_2015ae", "Brocklehurst_etal_2015af", "Bronzati_etal_2012a", "Bronzati_etal_2015ab", "Brusatte_etal_2009ba", "Campbell_etal_2016ab", "Carr_et_Williamson_2004a", "Carr_etal_2017ab", "Frederickson_et_Tumarkin-Deratzian_2014aa", "Frederickson_et_Tumarkin-Deratzian_2014ab", "Frederickson_et_Tumarkin-Deratzian_2014ac", "Frederickson_et_Tumarkin-Deratzian_2014ad", "Garcia_etal_2006a", "Gatesy_etal_2004ab", "Grellet-Tinner_2006a", "Grellet-Tinner_et_Chiappe_2004a", "Grellet-Tinner_et_Makovicky_2006a", "Knoll_2008a", "Kurochkin_1996a", "Lopez-Martinez_et_Vicens_2012a", "Lu_etal_2014aa", "Norden_etal_inpressa", "Pisani_etal_2002a", "Ruiz-Omenaca_etal_1997a", "Ruta_etal_2003ba", "Ruta_etal_2003bb", "Ruta_etal_2007a", "Selles_et_Galobart_2016a", "Sereno_1993a", "Sidor_2001a", "Skutschas_etal_inpressa", "Tanaka_etal_2011a", "Toljagic_et_Butler_2013a", "Tsuihiji_etal_2011aa", "Varricchio_et_Jackson_2004a", "Vila_etal_2017a", "Wilson_2005aa", "Wilson_2005ab", "Zelenitsky_et_Therrien_2008a")
#HigherTaxaToCollapse = c()
#SpeciesToExclude = c()
#MissingSpecies = "exclude"
#Interval = NULL
#VeilLine = TRUE
#IncludeSpecimenLevelOTUs = TRUE
#BackboneConstraint = NULL
#MonophylyConstraint = NULL
#RelativeWeights = c(0, 100, 10, 1)
#WeightCombination = "sum"
#ReportContradictionsToScreen = FALSE

#Primates <- Metatree(MRPDirectory = MRPDirectory, XMLDirectory = XMLDirectory, TargetClade = TargetClade, InclusiveDataList = c(), ExclusiveDataList = ExclusiveDataList, HigherTaxaToCollapse = HigherTaxaToCollapse, SpeciesToExclude = SpeciesToExclude, MissingSpecies = MissingSpecies, Interval = Interval, VeilLine = VeilLine, IncludeSpecimenLevelOTUs = IncludeSpecimenLevelOTUs, BackboneConstraint = BackboneConstraint, MonophylyConstraint = MonophylyConstraint, RelativeWeights = RelativeWeights, WeightCombination = WeightCombination, ReportContradictionsToScreen = ReportContradictionsToScreen)
