#' Palaeobiology Database Tree Builder
#'
#' @description
#'
#' Using the Paleobiology Database taxonomy builds a phylogenetic tree.
#'
#' @param taxon_nos Either a vector of Paleobiology Database taxon numbers that will serve as tips, or a single number that will define the clade requested. Number(s) must match Paleobiology Database records.
#' @param taxon_names Either a vector of Paleobiology Database taxon names that will serve as tips, or a single name that will define the clade requested. Name(s) must match Paleobiology Database records.
#' @param original Option to be passed to \link{PaleobiologyDBDescendantFinder} or \link{PaleobiologyDBTaxaQuerier}.
#' @param interval Option to be passed to \link{PaleobiologyDBDescendantFinder} or \link{PaleobiologyDBTaxaQuerier}.
#' @param extant Option to be passed to \link{PaleobiologyDBDescendantFinder} or \link{PaleobiologyDBTaxaQuerier}.
#' @param stopfororphans Option to be passed to \link{PaleobiologyDBTaxaQuerier}.
#' @param validonly Option to be passed to \link{PaleobiologyDBDescendantFinder}.
#' @param returnrank Option to be passed to \link{PaleobiologyDBDescendantFinder}. Default is "3" (species level).
#' @param breaker Option to be passed to \link{PaleobiologyDBDescendantFinder} or \link{PaleobiologyDBTaxaQuerier}.
#' @param plot.tree Logical whether or not to produce a plot of the resulting tree alongside the output (default is FALSE).
#' @param TimeScale Logical indicating whether or not to timescale the tree. (NOT OPERATIONAL YET!)
#'
#' @details
#'
#' Taxonomies such as those in the Paleobiology Database (\code{paleobiodb.org}; queriable via the API, Peters and McLennen 2016) can also be represented as phylogenetic trees. These can then be used either to aid in reconstructing meta-analytical phylogenies (e.g., Lloyd et al. 2016), or as surrogates for phylogenies where none exist (Soul and Friedman 2015).
#'
#' The function presented here can either be handed a single higher taxon (e.g., Dinosauria) or a series of (presumably) lower level taxa (genera, species) for which a taxonomic "phylogeny" is desired. Note that all names must appear in the Paleobiology Database, and that it is recommended you actual use the taxonomic number desired to avoid any potential homonym issues, either between animals and plants or lower- and higher-level taxa.
#'
#' Internally the function will call \link{PaleobiologyDBDescendantFinder} and/or \link{PaleobiologyDBTaxaQuerier} and most options refer to these functions. Note that this means the function will typically take several seconds to run so do not expect an immediate results, even if the desired clade is very small.
#'
#' The resulting tree is in ape format with node labels, but currently no ages or branch lengths (these are planned future additions).
#'
#' Note also that the resulting tree should be treated with caution as it will rely solely on what is currently accessible in the Paleobiology Database (\code{paleobiodb.org}). As this database is incomplete it will not necessarily include all valid species (i.e., it logically cannot include taxa not entered in the database). Similarly, the resulting taxonomy is a synthesis of the set of opinions currently entered into the database for the taxa concerned.
#'
#' Finally, note that typically the resulting tree will contain multiple multifurcations (polytomies) and hence may not be usable in a lot of phylogenetic applications without first resolving these to form a bifurcating tree. (Interested users should consult the functions \code{paleotree::resolveTreeChar} or \code{paleotree::timeLadderTree}.) Caution should also be applied in how this is done, e.g., see Bell and Lloyd (2015).
#'
#' @return
#'
#' A phylo object in ape format with node labels. Note that sometimes multiple labels may be valid for a node and if so these are separated by "_et_".
#'
#' @author
#'
#' Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Bell, M. A. and Lloyd, G. T., 2015. strap: an R package for plotting phylogenies against stratigraphy and assessing their stratigraphic congruence. \emph{Palaeontology}, \bold{58}, 379-389.
#'
#' Lloyd, G. T., Bapst, D. W., Friedman, M. and Davis, K. E., 2016. Probabilistic divergence time estimation without branch lengths: dating the origins of dinosaurs, avian flight, and crown birds. \emph{Biology Letters}, \bold{12}, 20160609.
#'
#' Peters, S. E. and McClennen, M., 2016. The Paleobiology Database application programming interface. \emph{Paleobiology}, \bold{42}, 1-7.
#'
#' Soul, L. C. and Friedman, M., 2015, Taxonomy and phylogeny can yield comparable results in comparative paleontological analyses. \emph{Systematic Biology}, \bold{64}, 608-620.
#'
#' @seealso
#'
#' See also the \code{makePBDBtaxonTree} function in the \code{paleotree} package.
#'
#' @examples
#'
#' # Build a phylogenetic tree out of some famous dinosaurs:
#' PaleobiologyDBTreeBuilder(taxon_names =
#'   c("Tyrannosaurus rex", "Triceratops horridus",
#'   "Brontosaurus excelsus"), plot.tree = TRUE)
#'
#' @export PaleobiologyDBTreeBuilder
PaleobiologyDBTreeBuilder <- function(taxon_nos = NULL, taxon_names = NULL, original = FALSE, interval = NULL, extant = "include", stopfororphans = TRUE, validonly = TRUE, returnrank = "3", breaker = 100, plot.tree = FALSE, TimeScale = FALSE) {

  # FOR FUTURE ADD TIMESCALING!
  
  # ROOTING? ARE TREES ROOTED BY TIME OR WHAT?
  # NEED CHECKS NAMES ARE VALID AND THAT RANKS ARE AT LEAST APPROXIMATELY EQUIVALENT (SPECIES OR GENERA?)
  # SINGLE TAXON INPUT MEANS FIND ALL TAXA ASSIGNED TO IT MULTIPLE MEANS USE THEM AS INPUT?
  # FOR MULTIPLE CLADE LABELS OFFER OPTION TO ONLY USE LEAST INCLUSIVE (LOWEST LEVEL) ONE (OR JUST DO THIS)
  
  #taxon_nos = NULL
  #taxon_names = "Dinosauria"
  #original = TRUE
  #interval = c("Triassic", "Triassic")
  #extant = "include"
  #stopfororphans = TRUE
  #validonly = TRUE
  #returnrank = "3"
  #breaker = 100
  #plot.tree = FALSE
  #TimeScale = FALSE
  
  # Check at least one of numbers or names has been set:
  if(is.null(taxon_nos) && is.null(taxon_names)) stop("Must define at least one taxon name or taxon number")
  
  # If both numbers and names have been set:
  if(!is.null(taxon_nos) && !is.null(taxon_names)) {
    
    # Set names to NULL:
    taxon_names <- NULL
    
    # Tell user that this is what has happened:
    cat("Both number(s) and name(s) have been set. Proceeding with only number(s).")
    
  }
  
  # If using taxon names:
  if(!is.null(taxon_names)) {
    
    # If a single taxon name (meaning all children are being requested:
    if(length(taxon_names) == 1) {
      
      # Set taxon numbers as descendants of taxon name:
      taxon_nos <- gsub("txn:|var:", "", unname(unlist(lapply(apply(PaleobiologyDBDescendantFinder(taxon_nos = "1", taxon_names = taxon_names, original = original, interval = interval, extant = extant, validonly = validonly, returnrank = returnrank, breaker = breaker)[, c("OriginalTaxonNo", "ResolvedTaxonNo")], 1, list), function(x) unlist(x)[!is.na(unlist(x))][1]))))
      
      # Now nullify taxon names:
      taxon_names <- NULL
      
    # If multiple taxon names:
    } else {
      
      # Set numbers for names:
      taxon_nos <- gsub("txn:|var:", "", unname(unlist(lapply(apply(PaleobiologyDBTaxaQuerier(taxon_nos = as.character(1:length(taxon_names)), taxon_names = taxon_names, original = original, interval = interval, extant = extant, stopfororphans = TRUE, breaker = breaker)[, c("OriginalTaxonNo", "ResolvedTaxonNo")], 1, list), function(x) unlist(x)[!is.na(unlist(x))][1]))))
      
    }
    
  # If using taxon numbers:
  } else {
    
    # If a single taxon number then get all children of that taxon and set as new taxon numbers:
    if(length(taxon_nos) == 1) taxon_nos <- gsub("txn:|var:", "", unname(unlist(lapply(apply(PaleobiologyDBDescendantFinder(taxon_nos = taxon_nos, original = original, interval = interval, extant = extant, validonly = validonly, returnrank = returnrank, breaker = breaker)[, c("OriginalTaxonNo", "ResolvedTaxonNo")], 1, list), function(x) unlist(x)[!is.na(unlist(x))][1]))))
    
  }
  
  # Build taxon matrix from what will be tips of tree using taxon numbers:
  TipMatrix <- TaxonMatrix <- metatree::PaleobiologyDBTaxaQuerier(taxon_nos = taxon_nos, original = original, interval = interval, extant = extant, stopfororphans = stopfororphans, breaker = breaker)
  
  # Can now set taxon names safely:
  taxon_names <- TipMatrix[, "TaxonName"]
  
  # Build a parent pool to sample from:
  ParentPool <- sort(unique(gsub("txn:", "", TipMatrix[, "ParentTaxonNo"])))
  
  # While there are still parents to find:
  while(length(ParentPool) > 0) {
    
    # Get parents as new taxa:
    NewTaxa <- metatree::PaleobiologyDBTaxaQuerier(taxon_nos = ParentPool, breaker = breaker)
    
    # Add new taxa to matrix:
    TaxonMatrix <- rbind(TaxonMatrix, NewTaxa)
    
    # Update parent pool from new taxa:
    ParentPool <- sort(unique(gsub("txn:", "", NewTaxa[, "ParentTaxonNo"])))
    
  }
  
  # If obsolete variants are found then replace these with current versions:
  if(any(TaxonMatrix[, "TaxonValidity"] == "obsolete variant of")) TaxonMatrix[which(TaxonMatrix[, "TaxonValidity"] == "obsolete variant of"), c("ResolvedTaxonNo", "TaxonName")] <- TaxonMatrix[which(TaxonMatrix[, "TaxonValidity"] == "obsolete variant of"), c("AcceptedNumber", "AcceptedName")]
  
  # Clean out txn:/var: parts of taxon numbers:
  TaxonMatrix[, c("OriginalTaxonNo", "ResolvedTaxonNo", "ParentTaxonNo")] <- gsub("txn:|var:", "", TaxonMatrix[, c("OriginalTaxonNo", "ResolvedTaxonNo", "ParentTaxonNo")])
  TipMatrix[, c("OriginalTaxonNo", "ResolvedTaxonNo", "ParentTaxonNo")] <- gsub("txn:|var:", "", TipMatrix[, c("OriginalTaxonNo", "ResolvedTaxonNo", "ParentTaxonNo")])
  
  # Build unique parent-child matrix of names:
  ParentChildMatrix <- do.call(rbind, strsplit(unique(paste(TaxonMatrix[, "TaxonName"], TaxonMatrix[match(TaxonMatrix[, "ParentTaxonNo"], TaxonMatrix[, "ResolvedTaxonNo"]), "TaxonName"], sep = "%%%%")), split = "%%%%"))
  
  # Remove NA lines:
  ParentChildMatrix <- ParentChildMatrix[-which(ParentChildMatrix[, 2] == "NA"), ]
  
  # Subfunction to get lineage for a single taxon from a parent child list:
  GetFullLineage <- function(BaseTaxon, ParentChildMatrix) {
    
    # Populate lineage with base taxon:
    Lineage <- BaseTaxon
    
    # Keep adding parents until hit life:
    while(Lineage[length(Lineage)] != "Life") Lineage <- c(Lineage, ParentChildMatrix[ParentChildMatrix[, 1] == Lineage[length(Lineage)], 2])
    
    # Return full lineage:
    return(Lineage)
    
  }
  
  # Get full lineages for all input taxa:
  AllLineages <- lapply(as.list(taxon_names), GetFullLineage, ParentChildMatrix = ParentChildMatrix)
  
  # Identify target clade (first taxon in intersect across all taxa):
  TargetClade <- Reduce(intersect, AllLineages)[1]
  
  # Prune redundant parts of lineages (i.e., higher relationships above):
  AllLineages <- lapply(AllLineages, function(x) x[1:which(x == TargetClade)])
  
  # Create RLE of clade names without target clade:
  RLESansTargetClade <- rle(sort(unlist(lapply(AllLineages, function(x) x[1:(length(x) - 1)]))))
  
  # Build all potential clades (may be empty):
  AllPotentialClades <- RLESansTargetClade$values[RLESansTargetClade$lengths > 1]
  
  # If there are potential clades:
  if(length(AllPotentialClades) > 0) {
    
    #
    CladeComponents <- lapply(as.list(AllPotentialClades), function(x) unlist(lapply(AllLineages, function(y) if(any(y == x)) y[1])))
    
    #
    names(CladeComponents) <- AllPotentialClades
    
    #
    CollapsedNames <- unlist(lapply(CladeComponents, paste, collapse = "%%%%"))
    
    #
    CladeComponents <- lapply(as.list(unique(CollapsedNames)), function(x) {Name <- paste(names(which(CollapsedNames == x)), collapse = "_et_"); Clade <- CladeComponents[which(CollapsedNames == x)[1]]; names(Clade) <- Name; Clade})
    
    # Reorder from smallest clade to largest:
    CladeComponents <- CladeComponents[order(unlist(lapply(CladeComponents, lapply, length)))]
    
    # Remove spaces and replace with underscores:
    CladeComponents <- lapply(CladeComponents, lapply, gsub, pattern = " ", replacement = "_")
    
    # Remove spaces and replace with underscores in taxon_names too:
    taxon_names <- gsub(" ", "_", taxon_names)
    
    # For each clade component in turn (from least to most inclusive):
    for(i in 1:length(CladeComponents)) {
      
      # Find the taxa in current clade:
      TaxaInCurrentCladeComponent <- CladeComponents[[i]][[1]]
      
      # Build new Newick component out of taxa in clade and name of clade:
      NewNewickComponent <- paste("(", paste(TaxaInCurrentCladeComponent, collapse = ","), ")", names(CladeComponents[[i]]), sep = "")
      
      # Update taxon names with newly collapsed clade:
      taxon_names <- c(setdiff(taxon_names, TaxaInCurrentCladeComponent), NewNewickComponent)
      
      # List any itesm in calde components that ned to be udated to collapsed (clades):
      ListItemsToUpdate <- which(unlist(lapply(CladeComponents[1:length(CladeComponents)], function(x) sum(unlist(lapply(as.list(TaxaInCurrentCladeComponent), '==', x[[1]]))))) == length(TaxaInCurrentCladeComponent))
      
      # If there are items in the list that need updating (i..e, because they include taxa that need to be collapsed intoa clade):
      if(length(ListItemsToUpdate) > 0) {
        
        # For each item in list update individual taxa with single clade:
        for(j in ListItemsToUpdate) CladeComponents[[j]][[1]] <- c(setdiff(CladeComponents[[j]][[1]], TaxaInCurrentCladeComponent), NewNewickComponent)
        
      }
      
    }
    
    # Make output into ape tree:
    Tree <- ape::read.tree(text = paste("(", paste(taxon_names, collapse = ","), ")", TargetClade, ";", sep = ""))
    
  # If there are no clades:
  } else {
    
    # Make star tree as output:
    Tree <- ape::read.tree(text = paste("(", paste(gsub(" ", "_", taxon_names), collapse = ","), ")", TargetClade, ";", sep = ""))
    
  }
  
  # If time-scaling the tree:
  if(TimeScale) {
    
    # CHECK NAMES STILL MATCH TIP MATRIX HERE! AND WHAT TO DO IF THEY DO NOT?

    # Collapse tip matrix to just tree taxa:
    TipMatrix <- TipMatrix[match(Tree$tip.label, gsub(" ", "_", TipMatrix[, "TaxonName"])), ]
    
    # Record occurrence numbers to use for checking for tip data:
    OccurrenceNumbers <- unname(unlist(lapply(apply(TipMatrix[, 1:2], 1, as.list), function(x) unlist(x)[!is.na(unlist(x))][1])))
    
    # Get occurrence data for tips:
    TipOccurrences <- PaleobiologyDBOccurrenceQuerier(taxon_nos = OccurrenceNumbers, original = TRUE, breaker = breaker)
    
  }
  
  # If plot is requested then plot tree with node labels included:
  if(plot.tree) {ape::plot.phylo(Tree); ape::nodelabels(gsub("_", "\n", Tree$node.label))}
  
  # Return output:
  return(Tree)
  
}

#PaleobiologyDBTreeBuilder(taxon_names = c("Mei long", "Yi qi", "Kol ghuva", "Talos sampsoni", "Scansoriopteryx heilmanni", "Borogovia gracilicrus", "Zanabazar junior", "Velociraptor mongoliensis", "Velociraptor osmolskae", "Atrociraptor marshalli"), plot.tree = TRUE)
#PaleobiologyDBTreeBuilder(taxon_names = c("Mei long", "Yi qi", "Zby atlanticus"), plot.tree = TRUE)
#PaleobiologyDBTreeBuilder(taxon_names = c("Tyrannosaurus rex", "Triceratops horridus", "Brontosaurus excelsus"), plot.tree = TRUE)

#PaleobiologyDBTreeBuilder(taxon_names = "Thyreophora", plot.tree = TRUE)
#PaleobiologyDBTreeBuilder(taxon_names = "Ceratopsia", plot.tree = TRUE)
#PaleobiologyDBTreeBuilder(taxon_names = "Titanosauria", plot.tree = TRUE)

#Dromaeosauridae <- PaleobiologyDBTreeBuilder(taxon_nos = "38561")
#Trilobita <- PaleobiologyDBTreeBuilder(taxon_nos = "19100")
#Parareptilia <- PaleobiologyDBTreeBuilder(taxon_nos = "99778")
