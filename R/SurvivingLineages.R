#' Finds unique surviving lineages of an extinction
#'
#' @description
#'
#' Finds unique surviving lineages of an extinction given a tree and vector of binary (extinction or survival) values.
#' 
#' @param Tree A single tree in ape's phylo format.
#' @param SurvivorVictim A numeric vector of victims (0) and survivors (1) with names matching to Tree$tip.label.
#'
#' @details
#'
#' Given a tree that crosses some extinction boundary, with both extinct and surviving species as tips, we might be interested in knowing how many individual lineages survived. This function finds unique lineages (sets of tips) that survive the boundary.
#'
#' @return
#'
#' A list of surviving lineages (may be clades - vectors of length greater than one) or individual tips (vectors of length one).
#'
#' @author
#'
#' Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Make sure random examples work:
#' set.seed(17)
#'
#' # Generate a random 10 tip tree:
#' Tree <- ape::rtree(10)
#'
#' # Remove branch lengths (not used here):
#' Tree$edge.length <- NULL
#'
#' # Create random vector of survivors (!) and victims (0):
#' SurvivorVictim <- sample(c(0, 1), size = 10, replace = TRUE)
#'
#' # Assign tip names to survivors and victims:
#' names(SurvivorVictim) <- Tree$tip.label
#'
#' # Run surviving lineages function on data:
#' SurvivingLineages(Tree, SurvivorVictim)
#'
#' # Visually confirm results:
#' plot(Tree)
#' ape::tiplabels(SurvivorVictim)
#'
#' # Generate multiple trees with same tip names:
#' Trees <- ape::rmtree(10, 10)
#'
#' # Calculate just total number of survivors on each tree:
#' unlist(lapply(Trees, function(x) length(SurvivingLineages(x, SurvivorVictim = SurvivorVictim))))
#'
#' @export SurvivingLineages
SurvivingLineages <- function(Tree, SurvivorVictim) {
  
  # MORE DATA CHECKS AT TOP
  
  # Check there are both survivors and vuictims and stop and warn user if not:
  if(length(unique(SurvivorVictim)) != 2) stop("SurvivorVictim must consist of both 0 (victim) and 1 (survivor) values.")
  
  # Make vector of internal nodes:
  InternalNodesVector <- (ape::Ntip(Tree) + 1):(ape::Ntip(Tree) + ape::Nnode(Tree))
  
  # Create list of all survivor clades:
  AllSurvivorClades <- lapply(as.list(InternalNodesVector), function(x) {Descendants <- SurvivorVictim[Tree$tip.label[strap::FindDescendants(n = x, tree = Tree)]]; if(all(Descendants == 1)) return(names(Descendants))})
  
  # Collapse to just the surviving clades (may end up an empty list if there are none):
  AllSurvivorClades <- AllSurvivorClades[which(unlist(lapply(AllSurvivorClades, length)) > 0)]
  
  # Create emoty vector of nested clades:
  NestedClades <- vector(mode = "numeric")
  
  # If there are all survivor clades
  if(length(AllSurvivorClades) > 0) {
    
    # If some of those clades are nested inside each other:
    if(any(duplicated(unlist(AllSurvivorClades)))) {
      
      # For each surviving clade:
      for(i in 1:length(AllSurvivorClades)) {
        
        # Find any overlapping clades:
        OverlappingClades <- base::setdiff(which(unlist(lapply(as.list(AllSurvivorClades), function(x) base::any(base::duplicated(base::c(x, AllSurvivorClades[[i]])))))), i)
        
        # If there is an overlapping clade:
        if(base::length(OverlappingClades) > 0) {
          
          # If ith clade is nested inside at least one other then add it to nested clades variable for deleting later:
          if(base::any(lapply(AllSurvivorClades[OverlappingClades], function(x) length(setdiff(x, AllSurvivorClades[[i]]))) > 0)) NestedClades <- c(NestedClades, i)
          
        }
        
      }
      
    }
    
  }
  
  # If there are nested clades to remove then remove them:
  if(length(NestedClades) > 0) AllSurvivorClades <- AllSurvivorClades[-NestedClades]
  
  # Add sinlge lineages to all survivor clades (if any):
  AllSurvivorClades <- c(AllSurvivorClades, as.list(setdiff(names(which(SurvivorVictim == 1)), unlist(AllSurvivorClades))))
  
  # Return list of all surviving lineages:
  return(AllSurvivorClades)
  
}
