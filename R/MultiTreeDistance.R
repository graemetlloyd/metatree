#' Pairwise distance matrix for a set of trees
#'
#' Given a set of input trees, calculates all pair-wise tree distances.
#'
#' This is effectively an attempt to both do a multi-tree version of the paleotree function
#' \code{\link[paleotree]{treeContradiction}} and a faster version of phytools' \code{\link[phytools]{multiRF}} function when
#' dealing with very large number of trees or tips (hundreds of trees or tips).
#' The function works by converting the partitions into matrices and operating on those rather than the tree objects.
#' This can substantially reduce the memory required and hence sped up the calculations, but only when data sets are large.
#' 

#' @param trees An object of class 'multi.phylo'.
#' @param distance The type of disatnce to use, must be one of "RF" (RObinson-Foulds) or "contradiction" (tree contradiction).
#' @param scale Whether or not to scale the distance (zero to one).
#'

#' @return A pairwise tree distance matrix.
#'

#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'

#' @references
#' This contradiction difference measure was introduced in:
#'
#' Bapst, D. W., H. A. Schreiber, and S. J. Carlson. 2018. Combined Analysis of Extant Rhynchonellida
#' (Brachiopoda) using Morphological and Molecular Data. \emph{Systematic Biology} 67(1):32-48. doi: 10.1093/sysbio/syx049

# Bapst paper (see paleotree package and treeContradiction function)

#' @seealso
#' A less-optimized function for obtaining the contradiction difference measure
#' is \code{\link[paleotree]{treeContradiction}}, found in package \code{paleotree}.

#'
#' @examples
#'
#' # Nothing yet
#'
#' @export MultiTreeDistance
MultiTreeDistance <- function(trees, distance = "contradiction", scale = TRUE) {
  
  # Unrooted trees are problematic as they can never be scaled to one.
  # I.e., they have a polytomy and hence without maximising bifurcations they can never realise the maximum distance.
  # Further this cannot be correceted for easily in the rescale set as this simply deals with the issue of the basal partition (everything) which all trees will have.
  
  # Check input format is multiPhylo:
  if(!inherits(trees, "multiPhylo")) stop("trees are not of class multiPhylo")
  
  # Check all trees have same number of tips:
  if(length(unique(unlist(lapply(trees, Ntip)))) != 1) stop("Trees have different numbers of tips")
  
  # Check all trees have same tip labels:
  if(length(unique(unlist(lapply(lapply(lapply(trees, '[[', "tip.label"), sort), paste, collapse = "%%")))) != 1) stop("Trees do not all have same complement of tip labels")
  
  # Check distance is a valid type:
  if(distance != "RF" && distance != "contradiction") stop("distance must be one of \"RF\" or \"contradiction\"")
  
  # Little function to replace taxon names with numbers in order:
  TipLabelSimplifier <- function(x) {
    
    # Replace tip labels with numbers:
    x$tip.label <- as.character(match(sort(trees[[1]]$tip.label), x$tip.label))
    
    # Return relabelled tree:
    return(x)
    
  }
  
  # Simplify taxon names to just numbers:
  trees <- lapply(trees, TipLabelSimplifier)
  
  # Restore multiPhylo class:
  class(trees) <- "multiPhylo"
  
  # Get total number of trees:
  NTrees <- length(trees)
  
  # Get total number of tips:
  NTips <- Ntip(trees[[1]])
  
  # Create all zero distance matrix (will fill later):
  OutputMatrix <- matrix(data = 0, nrow = NTrees, ncol = NTrees)
  
  # Get partitions of tree as list (stolen from Liam Revell's multiRF function):
  PropPartAsList <- function(pp) lapply(pp, function(x, pp) sort(attr(pp, "labels")[x]), pp = pp)
  
  # Create all partitions list:
  AllPartitionsList <- lapply(unclass(trees), function(x) PropPartAsList(prop.part(x)))
  
  # Sort so can match later without issues of a partition 123 being seens as different to 132:
  AllPartitionsList <- lapply(AllPartitionsList, lapply, sort)
  
  # Collapse paritions to single string with separator:
  AllPartitionsList <- lapply(AllPartitionsList, lapply, paste, collapse = "%%")
  
  # Comvert sublists to vectors:
  AllPartitionsList <- lapply(AllPartitionsList, unlist)
  
  # Get all unique partitions across all trees:
  AllUniquePartitions <- unique(unlist(AllPartitionsList))
  
  # If using Bapst's contradiction metric:
  if(distance == "contradiction") {
    
    # Get paritions to check (drops root which will not contradict anything):
    PartitionsToCheck <- setdiff(AllUniquePartitions, paste(sort(trees[[1]]$tip.label), collapse = "%%"))
    
    # Count of numbe rof paritions to check:
    NPartitionsToCheck <- length(PartitionsToCheck)
    
    # Create empty partiton by partition contradcitions matrix (1 indicating a contradiction between the ith and jth partitions):
    ContradictionMatrix <- matrix(0, nrow = NPartitionsToCheck, ncol = NPartitionsToCheck, dimnames = list(PartitionsToCheck, PartitionsToCheck))
    
    # Create list object of partitons to check:
    PartitionsToCheckList <- strsplit(PartitionsToCheck, split = "%%")
    
    # Add anems of partitons to list:
    names(PartitionsToCheckList) <- PartitionsToCheck
    
    # Order partitokn list from shortest to longest (speeds up later on as means largets number of clades are looked at with smallest number of taxa):
    PartitionsToCheckList <- PartitionsToCheckList[order(unlist(lapply(PartitionsToCheckList, length)))]
    
    # For each partioon to check (from samllest to largest):
    for(i in 1:(length(PartitionsToCheckList) - 1)) {
      
      # Get length (number of taxa) in current partition:
      CurrentPartitionLength <- length(PartitionsToCheckList[[i]])
      
      # Get length of interscetion with all remaning partitions:
      IntersectLengths <- unlist(lapply(lapply(PartitionsToCheckList[(i + 1):length(PartitionsToCheckList)], intersect, PartitionsToCheckList[[i]]), length))
      
      # Contradictory clades (i.e., will have an intersection of greater than zero and less than N, where N is all taxa in current partition):
      Contradictions <- names(IntersectLengths[intersect(which(IntersectLengths < CurrentPartitionLength), which(IntersectLengths > 0))])
      
      # If contradictins were found then store them in a matrix:
      if(length(Contradictions) > 0) ContradictionMatrix[names(PartitionsToCheckList)[i], Contradictions] <- ContradictionMatrix[Contradictions, names(PartitionsToCheckList)[i]] <- 1
      
    }
    
  }
  
  # Create empty partitions matrix:
  PartitionsMatrix <- matrix(0, ncol = length(AllUniquePartitions), nrow = NTrees, dimnames = list(c(), AllUniquePartitions))
  
  # Populate partition matrix for each tree:
  for(i in 1:NTrees) PartitionsMatrix[i, AllPartitionsList[[i]]] <- 1
  
  # Remove uselss root partition as will always be present and always appear in all trees:
  PartitionsMatrix <- PartitionsMatrix[, -which(colnames(PartitionsMatrix) == paste(sort(trees[[1]]$tip.label), collapse = "%%"))]
  
  # For every tree but the last:
  for(i in 1:(NTrees - 1)) {
    
    # Bind ith row to ith + 1 through N rows:
    BoundRows <- lapply(lapply(lapply(lapply(as.list(apply(PartitionsMatrix[(i + 1):NTrees, , drop = FALSE], 1, paste, collapse = "%%")), strsplit, split = "%%"), unlist), as.numeric), cbind, PartitionsMatrix[i, ])
    
    # If using RF distances get distances and store:
    if(distance == "RF") OutputMatrix[i, (i + 1):NTrees] <- OutputMatrix[(i + 1):NTrees, i] <- unlist(lapply(lapply(lapply(BoundRows, apply, 1, sum), '==', 1), sum))
    
    # If using tree contradiction:
    if(distance == "contradiction") {
      
      # Get contradictions sub matrix:
      ContradictionSubmatrix <- lapply(lapply(lapply(BoundRows, function(x) x[apply(x == 1, 1, sum) == 1, , drop = FALSE]), function(x) list(rownames(x)[x[, 1] == 1], rownames(x)[x[, 2] == 1])), function(x) ContradictionMatrix[x[[1]], x[[2]], drop = FALSE])
      
      # Sum rows and columns that contradict to get score for pair of trees and store:
      OutputMatrix[i, (i + 1):NTrees] <- OutputMatrix[(i + 1):NTrees, i] <- unlist(lapply(lapply(ContradictionSubmatrix, apply, 1, max), sum)) + unlist(lapply(lapply(ContradictionSubmatrix, apply, 2, max), sum))
      
    }
    
  }
  
  # If scaling trees divide through by maximum possible value:
  if(scale) OutputMatrix <- OutputMatrix / ((NTips - 2) * 2)
  
  # Return distacne matrix:
  return(OutputMatrix)
  
}
