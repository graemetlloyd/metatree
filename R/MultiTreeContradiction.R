#' Produces tree contradiction distance matrix
#'
#' Given a set of input trees calculates all pariwise tree contradictions.
#'
#' This is effectively an attempt to do for the paleotree function "treeContradiction" what phytools "multiRF" does for Robinson-Foulds, i.e., operate on all pairwise comparisons instead of a single one.
#'
#' @param trees An object of class 'multi.phylo'.
#' @param rescale Whether or not to rescale contradictions (zero to one).
#'
#' @return A contradiction (distance) matrix.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'

#' @references
#' This contradiction difference measure was introduced in:
#'
#' Bapst, D. W., H. A. Schreiber, and S. J. Carlson. 2018. Combined Analysis of Extant Rhynchonellida
#' (Brachiopoda) using Morphological and Molecular Data. \emph{Systematic Biology} 67(1):32-48. doi: 10.1093/sysbio/syx049

#' @seealso
#' A less-optimized function for obtaining the contradiction difference measure
#' is \code{\link[paleotree]{treeContradiction}}, found in package \code{paleotree}.



#' @examples
#'
#' # Nothing yet
#'
#' @export MultiTreeContradiction
MultiTreeContradiction <- function(trees, rescale = TRUE) {
  
  # Test contradiction subfunction (from Dave Bapst):
  testContradiction <- function(namesA, namesB) {
    
    #
    matchA <- namesA %in% namesB
    
    #
    matchB <- namesB %in% namesA
    
    #
    if(any(matchB)) {
      
      #
      Output <- !(all(matchA) | all(matchB))
      
      #
    } else {
      
      #
      Output <- FALSE
      
    }
    
    #
    return(Output)
    
  }
  
  # Number of contradictions subfunction (from Dave Bapst):
  nContradiction <- function(partA, partB) {
    
    #
    partContra <- sapply(partA, function(x) any(sapply(partB, function(y) testContradiction(x, y))))
    
    #
    Output <- sum(partContra)
    
    #
    return(Output)
    
  }
  
  # Get partitions of tree as list (Liam Revell's function):
  PropPartAsList <- function(pp) lapply(pp, function(x, pp) sort(attr(pp, "labels")[x]), pp = pp)
  
  # Get tree contradiction (Dave Bapst's function):
  TreeContradiction <- function(part1, part2) {
    
    # Get number of contradictions comparing parts 1 and 2:
    nContra1 <- nContradiction(part1, part2)
    
    # Get number of contradictions comparing parts 2 and 1:
    nContra2 <- nContradiction(part2, part1)
    
    # Sum values and store as output:
    Output <- nContra1 + nContra2
    
    # Return output:
    return(Output)
    
  }
  
  # Check input format is multiPhylo:
  if(!inherits(trees, "multiPhylo")) stop("trees are not of class multiPhylo")
  
  # Check all trees have same number of tips:
  if(length(unique(unlist(lapply(trees, Ntip)))) != 1) stop("Trees have different numbers of tips")
  
  # Check all trees have same tip labels:
  if(length(unique(unlist(lapply(lapply(lapply(trees, '[[', "tip.label"), sort), paste, collapse = "%%")))) != 1) stop("Trees do not all have same complement of tip labels")
  
  # Get total number of trees:
  NTrees <- length(trees)
  
  # Create all zero distance matrix (will fill later):
  OutputMatrix <- matrix(data = 0, nrow = NTrees, ncol = NTrees)
  
  # Create all partitions list:
  AllPartitionsList <- lapply(unclass(trees), function(x) PropPartAsList(prop.part(x)))
  
  # For the ith tree:
  for(i in 1:(NTrees - 1)) {
    
    # For the jth tree:
    for(j in (i + 1):NTrees) {
      
      # Get tree contradiction (sensu Bapst) and store:
      OutputMatrix[i, j] <- OutputMatrix[j, i] <- TreeContradiction(part1 = AllPartitionsList[[i]], part2 = AllPartitionsList[[j]])
      
    }
    
  }
  
  # If rescaling output:
  if(rescale) {
    
    # Number of possible nodes that could contradict on one unrooted tree:
    nPossNodes <- Ntip(trees[[1]]) - 2
    
    # Rescale matrix (0-1 range):
    OutputMatrix <- OutputMatrix / (2 * nPossNodes)
    
  }
  
  # Return distance matrix:
  return(OutputMatrix)
  
}
