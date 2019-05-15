#' Converts tree(s) to an MRP matrix
#' 
#' Given an input tree or trees creates an MRP matrix.
#' 
#' Matrix Representtaion with parsimony (Baum 1992; Ragan 1992).
#' 
#' @param Trees An object of class 'phylo' or 'multi.phylo'.
#'
#' @return An MRP matrix in \link{ReadMorphNexus} format.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Baum (1992)
#'
#' Ragan (1992)
#'
#' @examples
#' 
#' # Nothing yet
#'
#' @export Tree2MRP
Tree2MRP <- function(Trees) {
    
    # Allow Purvis format?
    # Check tip names match across all trees!
    
    # If a single tree convert to a list so rest of code works:
    if(class(Trees) == "phylo") Trees <- as.list(Trees)
    
    # Get all unique MRP vectors across all trees:
    MRPVectors <- unique(unlist(lapply(Trees, function(x) lapply(as.list((ape::Ntip(x) + 1):(ape::Ntip(x) + ape::Nnode(x))), function(y) {InClade <- x$tip.label[strap::FindDescendants(y, x)]; MRPVector <- rep(0, ape::Ntip(x)); names(MRPVector) <- x$tip.label; MRPVector[InClade] <- 1; paste(MRPVector[sort(names(MRPVector))], collapse = "")}))))
    
    # Build into MRP matrix:
    MRPMatrix <- do.call(cbind, lapply(MRPVectors, function(x) strsplit(x, split = "")[[1]]))
    
    # All zero outgroup at top:
    MRPMatrix <- rbind(rep(0, ncol(MRPMatrix)), MRPMatrix)
    
    # Add taxa as rwonames:
    rownames(MRPMatrix) <- c("allzero", sort(Trees[[1]]$tip.label))
    
    # Collapse to just variable characters:
    MRPMatrix <- MRPMatrix[, unlist(lapply(apply(MRPMatrix, 2, list), function(x) length(unique(unlist(x))))) == 2]
    
    # Format as Claddis matrix:
    MRPMatrix <- Claddis::MakeMorphMatrix(MRPMatrix)
    
    # Return matrix as output:
    return(MRPMatrix)
    
}
