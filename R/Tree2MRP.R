#' Converts tree(s) to an MRP matrix
#'
#' @description
#'
#' Given an input tree or trees creates an MRP matrix.
#'
#' @param Trees An object of class 'phylo' or 'multi.phylo'.
#' @param AddAllZero Whether (TRUE; default) to add an allzero outgroup taxon or not (FALSE).
#'
#' @details
#'
#' The \link{Metatree} function requires samples of tree(s) as input, but these must already be represented using Matrix Representtaion with Parsimony (MRP; Baum 1992; Ragan 1992). This encoding works by taking each bipartition (clade, or internal node, excluding the root) and coding taxa present \emph{inside} the clade as 1 ("derived" state) and those outside the clade as 0 "primitive" state). Thus the tree:
#'
#' \preformatted{ /-----A
#'  |
#' R| /---B
#'  \-|
#'   X| /-C
#'    \-|
#'     Y\-D}
#'
#' Would be encoded as the matrix:
#'
#' \preformatted{  XY
#' A 00
#' B 10
#' C 11
#' D 11}
#'
#' Or, if using \code{AddAllZero = TRUE} the root (R) can also be encoded:
#'
#' \preformatted{        RXY
#' allzero 000
#' A       100
#' B       110
#' C       111
#' D       111}
#'
#' Note that the function will remove any duplicate characters automatically, meaning if nodes reappear across a sample of trees these will be collapsed to a single MRP character (see example below), as suggested by Lloyd et al. (2016).
#'
#' Currently an option to perform Purvis MRP (Purvis 1995) is not available.
#'
#' @return
#'
#' An MRP matrix in \link{Claddis::ReadMorphNexus} format.
#'
#' @author
#'
#' Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Baum, B. R., 1992. Combining trees as a way of combining data sets for phylogenetic inference, and the desirability of combining gene trees. \emph{Taxon}, \bold{41}, 3-10.
#'
#' Lloyd, G. T., Bapst, D. W., Friedman, M. and Davis, K. E., 2016. Probabilistic divergence time estimation without branch lengths: dating the origins of dinosaurs, avian flight, and crown birds. \emph{Biology Letters}, \bold{12}, 20160609.
#'
#' Purvis, A., 1995. A modification to Baum and Ragan's method for combining phylogenetic trees. \emph{Systematic Biology}, \bold{44}, 251-255.
#'
#' Ragan, M., 1992. Phylogenetic inference based on matrix representation of trees. \emph{Molecular Phylogenetics and Evolution}, \bold{1}, 113-126.
#'
#' @examples
#' 
#' # Generate a set of two example trees:
#' ExampleTrees <- ape::read.tree(text = c("(A,(B,(C,D)));", "(A,(C,(B,D)));"))
#'
#' # Convert to MRP and show just matrix:
#' Tree2MRP(ExampleTrees)$Matrix_1$Matrix
#'
#' # To confirm this collapses duplicate nodes show individual tree results:
#' Tree2MRP(ExampleTrees[[1]])$Matrix_1$Matrix
#' Tree2MRP(ExampleTrees[[2]])$Matrix_1$Matrix
#'
#' @export Tree2MRP
Tree2MRP <- function(Trees, AddAllZero = TRUE) {
    
    # Allow Purvis format?
    # Check tip names match across all trees!
    # Allow weighting by frequency
    
    # If a single tree convert to a list so rest of code works:
    if(class(Trees) == "phylo") Trees <- list(Trees)
    
    # Get all unique MRP vectors across all trees:
    MRPVectors <- unique(unlist(lapply(Trees, function(x) lapply(as.list((ape::Ntip(x) + 1):(ape::Ntip(x) + ape::Nnode(x))), function(y) {InClade <- x$tip.label[strap::FindDescendants(y, x)]; MRPVector <- rep(0, ape::Ntip(x)); names(MRPVector) <- x$tip.label; MRPVector[InClade] <- 1; paste(MRPVector[sort(names(MRPVector))], collapse = "")}))))
    
    # Build into MRP matrix:
    MRPMatrix <- do.call(cbind, lapply(MRPVectors, function(x) strsplit(x, split = "")[[1]]))
    
    # All zero outgroup at top:
    if(AddAllZero) MRPMatrix <- rbind(rep(0, ncol(MRPMatrix)), MRPMatrix)
    
    # Add taxa as rownames:
    if(AddAllZero) rownames(MRPMatrix) <- c("allzero", sort(Trees[[1]]$tip.label))
    
    # Collapse to just variable characters:
    MRPMatrix <- MRPMatrix[, unlist(lapply(apply(MRPMatrix, 2, list), function(x) length(unique(unlist(x))))) == 2]
    
    # Format as Claddis matrix:
    MRPMatrix <- Claddis::MakeMorphMatrix(MRPMatrix)
    
    # Return matrix as output:
    return(MRPMatrix)
    
}
