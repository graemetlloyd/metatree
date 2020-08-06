#' Collapses an MRP matrix with duplicate taxa
#'
#' @description
#'
#' Collapses an MRP matrix with duplicate taxa until all taxa are unique.
#'
#' @param MRPMatrix An MRP matrix in the format imported by \code{Claddis::read_nexus_matrix}.
#'
#' @details
#'
#' Duplicate taxon names are problematic for downstream analysis. Generally speaking users should ensure all OTU names are unique. However, in some cases it might make more sense to collapse the matrix until all taxa are unique. Typically this will involve reducing the row count (unique taxa means fewer taxa), but increasing the character count (duplicate taxa means more unique codings).
#'
#' This function automates this process and is used internally by the \link{Metatree} function, but is also made available for individual use here.
#'
#' @author
#'
#' Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @keywords
#
#' MRP,NEXUS
#'
#' @seealso
#'
#' For another means by which MRP matrices require collapsing see \link{PisaniMRPPrune}.
#'
#' @examples
#'
#' # Build an example matrix:
#' ExampleMRP <- Claddis::build_cladistic_matrix(
#'   character_taxon_matrix = matrix(as.character(c(0, 0,
#'   0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)),
#'   nrow = 6, dimnames = list(LETTERS[1:6], c()), byrow = TRUE)
#' )
#'
#' # Show taxon names are currently unique by inspection:
#' rownames(ExampleMRP$matrix_1$matrix)
#'
#' # Now rename taxon F as taxon E:
#' rownames(ExampleMRP$matrix_1$matrix)[6] <- "E"
#'
#' # Show taxon E is duplicated by inspection:
#' rownames(ExampleMRP$matrix_1$matrix)
#'
#' # Now collapse to just unique taxa using function:
#' ExampleMRP <- CollapseDuplicateTaxonMRP(ExampleMRP)
#'
#' # Confirm that matrix has now been collapsed to five unique taxa
#' # (and four unique characters) by inspection:
#' ExampleMRP$matrix_1$matrix
#'
#' @export CollapseDuplicateTaxonMRP
CollapseDuplicateTaxonMRP <- function(MRPMatrix) {
  
  # Check there are no step matrices:
  if(any(!is.null(MRPMatrix$Topper$StepMatrices))) stop("Function not meant for use with step matrices.")
  
  # Check is only a single block and stop and warn user if not:
  if(length(MRPMatrix) > 2) stop("MRP matrices should consist of a single block.")
  
  # Check is only zeroes and ones and stop and warn user if not:
  if(length(setdiff(unique(as.vector(MRPMatrix$matrix_1$matrix)), c("0", "1"))) > 0) stop("MRP matrix must consist on only zeroes and ones.")
  
  # Check for zero weight characters and stop and warn user if found:
  if(any(MRPMatrix$matrix_1$weights == 0)) stop("Function not intended for zero weight characters.")
  
  # Check taxa are actually duplicated:
  if(!any(duplicated(rownames(MRPMatrix$matrix_1$matrix)))) stop("No taxa are actually duplicated.")
  
  # Store unique taxon names in order (to reorder matrix at end):
  UniqueTaxonNames <- unique(rownames(MRPMatrix$matrix_1$matrix))
  
  # Build a vector of unique duplicated names:
  DuplicatedNames <- sort(unique(rownames(MRPMatrix$matrix_1$matrix)[duplicated(rownames(MRPMatrix$matrix_1$matrix))]))
  
  # Subfunction to collapse a single duplicated taxon:
  CollapseDuplicateTaxonMRPMatrix <- function(MRPMatrix, DuplicatedTaxon) {
    
    # Find rows corresponding to duplicated taxa:
    DuplicateRows <- which(rownames(MRPMatrix$matrix_1$matrix) == DuplicatedTaxon)
    
    # Isolate block corresponding to duplicated taxa:
    DuplicatedBlock <- MRPMatrix$matrix_1$matrix[DuplicateRows, , drop = FALSE]
    
    # Get unique characters to use:
    CharactersToUse <- unlist(mapply(function(x, y) rep(x, y), x = as.list(1:ncol(DuplicatedBlock)), y = as.list(apply(DuplicatedBlock, 2, function(x) length(unique(x))))))
    
    # Build new block (without duplicated taxon):
    NewBlock <- MRPMatrix$matrix_1$matrix[-DuplicateRows, CharactersToUse, drop = FALSE]
    
    # Add duplicated taxon to new block:
    NewBlock <- rbind(NewBlock, unlist(apply(DuplicatedBlock, 2, function(x) as.list(unique(x)))))
    
    # Add duplicated taxon as row name:
    rownames(NewBlock)[nrow(NewBlock)] <- DuplicatedTaxon
    
    # Replace current block with new block in MRP matrix:
    MRPMatrix$matrix_1$matrix <- NewBlock
    
    # Update ordering in MRP matrix:
    MRPMatrix$matrix_1$ordering <- MRPMatrix$matrix_1$ordering[CharactersToUse]
    
    # Update weights in MRP matrix:
    MRPMatrix$matrix_1$weights <- MRPMatrix$matrix_1$weights[CharactersToUse]
    
    # Update minimum values in MRP matrix:
    MRPMatrix$matrix_1$minimum_values <- MRPMatrix$matrix_1$minimum_values[CharactersToUse]
    
    # Update maximum values in MRP matrix:
    MRPMatrix$matrix_1$maximum_values <- MRPMatrix$matrix_1$maximum_values[CharactersToUse]
    
    # Return modified MRP matrix:
    return(MRPMatrix)
   
  }
  
  # For each duplicated taxon collapse matrix appropriately (must be done sequentially):
  for(i in DuplicatedNames) MRPMatrix <- CollapseDuplicateTaxonMRPMatrix(MRPMatrix = MRPMatrix, DuplicatedTaxon = i)
  
  # Resort matrix by original taxon ordering:
  MRPMatrix$matrix_1$matrix <- MRPMatrix$matrix_1$matrix[UniqueTaxonNames, , drop = FALSE]
  
  # Return MRP matrix:
  return(MRPMatrix)
  
}
