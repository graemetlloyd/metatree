#' Collapses an MRP matrix with duplicate taxa
#'
#' @description
#'
#' Collapses an MRP matrix with duplicate taxa until all taxa are unique.
#'
#' @param MRPMatrix An MRP matrix in the format improted by \code{Claddis::ReadMorphNexus}.
#'
#' @details
#'
#' Duplicate taxon names are problematic for downstream anslysis. Generally speaking users should ensure all OTU names are unique. However, in some cases it might make more sense to collapse the matrix until all taxa are unique. Typically this will involve reducing the row count (unique taxa means fewer taxa), but increasing the character count (duplicate taxaon means more unique codings)..
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @keywords NEXUS
#'
#' @examples
#'
#' # Nothing yet.
#'
#' @export CollapseDuplicateTaxonMRP
CollapseDuplicateTaxonMRP <- function(MRPMatrix) {
  
  # Check there are no step matrices:
  if(any(!is.null(MRPMatrix$Topper$StepMatrices))) stop("Function not meant for use with step matrices.")
  
  # Check is only a single block and stop and warn user if not:
  if(length(MRPMatrix) > 2) stop("MRP matrices should consist of a single block.")
  
  # Check is only zeroes and ones and stop and warn user if not:
  if(length(setdiff(unique(as.vector(MRPMatrix$Matrix_1$Matrix)), c("0", "1"))) > 0) stop("MRP matrix must consist on only zeroes and ones.")
  
  # Check for zero weight characters and stop and warn user if found:
  if(any(MRPMatrix$Matrix_1$Weights == 0)) stop("Function not intended for zero weight characters.")
  
  # Check taxa are actually duplicated:
  if(!any(duplicated(rownames(MRPMatrix$Matrix_1$Matrix)))) stop("No taxa are actually duplicated.")
  
  # Store unique taxon names in order (to reorder matrix at end):
  UniqueTaxonNames <- unique(rownames(MRPMatrix$Matrix_1$Matrix))
  
  # Build a vector of unique duplicated names:
  DuplicatedNames <- sort(unique(rownames(MRPMatrix$Matrix_1$Matrix)[duplicated(rownames(MRPMatrix$Matrix_1$Matrix))]))
  
  # Subfunction to collapse a single duplicated taxon:
  CollapseDuplicateTaxonMRPMatrix <- function(MRPMatrix, DuplicatedTaxon) {
    
    # Find rows corresponding to duplicated taxa:
    DuplicateRows <- which(rownames(MRPMatrix$Matrix_1$Matrix) == DuplicatedTaxon)
    
    # Isolate block corresponding to duplicated taxa:
    DuplicatedBlock <- MRPMatrix$Matrix_1$Matrix[DuplicateRows, ]
    
    # Get unique characters to use:
    CharactersToUse <- unlist(mapply(function(x, y) rep(x, y), x = as.list(1:ncol(DuplicatedBlock)), y = as.list(apply(DuplicatedBlock, 2, function(x) length(unique(x))))))
    
    # Build new block (without duplicated taxon):
    NewBlock <- MRPMatrix$Matrix_1$Matrix[-DuplicateRows, CharactersToUse]
    
    # Add duplicatd taxon to new block:
    NewBlock <- rbind(NewBlock, unlist(apply(DuplicatedBlock, 2, function(x) as.list(unique(x)))))
    
    # Add duplicated taxon as row name:
    rownames(NewBlock)[nrow(NewBlock)] <- DuplicatedTaxon
    
    # Replace current block with new block in MRP matrix:
    MRPMatrix$Matrix_1$Matrix <- NewBlock
    
    # Update ordering in MRP matrix:
    MRPMatrix$Matrix_1$Ordering <- MRPMatrix$Matrix_1$Ordering[CharactersToUse]
    
    # Update weights in MRP matrix:
    MRPMatrix$Matrix_1$Weights <- MRPMatrix$Matrix_1$Weights[CharactersToUse]
    
    # Update minimum values in MRP matrix:
    MRPMatrix$Matrix_1$MinVals <- MRPMatrix$Matrix_1$MinVals[CharactersToUse]
    
    # Update maximum values in MRP matrix:
    MRPMatrix$Matrix_1$MaxVals <- MRPMatrix$Matrix_1$MaxVals[CharactersToUse]
    
    # Return modified MRP matrix:
    return(MRPMatrix)
   
  }
  
  # For each duplicated taxon collapse matrix appropriately (must be done sequentially):
  for(i in DuplicatedNames) MRPMatrix <- CollapseDuplicateTaxonMRPMatrix(MRPMatrix = MRPMatrix, DuplicatedTaxon = i)
  
  # Resort matrix by original taxon ordering:
  MRPMatrix$Matrix_1$Matrix <- MRPMatrix$Matrix_1$Matrix[UniqueTaxonNames, ]
  
  # Return MRP matrix:
  return(MRPMatrix)
  
}
