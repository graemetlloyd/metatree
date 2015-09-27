#' Collapse MRP matrix
#' 
#' Collapses an MRP matrix to unique characters
#' 
#' Collapses an MRP matrix to just unique characters followung the protocol aid out in Pisani et al. (2002).
#' 
#' @param matrix An MRP matrix in \link{ReadMorphNexus} format.
#'
#' @return An MRP matrix in \link{ReadMorphNexus} format where all characters are unique.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @reference
#'
#' Pisani et al. (2002)
#'
#' @examples
#' 
#' # Nothing yet.
#'
#' @export MRPCollapse
MRPCollapse <- function(matrix) {
    
    matrix <- as.matrix(matrix)
    
    taxon.names <- rownames(matrix)
    
    transposed.matrix <- t(matrix)
    
    unq.rows <- vector(mode = "numeric")
    
    for(i in length(transposed.matrix[, 1]):1) unq.rows <- unique(c(unq.rows, paste(transposed.matrix[i, ], collapse = "")))

    unq.rows <- sort(unique(unq.rows), decreasing = TRUE)
    
    if(length(unq.rows) > 0) {
        
        for(i in length(unq.rows):1) {
            
            if(sum(as.numeric(strsplit(unq.rows[i], "")[[1]])) < 2) {
                
                unq.rows <- unq.rows[-i]
                
            }
            
        }
        
        if(length(unq.rows) > 0) {
            
            transposed.matrix <- vector(mode = "character")
            
            for(i in 1:length(unq.rows)) {
                
                transposed.matrix <- rbind(transposed.matrix, strsplit(unq.rows[i], "")[[1]])
                
            }
            
            if(nrow(as.matrix(transposed.matrix)) > 1) {
                
                out <- t(transposed.matrix)
                
                rownames(out) <- taxon.names
                
            } else {
                
                out <- as.matrix(transposed.matrix)
                
                colnames(out) <- taxon.names
                
                out <- t(out)
                
            }
            
        } else {
            
            out <- "There are no informative characters after collapsing"
            
        }
        
    } else {
        
        out <- "There are no informative characters after collapsing"
        
    }
    
    return(out)
    
}
