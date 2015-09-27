#' Converts TNT MRP output to NEXUS format
#' 
#' Converts TNT formatted MRP matrix to
#' 
#' TNT (Goloboff et al. 2008) \\#NEXUS (Maddison ????)
#' 
#' @param filename The TNT input file name.
#'
#' @return A cladistic matrix in \link{ReadMorphNexus} format.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#' 
#' # Nothing yet
#'
#' @export TNTMRP2NEXUSMRP
TNTMRP2NEXUSMRP <- function(filename) {
    
    # NEEDS A BETTER NAME, NOT REALLY CONVERTING TO NEXUS
    
    # Read in raw file:
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE) # Read in NEXUS file
    start <- grep("ROOT", X)+1
    end <- grep(";", X[start:length(X)])[1]+start-2
    MATRIX <- X[start:end]
    while(length(grep("  ", MATRIX)) > 0) MATRIX <- gsub("  ", " ", MATRIX)
    characters <- names <- vector(mode="character")
    for(i in 1:length(MATRIX)) {
        names[i] <- strsplit(MATRIX[i], " ")[[1]][1]
        characters[i] <- strsplit(MATRIX[i], " ")[[1]][2]
    }
    char.block <- matrix(nrow=length(names), ncol=nchar(characters[1]))
    for(i in 1:length(MATRIX)) char.block[i, ] <- strsplit(characters[i], "")[[1]]
    char.block <- t(char.block)
    characters <- vector(mode="character")
    for(i in 1:length(char.block[, 1])) characters[i] <- paste(char.block[i, ], collapse = "")
    characters <- sort(characters)
    characters <- rle(characters)
    weights <- characters$lengths
    characters <- characters$values
    char.block <- char.block[1:length(characters), ]
    for(i in 1:length(characters)) char.block[i, ] <- strsplit(characters[i], "")[[1]]
    char.block <- t(char.block)
    rownames(char.block) <- names
    # Make into clad.matrix format:
    header <- ""
    ordering <- rep("unord", length(char.block[1, ]))
    max.vals <- rep(1, length(char.block[1, ]))
    min.vals <- rep(0, length(char.block[1, ]))
    result <- list(header, char.block, ordering, weights, max.vals, min.vals)
    names(result) <- c("header", "matrix", "ordering", "weights", "max.vals", "min.vals")
    return(result)
}
