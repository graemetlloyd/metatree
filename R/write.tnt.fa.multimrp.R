#' Counts the changes in a series of time bins
#' 
#' Given a vector of dates for a series of time bins and another for the times when a character change occurred will return the total number of changes in each bin.
#' 
#' Calculates the total number of evolutionary changes in a series of time bins. This is intended as an internal function for rate calculations, but could be used for other purposes (e.g., counting any point events in a series of time bins).
#' 
#' @param change.times A vector of ages in millions of years at which character changes are hypothesised to have occurred.
#' @param time.bins A vector of ages in millions of years of time bin boundaries in old-to-young order.
#'
#' @return A vector giving the number of changes for each time bin. Names indicate the maximum and minimum (bottom and top) values for each time bin.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#' 
#' # Create a random dataset of 100 changes:
#' change.times <- runif(100, 0, 100)
#' 
#' # Create time bins:
#' time.bins <- seq(100, 0, length.out=11)
#' 
#' # Get N changes for each bin:
#' ChangesInBins(change.times, time.bins)
#' 
#' @export ChangesInBins
write.tnt.fa.multimrp <- function(clad.matrix, file.name, n.iterations=20, n.trees=5000) {
    head <- clad.matrix$header # Separate out header
    mat <- clad.matrix$matrix # Separate out matrix
    ord <- clad.matrix$ordering # Separate out ordering
    wts <- clad.matrix$weights # Separate out weights
    ntaxa <- length(mat[, 1]) # Get number of taxa
    nchars <- length(mat[1, ]) # Get number of characters
    head <- gsub("'", "", head) # Remove any single quotes from header
    states <- sort(unique(as.vector(mat))) # Get states listed in matrix
    if(length(grep("&", states)) > 0) { # If polymorphisms are present
        while(length(grep("&", states)) > 0) { # For each polymorphism type
            temp <- strsplit(states[grep("&", states)[1]], "&")[[1]] # Split into separate components
            states <- states[-grep("&", states)[1]] # Remove type from list
            states <- unique(c(states, temp)) # Add components (if novel)
        }
    }
    if(length(states) > 8) { # If there are more than 8 states
        if(length(states) > 16) nstates <- 32 # Set nstates at 32
        if(length(states) < 17) nstates <- 16 # Set nstates at 16
    } else { # otherwise
        nstates <- 8 # Set nstates at 8 (default)
    }
    if(max(nchar(rownames(clad.matrix$matrix))) > 32) { # Ifthe longest taxon name exceeds 32 characters
        taxname.length <- max(nchar(rownames(clad.matrix$matrix)))+1 # Set maximum taxon name length as that of maximum taxon
    } else { # Otherwise
        taxname.length <- 32 # Just set it at 32 (TNT Ddefault)
    }
    if(head == "") { # If header text is absent:
        if(nstates == 8) toplines <- c("xread", paste(nchars, ntaxa, sep=" ")) # Make topline with no header and no states (8 or less)
        if(nstates == 16) toplines <- c("nstates num 16;", "xread", paste(nchars, ntaxa, sep=" ")) # Make topline with no header and states (greater than 8, but less than 17)
        if(nstates == 32) toplines <- c("nstates num 32;", "xread", paste(nchars, ntaxa, sep=" ")) # Make topline with no header and states (greater than 16, but less than 33)
    } else {
        if(nstates == 8) toplines <- c("xread", paste("'", paste(head, collapse="\n"), "'", sep=""), paste(nchars, ntaxa, sep=" ")) # Make topline with header and no states (8 or less)
        if(nstates == 16) toplines <- c("nstates num 16;", "xread", paste("'", paste(head, collapse="\n"), "'", sep=""), paste(nchars, ntaxa, sep=" ")) # Make topline with header and states (greater than 8, but less than 17)
        if(nstates == 32) toplines <- c("nstates num 32;", "xread", paste("'", paste(head, collapse="\n"), "'", sep=""), paste(nchars, ntaxa, sep=" ")) # Make topline with header and states (greater than 16, but less than 33)
    }
    if(taxname.length > 32) toplines <- c(paste("taxname +", taxname.length, ";", sep=""), toplines) # Add extra line for extra-long taxon names
    taxa <- rownames(mat) # Get taxon names
    spaces.toadd <- (max(nchar(taxa))-nchar(taxa))+2 # Get number of spaces to add to character names
    for(i in 1:length(taxa)) taxa[i] <- paste(taxa[i], paste(rep(" ", spaces.toadd[i]), collapse=""), sep="") # Add spaces to taxon names
    # Convert alphabetic characters to numerical values:
    mat <- gsub("10", "A", mat)
    mat <- gsub("11", "B", mat)
    mat <- gsub("12", "C", mat)
    mat <- gsub("13", "D", mat)
    mat <- gsub("14", "E", mat)
    mat <- gsub("15", "F", mat)
    mat <- gsub("16", "G", mat)
    mat <- gsub("17", "H", mat)
    mat <- gsub("18", "I", mat)
    mat <- gsub("19", "J", mat)
    mat <- gsub("20", "K", mat)
    mat <- gsub("21", "L", mat)
    mat <- gsub("22", "M", mat)
    mat <- gsub("23", "N", mat)
    mat <- gsub("24", "O", mat)
    mat <- gsub("25", "P", mat)
    mat <- gsub("26", "Q", mat)
    mat <- gsub("27", "R", mat)
    mat <- gsub("28", "S", mat)
    mat <- gsub("29", "T", mat)
    mat <- gsub("30", "U", mat)
    mat <- gsub("31", "V", mat)
    mat <- gsub("32", "W", mat)
    mat <- gsub("33", "X", mat)
    mat <- gsub("34", "Y", mat)
    mat <- gsub("35", "Z", mat)
    if(length(grep("&", mat)) > 0) { # If there are polymorphic characters
        for(i in grep("&", mat)) { # For each polymorphic character
            mat[i] <- paste("[", gsub("&", "", mat[i]), "]", sep="") # Convert to tnt format
        }
    }
    if(length(grep(TRUE, is.na(mat))) > 0) { # If missing characters are present
        missing <- grep(TRUE, is.na(mat)) # List missing values
        mat[missing] <- "?" # Replace with question mark
    }
    mtrx <- vector(mode="character") # Character matrix storage vector
    for(i in 1:ntaxa) { # For each taxon
        mtrx[i] <- paste(taxa[i], paste(mat[i, ], collapse=""), sep="") # Store taxon lines in tnt format
    }
    ord <- gsub("unord", "-", ord) # Replace UNORD with - for non-additive characters
    ord <- gsub("ord", "+", ord) # Replace ORD with + for additive characters
    ord.wts <- cbind(ord, wts) # Combine into single character block
    ord.wts <- paste(ord.wts[, 1], ord.wts[, 2], sep="[/") # Separate a la tnt format
    if(length(grep(TRUE, wts == 0)) > 0) { # Case if zero weight (inactive characters) are present
        for(i in grep(TRUE, wts == 0)) ord.wts[i] <- gsub("\\[/0", "]/1", ord.wts[i]) # Make zero weights into inactive characters
    }
    charnos <- c(1:length(ord.wts))-1 # Get character numbers (starting from 0)
    ord.wts <- cbind(ord.wts, charnos) # Add character numbers to character block
    ord.wts <- paste(ord.wts[, 1], ord.wts[, 2], sep=" ") # Convert into vector
    ord.wts <- paste(ord.wts, collapse=" ") # Convert into single text string
    ccode <- paste("ccode ", ord.wts, " ;", sep="") # Make character block
    # Everything above is as in write.tnt(); below is the novel part:
    if(ntaxa <= 24) { # If there are few enough taxa for an exact solution
        anal <- c("collapse [;", "ienum;") # Use implicit enumeration
    } else { # Otherwise
        anal <- c("rseed*;", "hold 999;", "xmult=rss fuse 50 drift 50 ratchet 50;", "mult 50 =tbr drift;", "tsave scratch.tre;", "save;", "tsave /;") # First iteration with new tech
        anal <- c(anal, rep(c("rseed*;", "hold 999;", "xmult=rss fuse 50 drift 50 ratchet 50;", "mult 50 =tbr drift;", "tsave scratch.tre +;", "save;", "tsave /;"), 19)) # Iterations 2-20
        anal <- c(anal, c(paste("hold ", n.trees, ";", sep=""), "shortread scratch.tre;", "bbreak=tbr;")) # Read in new technology trees and finish with a heuristic (tbr) search for all mpts
    }
    out.file <- strsplit(strsplit(file.name, "/")[[1]][length(strsplit(file.name, "/")[[1]])], "\\.")[[1]][1] # Get stripped file name for use in export lines
    search.lines <- c(toplines, mtrx, ";", ccode, anal)
    mrp.line <- vector(mode="character")
    for(i in 1:n.iterations) {
        mrp.line <- c(mrp.line, c(search.lines, "mrp;", paste("export ", out.file, "mrp_", i, ".nex;", sep="")))
    }
    tnt.file <- c("mxram 4096;", mrp.line, "proc/;")
    write(tnt.file, file.name)
}
