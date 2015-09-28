#' Counts the changes in a series of time bins
#' 
#' Given a vector of dates for a series of time bins and another for the times when a character change occurred will return the total number of changes in each bin.
#' 
#' Calculates the total number of evolutionary changes in a series of time bins. This is intended as an internal function for rate calculations, but could be used for other purposes (e.g., counting any point events in a series of time bins).
#' 
#' @param xmlwd Text.
#' @param exclude.list Text.
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
#' @export GetPseudoIndependentPhyloFiles
GetPseudoIndependentPhyloFiles <- function(xmlwd, exclude.list) {
    
    # NB: Is "Pseudo" as taxonomic reconciliation needed to be absolutely sure of which data sets to take as independent
    
    # Set working directory:
    setwd(xmlwd)
    
    # Get list of XML files:
    file.list <- list.files()
    
    # Exclude supertrees and stratigraphic data:
    file.list <- setdiff(file.list, exclude.list)
    
    # Create empty character vectors:
    statements <- filenames <- vector(mode="character")
    
    # Create empty numeric vectors:
    Ntax <- Nchar <- vector(mode="numeric")
    
    # For each file:
    for(i in 1:length(file.list)) {
        
        # Read in file:
        x <- readLines(file.list[i])
        
        # Add file name to list:
        filenames[i] <- gsub(".xml", "", file.list[i])
        
        # Add number of taxa to list:
        Ntax[i] <- as.numeric(strsplit(x[grep("<Taxa", x)], "\"")[[1]][2])
        
        # Add number of characters to list:
        Nchar[i] <- as.numeric(strsplit(x[grep("<Morphological |<Molecular ", x)], "\"")[[1]][2])
        
        # Add names to lists:
        names(Ntax)[i] <- names(Nchar)[i] <- filenames[i]
        
        # If there is a parent and/or a sibling:
        if(length(grep("<Parent>|<Sibling>", x)) > 0) {
            
            # Record parent statement:
            if(length(grep("<Parent>", x)) > 0) statements <- c(statements, paste(strsplit(x[grep("<Parent>", x)], "<Parent>|</Parent>")[[1]][2], "is the parent of", filenames[i]))
            
            # Record sibling statement:
            if(length(grep("<Sibling>", x)) > 0) statements <- c(statements, paste(strsplit(x[grep("<Sibling>", x)], "<Sibling>|</Sibling>")[[1]][2], "is the sibling of", filenames[i]))
            
            # If there is no parent or sibling:
        } else {
            
            # Record independent statement:
            statements <- c(statements, paste(filenames[i], "is independent"))
            
        }
        
    }
    
    # Get initial list of linked files:
    linked.files <- unique(sort(unlist(strsplit(statements, " is the parent of | is independent| is the sibling of ")))[duplicated(sort(unlist(strsplit(statements, " is the parent of | is independent| is the sibling of "))))])
    
    # Get list of unique files:
    unique.files <- setdiff(gsub(" is independent", "", statements[grep("is independent", statements)]), linked.files)
    
    # Update linked file list:
    linked.files <- setdiff(gsub(".xml", "", file.list), unique.files)
    
    # Remove independent data sets from statement list:
    if(length(grep("is independent", statements)) > 0) statements <- statements[-grep("is independent", statements)]
    
    # Make matrix of parent statements:
    parentofs <- matrix(unlist(strsplit(statements[grep(" parent of ", statements)], " is the parent of ")), ncol=2, byrow=TRUE)
    
    # Make matrix of sibling statements:
    siblingofs <- matrix(unlist(strsplit(statements[grep(" sibling of ", statements)], " is the sibling of ")), ncol=2, byrow=TRUE)
    
    # Order parent matrix:
    parentofs <- parentofs[order(parentofs[, 1]), ]
    
    # Sort sibling matrix:
    siblingofs <- t(apply(siblingofs, 1, sort))
    
    # For each parent:
    for(i in unique(parentofs[, 1])) {
        
        # Record all children:
        newsiblings <- parentofs[grep(TRUE, parentofs[, 1] == i), 2]
        
        # If there is more than one child:
        if(length(newsiblings) > 1) {
            
            # Add children as siblings to sibling matrix:
            for(j in 2:length(newsiblings)) siblingofs <- rbind(siblingofs, sort(newsiblings[c((j - 1), j)]))
            
        }
        
    }
    
    # Add list of siblings for which a parent statement may be missing:
    checkforparents <- unique(sort(unique(siblingofs)))[grep(TRUE, is.na(match(unique(sort(unique(siblingofs))), parentofs)))]
    
    # For each taxon in the list:
    for(i in checkforparents) {
        
        # Record the child:
        child <- i
        
        # Record the rows the child occurs in:
        rows <- grep(TRUE, siblingofs == i) %% nrow(siblingofs)
        
        # Little correction for last row:
        rows[rows == 0] <- nrow(siblingofs)
        
        # If there are more siblings to find:
        while(length(unique(as.vector(siblingofs[rows, ]))) > length(i)) {
            
            # Update siblings:
            i <- sort(unique(as.vector(siblingofs[rows, ])))
            
            # Update rows:
            for(j in i) rows <- sort(unique(c(rows, grep(TRUE, siblingofs == j) %% nrow(siblingofs))))
            
            # Little correction for last row (again):
            rows[rows == 0] <- nrow(siblingofs)
            
        }
        
        # If parent statement is missing:
        if(length(unique(parentofs[sort(match(i, parentofs[, 2])), 1])) > 0) {
            
            # Add parent statement:
            parentofs <- rbind(parentofs, c(unique(parentofs[sort(match(i, parentofs[, 2])), 1]), child))
            
            # Update ordering of parent matrix:
            parentofs <- parentofs[order(parentofs[, 1]), ]
            
        }
        
    }
    
    # For each potential aunt-niece sibling relationship list:
    for(i in unique(parentofs[, 1])) {
        
        # If there are such relationships:
        if(length(sort(c(match(i, siblingofs[, 1]), match(i, siblingofs[, 2])))) > 0) {
            
            # Record parent siblings (aunts):
            parent.siblings <- setdiff(unique(as.vector(siblingofs[unique(sort(c(match(i, siblingofs[, 1]), match(i, siblingofs[, 2])))), ])), i)
            
            # Record child siblings (nieces):
            children <- parentofs[grep(TRUE, parentofs[, 1] == i), 2]
            
            # For each aunt:
            for(j in parent.siblings) {
                
                # For each niece; add to siblings list:
                for(k in children) siblingofs <- rbind(siblingofs, sort(c(j, k)))
                
            }
            
        }
        
    }
    
    # Empty vector to store redundant parents:
    parents.to.delete <- vector(mode="character")
    
    # For each unique parent:
    for(i in unique(parentofs[, 1])) {
        
        # If present:
        if(!is.na(Nchar[i])) {
            
            # If superceded then add to delete list:
            if(any(apply(rbind(Ntax[parentofs[grep(TRUE, parentofs[, 1] == i), 2]] >= Ntax[i], Nchar[parentofs[grep(TRUE, parentofs[, 1] == i), 2]] >= Nchar[i]), 2, sum) == 2)) parents.to.delete <- c(parents.to.delete, i)
            
            # If data set is absent (excluded or not yet coded):
        } else {
            
            # Add to delete list:
            parents.to.delete <- c(parents.to.delete, i)
            
        }
        
    }
    
    # If there are parent data sets to remove:
    if(length(parents.to.delete) > 0) {
        
        # For each parent:
        for(i in parents.to.delete) {
            
            # Remove from the parent list:
            parentofs <- parentofs[-grep(TRUE, parentofs[, 1] == i), ]
            
        }
        
    }
    
    # Treat remaining parent child relationships as siblings by adding them to the sibling-ofs matrix:
    siblingofs <- rbind(siblingofs, parentofs)
    
    # If there are still dead files (absent through exclusion):
    if(length(grep(TRUE, is.na(match(unique(as.vector(siblingofs)), filenames)))) > 0) {
        
        # List dead files:
        dead.files.to.delete <- unique(as.vector(siblingofs))[grep(TRUE, is.na(match(unique(as.vector(siblingofs)), filenames)))]
        
        # Remove from sibling matrix:
        siblingofs <- siblingofs[-unique(sort(c(match(dead.files.to.delete, siblingofs[, 1]), match(dead.files.to.delete, siblingofs[, 2])))), ]
        
    }
    
    # Create empty vector to store sibling clusters:
    clusters <- vector(mode="character")
    
    # For each sibling:
    for(i in unique(as.vector(siblingofs))) {
        
        # Record sibling:
        child <- i
        
        # Record rows in which sibling occurs:
        rows <- grep(TRUE, siblingofs == i) %% nrow(siblingofs)
        
        # Little correction for last row:
        rows[rows == 0] <- nrow(siblingofs)
        
        # Collapse to unique rows:
        rows <- sort(unique(rows))
        
        # As long as there are more siblings in the cluster:
        while(length(unique(as.vector(siblingofs[rows, ]))) > length(i)) {
            
            # Update sibling list:
            i <- sort(unique(as.vector(siblingofs[rows, ])))
            
            # Update row list:
            for(j in i) rows <- sort(unique(c(rows, grep(TRUE, siblingofs == j) %% nrow(siblingofs))))
            
            # Correction for last row:
            rows[rows == 0] <- nrow(siblingofs)
            
            # Collapse to unique rows:
            rows <- sort(unique(rows))
            
        }
        
        # Add cluster to list:
        clusters <- unique(c(clusters, paste(sort(unique(as.vector(siblingofs[rows, ]))), collapse="%%")))
        
    }
    
    # Create empty vector to store files to delete (dependent and not retained):
    files.to.delete <- vector(mode="character")
    
    # For each cluster:
    for(i in clusters) {
        
        # Break out cluster members:
        members <- strsplit(i, "%%")[[1]]
        
        # Record publication years:
        pub.years <- as.numeric(gsub("[a-z]", "", gsub("inpress", strsplit(as.character(Sys.Date()), "-")[[1]][1], unlist(strsplit(members, "_"))[cumsum(unlist(lapply(strsplit(members, "_"), length)))])))
        
        # Record number of taxa:
        Ntaxa <- Ntax[members]
        
        # Record number of taxa:
        Nchars <- Nchar[members]
        
        # If there is a single data set with the most taxa:
        if(length(grep(TRUE, Ntaxa == max(Ntaxa)) == 1)) {
            
            # Add all other data sets to delete list:
            files.to.delete <- c(files.to.delete, members[-grep(TRUE, Ntaxa == max(Ntaxa))])
            
            # If two or more data sets have equally the most taxa:
        } else {
            
            # Update members to just those with the most taxa:
            members <- members[grep(TRUE, Ntaxa == max(Ntaxa))]
            
            # Update publication years to just those with the most taxa:
            pub.years <- pub.years[grep(TRUE, Ntaxa == max(Ntaxa))]
            
            # Update number of characters to just those with the most taxa:
            Nchars <- Nchars[grep(TRUE, Ntaxa == max(Ntaxa))]
            
            # Update number of taxa just those with the most taxa:
            Ntaxa <- Ntaxa[grep(TRUE, Ntaxa == max(Ntaxa))]
            
            # If a single data set has the most characters (and taxa):
            if(length(grep(TRUE, Nchars == max(Nchars)) == 1)) {
                
                # Add all other files to delete list:
                files.to.delete <- c(files.to.delete, strsplit(i, "%%")[[1]][grep(FALSE, strsplit(i, "%%")[[1]] == members[grep(TRUE, Nchars == max(Nchars))])])
                
                # If two or more data sets have both the most taxa and the most characters:
            } else {
                
                # Update members list:
                members <- members[grep(TRUE, Nchars == max(Nchars))]
                
                # Update publication years:
                pub.years <- pub.years[grep(TRUE, Nchars == max(Nchars))]
                
                # Update number of taxa:
                Ntaxa <- Ntaxa[grep(TRUE, Nchars == max(Nchars))]
                
                # Update number of characters:
                Nchars <- Nchars[grep(TRUE, Ntaxa == max(Nchars))]
                
                # If a single data set has the most taxa and characters and is the most recently published:
                if(length(grep(TRUE, pub.years == max(pub.years)) == 1)) {
                    
                    # Delete all other data sets:
                    files.to.delete <- c(files.to.delete, strsplit(i, "%%")[[1]][grep(FALSE, strsplit(i, "%%")[[1]] == members[grep(TRUE, pub.years == max(pub.years))])])
                    
                    # If there is no grounds to pick a single data set as best:
                } else {
                    
                    # Delete all, but an arbitrarily selected most taxa/characters/recently published data set:
                    files.to.delete <- c(files.to.delete, members[grep(FALSE, strsplit(i, "%%")[[1]] == strsplit(i, "%%")[[1]][grep(TRUE, pub.years == max(pub.years))])][1])
                    
                }
                
            }
            
        }
        
    }
    
    # Record files to keep:
    files.to.keep <- unique(sort(c(unique.files, setdiff(sort(unique(unlist(strsplit(clusters, "%%")))), files.to.delete))))
    
    # Return files to keep:
    return(files.to.keep)
    
}
