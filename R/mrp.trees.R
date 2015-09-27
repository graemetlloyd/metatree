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
mrp.trees <- function(trees) {
    ntrees <- length(summary(trees)[, 1]) # How many trees are there?
    strict <- consensus(trees) # Make strict consensus tree
    strict <- root(strict, outgroup=strict$tip.label[1], resolve.root=TRUE) # Root it
    n.furcations <- rle(sort(strict$edge[, 1]))$lengths # Number of branches emanating from each node
    names(n.furcations) <- rle(sort(strict$edge[, 1]))$values # Add node names
    bi.nodes <- as.numeric(names(n.furcations)[grep(TRUE, n.furcations == 2)]) # Get bifurcating nodes
    poly.nodes <- as.numeric(names(n.furcations)[grep(TRUE, n.furcations > 2)]) # Get polytomous nodes
    if(length(grep(TRUE, (Ntip(strict)+1) == bi.nodes)) == 1) bi.nodes <- bi.nodes[-grep(TRUE, (Ntip(strict)+1) == bi.nodes)] # Remove root node if found
    if(length(grep(TRUE, (Ntip(strict)+1) == poly.nodes)) == 1) poly.nodes <- poly.nodes[-grep(TRUE, (Ntip(strict)+1) == bi.nodes)] # Remove root node if found
    all.nodes <- sort(c(bi.nodes, poly.nodes)) # Combination of above
    mrp.all <- matrix(0, nrow=Ntip(strict), ncol=length(all.nodes)) # Make global MRP matrix (initially just for bifurcating nodes)
    rownames(mrp.all) <- sort(strict$tip.label) # Name rows according to taxa
    mrp.weights <- rep(1, length(all.nodes)) # Make global MRP weights (initially just for bifurcating nodes)
    for(i in 1:length(all.nodes)) { # For each bifurcating node in strict consensus
        mrp.all[strict$tip.label[FindDescendants(all.nodes[i], strict)], i] <- 1 # Fill intial MRP matrix
    }
    tips.to.drop <- vector(mode="character") # Tips not relevant to polytomies
    for(i in 1:Ntip(strict)) { # For each tip
        links <- strict$edge[grep(TRUE, strict$edge[, 2] == i), 1] # Start to record nodes along path to tip
        if(links != (Ntip(strict)+1)) { # If not already at root
            while(links[length(links)] != (Ntip(strict)+1)) { # keep going until you are
                links <- c(links, strict$edge[grep(TRUE, strict$edge[, 2] == links[length(links)]), 1]) # Record nodes along path
            }
        }
        if(max(n.furcations[as.character(links)]) == 2) tips.to.drop <- c(tips.to.drop, strict$tip.label[i]) # If a droppable tip add to the list
    }
    if(length(tips.to.drop) > 0) { # Only if there are tips to drop
        pruned.trees <- rmtree(ntrees, Ntip(drop.tip(trees[[1]], tips.to.drop))) # Create new MPT list
        for(i in 1:ntrees) pruned.trees[[i]] <- drop.tip(trees[[i]], tips.to.drop) # without taxa irrelevant to polytomies
        strict <- drop.tip(strict, tips.to.drop) # Remove from consensus too
    }
    nodestogo <- (Ntip(strict)+1):(Ntip(strict)+Nnode(strict)) # Nodes in strict consensus to be ignored in MPTs
    for(i in nodestogo) assign(paste("desc.", i, sep=""), strict$tip.label[FindDescendants(i, strict)])
    desc.groups <- vector(mode="character")
    for(i in 1:ntrees) { # For each MPT
        MPT.nodes <- (Ntip(pruned.trees[[i]])+1):(Ntip(pruned.trees[[i]])+Nnode(pruned.trees[[i]])) # Get node list
        for(j in nodestogo) { # For each node from consensus to be removed from MPT
            MPT.nodes <- MPT.nodes[-match(FindAncestor(get(paste("desc.", j, sep="")), pruned.trees[[i]]), MPT.nodes)] # Remove it
        }
        for(j in MPT.nodes) { # For each novel node in MPT
            desc.groups <- c(desc.groups, paste(sort(pruned.trees[[i]]$tip.label[FindDescendants(j, pruned.trees[[i]])]), collapse="XQ")) # Add result to descendant groups
        }
    }
    new.nodes <- rle(sort(desc.groups))$lengths # Get clade frequencies across all MPTs
    names(new.nodes) <- rle(sort(desc.groups))$values # Label with clades, i.e. taxon clusters
    mrp.all <- cbind(mrp.all, matrix(rep(0, length(new.nodes)*length(mrp.all[, 1])), ncol=length(new.nodes))) # Make space for results in mrp.all
    mrp.weights <- c(mrp.weights, rep(0, length(new.nodes))) # Make space for results in mrp.weights
    for(i in 1:length(new.nodes)) { # For each of the new nodes
        taxa <- strsplit(names(new.nodes[i]), "XQ")[[1]] # Get taxa that make up clade
        col <- (length(mrp.all[1, ])+1)-i # Find appropraite column in MRP matrix and weights
        mrp.all[taxa, col] <- 1 # Insert MRP data
        mrp.weights[col] <- new.nodes[i]/ntrees # Insert MRP results
    }
    out <- list(mrp.all, mrp.weights) # Create results list
    names(out) <- c("mrp", "weights") # Add names
    return(out) # Return results
}
