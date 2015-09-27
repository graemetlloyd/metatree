#' Scrapes the Paleobiology Database
#' 
#' Resolve input name with taxonomy in Paleobiology Database
#' 
#' Deprecated as replaced by API functions.
#' 
#' @param taxon The taxon name.
#' @param taxon.no Alternatively, the taxon number (supercedes taxon name if supplied).
#' @param occurrences Whether or not to also return information on occurrences (if available).
#' @param db The database to use, either "Paleobiodb" or "Fossilworks".
#'
#' @return A matrix.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#' 
#' # Nothing yet
PaleoDBWebScraper <- function(taxon, taxon.no = NULL, occurrences = FALSE, db = "Paleobiodb") {
    
    # Need to set this to avoid "foreign characters" (such as umlauts etc.) from screwing up the code!:
    Sys.setlocale('LC_ALL', 'C')
    
    # Sub function to grab html from PaleoDB:
    PaleoDBScraper <- function(web.address) {
        
        # Set x as NA in case of issues with server requiring repeated scraping attempts:
        x <- NA
        
        # Whilst x is not set (e.g., because of a server error):
        while(length(c(grep("</html>", x), grep("Please select a taxonomic name", x))) == 0) {
            
            # Keep trying to set x:
            try(x <- readLines(web.address), silent = TRUE)
            
        }
        
        # Return web page:
        return(x)
        
    }
    
    # Function to extract web addresses for each taxon given a list in html format from the PaleoDB:
    GetPaleoDBTaxonAddresses <- function(html, db = db) {
        
        # Initial set of taxon.numbers:
        taxon.numbers <- matrix(unlist(strsplit(html[grep("a href=\"\\?a=basicTaxon", html)], "taxon_no=")), byrow=TRUE, ncol=2)[, 2]
        
        # Modify further:
        taxon.numbers <- unlist(apply(matrix(taxon.numbers), 2, FUN=strsplit, split="\""))
        
        # For each line in taxon.numbers:
        for(i in length(taxon.numbers):1) {
            
            # If not a taxon number:
            if(length(grep("[0-9]", strsplit(taxon.numbers[i], "")[[1]])) < nchar(taxon.numbers[i])) {
                
                # Remove non-taxon number:
                taxon.numbers <- taxon.numbers[-i]
                
            }
            
        }
        
        # Turn taxon numbers into web addresses:
        if(db == "Paleobiodb") web.addresses <- paste("http://paleobiodb.org/cgi-bin/bridge.pl?a=basicTaxonInfo&taxon_no=", taxon.numbers, sep = "")
        
        # Turn taxon numbers into web addresses:
        if(db == "Fossilworks") web.addresses <- paste("http://fossilworks.org/cgi-bin/bridge.pl?a=basicTaxonInfo&taxon_no=", taxon.numbers, sep = "")
        
        # Output web addresses:
        return(web.addresses)
        
    }
    
    # Function to extract parent name from a PaleoDB taxon web page:
    GetTaxonName <- function(html) {
        
        # Identify the line with the taxon name:
        taxon.name.line <- html[grep("displayPanelHeader", html)]
        
        # Clean up line so first part relates to the taxon:
        taxon.name.line <- gsub("<span class=\"displayPanelHeader\">|&dagger;|</span>", "", taxon.name.line)
        
        # Case if taxon is stripped of all other information (i.e., no reference for name):
        if(length(strsplit(taxon.name.line, " ")[[1]]) == 1) {
            
            # Add name to list:
            taxon.name <- taxon.name.line
            
        # Case if taxon contains additional information (i.e., species name, reference etc.):
        } else {
            
            # Case if taxon is genus or species:
            if(paste(strsplit(taxon.name.line, "")[[1]][1:3], collapse="") == "<i>") {
                
                # Add name to list:
                taxon.name <- strsplit(taxon.name.line, "<i>|</i>")[[1]][2]
                
                # Find year named line:
                naming.year.line <- unlist(strsplit(strsplit(gsub("<a href=|</a>|>", "", strsplit(taxon.name.line, "<i>|</i>")[[1]][3]), "\"")[[1]], " "))
                
                # If there is a year named:
                if(length(intersect(grep("[0-9]{4}", naming.year.line), grep(TRUE, nchar(naming.year.line) == 4))) > 0) {
                    
                    # Record year named:
                    naming.year <- as.numeric(naming.year.line[intersect(grep("[0-9]{4}", naming.year.line), grep(TRUE, nchar(naming.year.line) == 4))])
                    
                }
                
            # Case if taxon is supra-generic:
            } else {
                
                # Case if taxon is ranked:
                if(strsplit(taxon.name.line, " ")[[1]][2] != "clade") {
                    
                    # Store taxon name:
                    taxon.name <- strsplit(taxon.name.line, " ")[[1]][2]
                    
                    # Find naming year line:
                    naming.year.line <- strsplit(gsub("<a href=|</a>|>", "", taxon.name.line), "\"| ")[[1]]
                    
                    # If year is present:
                    if(length(grep("[0-9]{4}", naming.year.line)) > 0) {
                        
                        # Record year:
                        naming.year <- as.numeric(naming.year.line[max(grep("[0-9]{4}", naming.year.line))])
                        
                    }
                    
                # Case if taxon is unranked clade:
                } else {
                    
                    # Store taxon name:
                    taxon.name <- strsplit(taxon.name.line, " ")[[1]][3]
                    
                    # Find naming year line:
                    naming.year.line <- strsplit(gsub("<a href=|</a>|>", "", taxon.name.line), "\"| ")[[1]]
                    
                    # If year is present:
                    if(length(intersect(grep("[0-9]{4}", naming.year.line), grep(TRUE, nchar(naming.year.line) == 4))) > 0) {
                        
                        # Record year:
                        naming.year <- as.numeric(naming.year.line[intersect(grep("[0-9]{4}", naming.year.line), grep(TRUE, nchar(naming.year.line) == 4))])
                        
                    }
                    
                }
                
            }
            
        }
        
        # Compile output data:
        out <- list(taxon.name, naming.year)
        
        # Add names:
        names(out) <- c("taxon.name", "year.named")
        
        # Return output:
        return(out)
        
    }
    
    # Function to get taxon number from a taxon PaleoDB page:
    GetTaxonNumber <- function(html) {
        
        # Store taxon number:
        taxon.number <- as.numeric(gsub("\">", "", strsplit(html[max(grep("<input type=\"hidden\" name=\"last_taxon\" value=\"", html))], "value=\"")[[1]][2]))
        
        # Return taxon number:
        return(taxon.number)
        
    }
    
    # Function to extract parent name from a PaleoDB taxon web page:
    GetParentName <- function(html) {
        
        # Case if the taxon is assigned:
        if(length(grep("Belongs to|Parent taxon:", html)) > 0) {
            
            # Identify the line with the taxon name:
            parent.name.line <- html[grep("Belongs to|Parent taxon:", html)]
            
            # Strip out some of the surrounding rubbish:
            parent.name.line <- strsplit(parent.name.line, "Belongs to |Parent taxon: | according to ")[[1]][2]
            
            # Set parent name:
            parent.name <- strsplit(gsub("<i>|</i>", "", parent.name.line), "\">|</")[[1]][2]
            
        # Case if taxon is unassigned:
        } else {
            
            # Set parent name as unassigned:
            parent.name <- "Unassigned"
            
        }
        
        # Output parent name:
        return(parent.name)
        
    }
    
    # Function to extract parent name from a PaleoDB taxon web page:
    GetParentNumber <- function(html) {
        
        # Case if the taxon is assigned:
        if(length(grep("Belongs to|Parent taxon:", html)) > 0) {
            
            # Identify the line with the taxon name:
            parent.name.line <- html[grep("Belongs to|Parent taxon:", html)]
            
            # Strip out some of the surrounding rubbish:
            parent.name.line <- strsplit(parent.name.line, "Belongs to |Parent taxon: | according to ")[[1]][2]
            
            # Set parent name:
            parent.number <- as.numeric(strsplit(strsplit(parent.name.line, "taxon_no=")[[1]][2], "\"")[[1]][1])
            
        # Case if taxon is unassigned:
        } else {
            
            # Set parent name as unassigned:
            parent.number <- NULL
            
        }
        
        # Output parent name:
        return(parent.number)
        
    }
    
    # Function to extract subtaxa from a PaleoDB taxon web page:
    GetSubtaxa <- function(html) {
        
        # If there are subtaxa:
        if(length(grep("Subtaxa", html)) > 0) {
            
            # Get subtaxa line:
            subtaxa.line <- html[grep("Subtaxa", html)]
            
            # Modify to help isolate subtaxa:
            subtaxa.line <- strsplit(strsplit(subtaxa.line, "Subtaxa: ")[[1]][2], ", ")[[1]]
            
            # Modify again to help isolate subtaxa:
            subtaxa.line <- matrix(unlist(strsplit(subtaxa.line, "taxon_no=")), ncol=2, byrow=TRUE)[, 2]
            
            # Record subtaxa:
            subtaxa <- sort(gsub("[0-9]|\"|</a>|<i>|</i>|>|</span>|</p>", "", subtaxa.line))
            
        # If there are not subtaxa:
        } else {
            
            # Set subtaxa as "None":
            subtaxa <- "None"
            
        }
        
        # Output subtaxa found:
        return(subtaxa)
        
    }
    
    # Function to get all subtaxa for a particular taxon:
    GetAllSubtaxa <- function(html, db=db) {
        
        # Get taxon number to identify nomen dubia (that cycle back to main page:
        taxon.number <- GetTaxonNumber(html)
        
        # Subfunction to get subtaxa numbers from PaleoDB taxon page:
        GetSubtaxaNumbers <- function(html) {
            
            # If the taxon is supraspecific:
            if(length(grep("Subtaxa: ", html)) > 0) {
                
                # Get subtaxa line:
                subtaxa.line <- strsplit(html[grep("Subtaxa: ", html)], "Subtaxa: ")[[1]][2]
                
                # Case if ther are no sub taxa:
                if(length(grep("<i>none</i>", subtaxa.line)) == 1) {
                    
                    # Create empty vector:
                    taxon.numbers <- vector(mode="numeric")
                    
                # Case if there are subtaxa:
                } else {
                    
                    # Modify to help zero in on taxon numbers:
                    subtaxa.line <- matrix(unlist(strsplit(strsplit(subtaxa.line, ", ")[[1]], "taxon_no=")), ncol=2, byrow=TRUE)[, 2]
                    
                    # Get initial set of taxon numbers:
                    taxon.numbers <- as.numeric(matrix(unlist(strsplit(subtaxa.line, "\"")), ncol=2, byrow=TRUE)[, 1])
                    
                }
                
            # If taxon is a species:
            } else {
                
                # Create empty vector:
                taxon.numbers <- vector(mode="numeric")
                
            }
            
            # Output taxon numbers:
            return(taxon.numbers)
            
        }
        
        # Get initial taxon numbers for subtaxa:
        taxon.numbers <- GetSubtaxaNumbers(html)
        
        # If there are subtaxa:
        if(length(taxon.numbers) > 0) {
            
            # Get taxon web addresses:
            if(db == "Paleobiodb") web.addresses <- paste("http://paleobiodb.org/cgi-bin/bridge.pl?a=basicTaxonInfo&taxon_no=", taxon.numbers, sep="")
            
            # Get taxon web addresses:
            if(db == "Fossilworks") web.addresses <- paste("http://fossilworks.org/cgi-bin/bridge.pl?a=basicTaxonInfo&taxon_no=", taxon.numbers, sep="")
            
            # Make first pass of subtaxa to gather taxon numbers:
            first.pass <- lapply(as.list(web.addresses), PaleoDBScraper)
            
            # Get subtaxa numbers:
            subtaxa.numbers <- setdiff(sort(unlist(lapply(first.pass, GetTaxonNumber))), taxon.number)
            
            # Get taxon web addresses:
            if(db == "Paleobiodb") web.addresses <- paste("http://paleobiodb.org/cgi-bin/bridge.pl?a=basicTaxonInfo&taxon_no=", subtaxa.numbers, sep="")
            
            # Get taxon web addresses:
            if(db == "Fossilworks") web.addresses <- paste("http://fossilworks.org/cgi-bin/bridge.pl?a=basicTaxonInfo&taxon_no=", subtaxa.numbers, sep="")
            
            # As long as there are new subtaxa web addresses to query:
            while(length(web.addresses) > 0) {
                
                # Grab html for subtaxon pages:
                new.pass <- lapply(as.list(web.addresses), PaleoDBScraper)
                
                # Grab new subtaxon numbers:
                new.subtaxa.numbers <- setdiff(unlist(lapply(new.pass, GetSubtaxaNumbers)), c(taxon.number, subtaxa.numbers))
                
                # If there are new subtaxa:
                if(length(new.subtaxa.numbers) > 0) {
                    
                    # Update web addresses:
                    if(db == "Paleobiodb") web.addresses <- paste("http://paleobiodb.org/cgi-bin/bridge.pl?a=basicTaxonInfo&taxon_no=", new.subtaxa.numbers, sep="")
                    
                    # Update web addresses:
                    if(db == "Fossilworks") web.addresses <- paste("http://fossilworks.org/cgi-bin/bridge.pl?a=basicTaxonInfo&taxon_no=", new.subtaxa.numbers, sep="")
                    
                    # Update subtaxa numbers:
                    subtaxa.numbers <- sort(c(subtaxa.numbers, setdiff(unlist(lapply(new.pass, GetSubtaxaNumbers)), c(taxon.number, subtaxa.numbers))))
                    
                # If there are no new subtaxa:
                } else {
                    
                    # Create empty vector that will break loop:
                    web.addresses <- vector(mode="character")
                    
                }
                
            }
            
        # Case if there are no subtaxa:
        } else {
            
            # Set subtaxa numbers as "None":
            subtaxa.numbers <- "None"
            
        }
        
        # Return numbers for all subtaxa:
        return(subtaxa.numbers)
        
    }
    
    # Get occurrences from PaleoDB taxon page:
    GetCollectionNumbers <- function(html, db=db) {
        
        # If there are no occurrences of this taxon:
        if(length(grep("there are no occurrences of", html)) == 1) {
            
            # Record that there are no occurrences of this taxon:
            collection.numbers <- "There are no listed occurrences for this taxon"
            
        # If there are occurrences:
        } else {
            
            # Case if there is a single occurrence:
            if(length(grep("Distribution: found only at ", html)) == 1) {
                
                # Store collection number:
                collection.numbers <- as.numeric(strsplit(strsplit(html[grep("Distribution: found only at ", html)], "collection_no=")[[1]][2], "\"")[[1]][1])
                
            # Case if there are no or multiple ocurrences:
            } else {
                
                # If there are no occurrences (because taxon is too big):
                if(length(grep("<p>Distribution:</p>", html)) == 0) {
                    
                    # Store warning as collection.numbers:
                    collection.numbers <- "There are no listed occurrences for this taxon"
                    
                # Case if there are multiple occurrences:
                } else {
                    
                    # Create vector to store occurrence web addresses:
                    occurrence.web.addresses <- vector(mode="character")
                    
                    # Found occurrence lines:
                    occurrence.lines <- html[c(grep("<p>Distribution:</p>", html):length(html))][grep("&bull;", html[c(grep("<p>Distribution:</p>", html):length(html))])]
                    
                    # Modify occurrence lines:
                    occurrence.lines <- strsplit(paste(occurrence.lines, collapse=", "), ", ")[[1]]
                    
                    # Modify occurrence further:
                    occurrence.lines <- occurrence.lines[grep("<a href=", occurrence.lines)]
                    
                    # For each occurrence:
                    for(i in 1:length(occurrence.lines)) {
                        
                        # Get web address for that occurrence:
                        if(db == "Paleobiodb") occ.address <- paste("http://paleobiodb.org/cgi-bin/bridge.pl?a=displayCollResults", strsplit(strsplit(occurrence.lines[i], "displayCollResults")[[1]][2], "\">")[[1]][1], sep="")
                        
                        # Get web address for that occurrence:
                        if(db == "Fossilworks") occ.address <- paste("http://fossilworks.org/cgi-bin/bridge.pl?a=displayCollResults", strsplit(strsplit(occurrence.lines[i], "collectionSearch")[[1]][2], "\">")[[1]][1], sep="")
                        
                        # Add to list of web addresses:
                        occurrence.web.addresses <- c(occurrence.web.addresses, occ.address)
                        
                    }
                    
                    # Solve issue with spaces:
                    occurrence.web.addresses <- gsub(" ", "%20", occurrence.web.addresses)
                    
                    # Create vector to store collection numbers:
                    collection.numbers <- vector(mode="character")
                    
                    # For each web address:
                    for(i in 1:length(occurrence.web.addresses)) {
                        
                        # Scrape initial web page to create occurrences list:
                        occs.list <- PaleoDBScraper(occurrence.web.addresses[i])
                        
                        # If there are multiple occurrences:
                        if(length(union(grep("There are [0-9]{1} matches", occs.list), grep("There are [0-9]{2} matches", occs.list))) == 1) {
                            
                            # If the number of occurrences spans multiple pages:
                            if(length(grep("View the next 30 matches", occs.list)) > 0) {
                                
                                # Convert occurrences list into a list!:
                                occs.list <- list(occs.list)
                                
                                # While the last item in the list
                                while(length(grep("View the next 30 matches|View the last [0-9]{1} matches|View the last [0-9]{2} matches", occs.list[[length(occs.list)]])) > 0) {
                                    
                                    # Store last set of occurrences as current occurrence list:
                                    current.occ.list <- occs.list[[length(occs.list)]]
                                    
                                    # Find the line with the link to the next page:
                                    link.line <- current.occ.list[grep("View the next 30 matches|View the last ", current.occ.list)]
                                    
                                    # Make the actual link to the next page:
                                    if(db == "Paleobiodb") new.link <- paste("http://paleobiodb.org/cgi-bin/bridge.pl", strsplit(gsub("<b>|</b>", "", link.line), "<a href=\"|\">View")[[1]][2], sep="")
                                    
                                    # Make the actual link to the next page:
                                    if(db == "Fossilworks") new.link <- paste("http://fossilworks.org/cgi-bin/bridge.pl", strsplit(gsub("<b>|</b>", "", link.line), "<a href=\"|\">View")[[1]][2], sep="")
                                    
                                    # Scrape the next page and replace current occurrence list:
                                    current.occ.list <- PaleoDBScraper(new.link)
                                    
                                    # Add current occurrence list to full list:
                                    occs.list[[length(occs.list)+1]] <- current.occ.list
                                    
                                }
                                
                                # Convert list back into a vector:
                                occs.list <- unlist(occs.list)
                                
                            }
                            
                            # Cut list down to just collection lines:
                            occs.list <- occs.list[grep("collection_no=", occs.list)]
                            
                            # Cut down again to collection number parts:
                            occs.list <- matrix(unlist(strsplit(occs.list, "collection_no=")), ncol=2, byrow=TRUE)[, 2]
                            
                            # Cut down again to further zero in on collection numbers:
                            occs.list <- matrix(unlist(strsplit(occs.list, "</a>")), ncol=2, byrow=TRUE)[, 1]
                            
                            # Store collection numbers:
                            if(db == "Paleobiodb") collection.numbers <- c(collection.numbers, matrix(unlist(strsplit(occs.list, ">")), ncol=2, byrow=TRUE)[, 2])
                            
                            # Store collection numbers:
                            if(db == "Fossilworks") collection.numbers <- c(collection.numbers, matrix(unlist(strsplit(occs.list, "\"")), ncol=2, byrow=TRUE)[, 1])
                            
                        # If there is just a single occurrence:
                        } else {
                            
                            # Add collection number to list (if present):
                            if(length(grep("PaleoDB collection", occs.list)) > 0) collection.numbers <- c(collection.numbers, strsplit(occs.list[grep("PaleoDB collection", occs.list)], "PaleoDB collection |:")[[1]][4])
                            
                        }
                        
                    }
                    
                }
                
            }
            
        }
        
        # case if there are no collection numbers:
        if(any(collection.numbers == "There are no listed occurrences for this taxon")) {
            
            # Reiterate warning (just a placeholder):
            collection.numbers <- "There are no listed occurrences for this taxon"
            
        # If there are collection numbers:
        } else {
            
            # Store collection numbers as numbers:
            collection.numbers <- sort(as.numeric(collection.numbers))
            
        }
        
        # Return collection.numbers:
        return(collection.numbers)
        
    }
    
    # Function to extract data from a PaleoDB collection web page:
    GetCollectionData <- function(html) {
        
        # Find when line:
        when.line <- strsplit(html[grep("When:", html)], "\\(|\\)")[[1]]
        
        # Modify when line:
        when.line <- when.line[length(when.line) - 1]
        
        # Get max and min dates:
        max.min <- as.numeric(gsub(" Ma", "", strsplit(when.line, " - ")[[1]]))
        
        # Find where line:
        where.line <- html[grep("Where:", html)]
        
        # Modify where line to close in on coordinates:
        where.line <- strsplit(where.line, "\\(|\\)")[[1]][2]
        
        # Modify where line further to close in on coordinates:
        where.line <- strsplit(where.line, ": paleocoordinates ")[[1]]
        
        # Define present coordinates:
        present.coords <- gsub("&deg;| N| E", "", strsplit(where.line[1], ", ")[[1]])
        
        # Define palaeo coordinates:
        palaeo.coords <- gsub("&deg;| N| E", "", strsplit(where.line[2], ", ")[[1]])
        
        # FOr latitude and then longitude:
        for(i in 1:2) {
            
            # Look for and reformat negative coordinates:
            if(length(grep(" W| S", present.coords[i])) > 0) present.coords[i] <- paste("-", gsub(" W| S", "", present.coords[i]), sep="")
            
            # Look for and reformat negative coordinates:
            if(length(grep(" W| S", palaeo.coords[i])) > 0) palaeo.coords[i] <- paste("-", gsub(" W| S", "", palaeo.coords[i]), sep="")
            
        }
        
        # Convert palaeo coordinates to numeric:
        present.coords <- as.numeric(present.coords)
        
        # Convert present coordinates to numeric:
        palaeo.coords <- as.numeric(palaeo.coords)
        
        # Identify taxon lines:
        taxon.lines <- html[intersect(grep("<td valign=\"top\">", html), grep("<a href", html))]
        
        # As long as there are taxa at the site:
        if(length(taxon.lines) > 0) {
            
            # Modify to zoom in on names:
            taxon.lines <- unlist(strsplit(unlist(strsplit(taxon.lines, ", ")), "a href"))[grep("taxon_n", unlist(strsplit(unlist(strsplit(taxon.lines, ", ")), "a href")))]
            
            # Modify further to zoom in on names:
            taxon.lines <- matrix(unlist(strsplit(taxon.lines, "</a>")), ncol=2, byrow=TRUE)[, 1]
            
            # Modify further to zoom in on names:
            taxon.names <- gsub("<i>|</i>", "", matrix(unlist(strsplit(taxon.lines, "\">")), ncol=2, byrow=TRUE)[, 2])
            
            # Remove all double spaces in taxon names:
            while(length(grep("  ", taxon.names)) > 0) taxon.names <- gsub("  ", " ", taxon.names)
            
            # Create vector to store binomila names:
            binomial <- vector(mode="character")
            
            # For each taxon name:
            for(i in 1:length(taxon.names)) {
                
                # Case if quotes in name:
                if(length(grep("\"", taxon.names[i])) > 0) {
                    
                    # Get number of quotes in name:
                    n.quotes <- nchar(taxon.names[i]) - nchar(gsub("\"", "", taxon.names[i]))
                    
                    # Get number of spaces in name:
                    n.spaces <- nchar(taxon.names[i]) - nchar(gsub(" ", "", taxon.names[i]))
                    
                    # Get number of qualifiers in name:
                    n.qualifiers <- ((nchar(taxon.names[i]) - nchar(gsub("\\? ", "", taxon.names[i]))) / 2) + ((nchar(taxon.names[i]) - nchar(gsub("aff. ", "", taxon.names[i]))) / 5) + ((nchar(taxon.names[i]) - nchar(gsub("cf. ", "", taxon.names[i]))) / 4)
                    
                    # Simple case of a single name in quotes:
                    if(n.spaces == 0) binomial <- rbind(binomial, c("quotes", gsub("\"", "", taxon.names[i]), "", "indet."))
                    
                    # Case of at least two parts to name:
                    if(n.spaces == 1 && n.quotes == 4) binomial <- rbind(binomial, c("quotes", strsplit(gsub("\"", "", taxon.names[i]), " ")[[1]][1], "quotes", strsplit(gsub("\"", "", taxon.names[i]), " ")[[1]][2]))
                    
                    # Case if quotes belong to informal ending:
                    if(n.spaces > 0 && length(grep("\"", strsplit(taxon.names[i], " ")[[1]][length(strsplit(taxon.names[i], " ")[[1]])])) > 0 && n.quotes == 2) {
                        
                        # Split name by quotes:
                        taxon.name.in.parts <- strsplit(taxon.names[i], "\"")[[1]]
                        
                        # Store qualifier for second part of name:
                        qualifier.second.name <- "Informal"
                        
                        # Store second part of name:
                        second.name <- taxon.name.in.parts[2]
                        
                        # Split first part of name by spaces:
                        taxon.name.in.parts <- strsplit(taxon.name.in.parts[1], " ")[[1]]
                        
                        # If a qualifier is present:
                        if(length(taxon.name.in.parts) == 2) {
                            
                            # Store with qualifier for first name:
                            binomial <- rbind(binomial, c(taxon.name.in.parts[1], taxon.name.in.parts[2], qualifier.second.name, second.name))
                            
                        # If no qualifier is present:
                        } else {
                            
                            # Store without qualifier for first name:
                            binomial <- rbind(binomial, c("", taxon.name.in.parts, qualifier.second.name, second.name))
                            
                        }
                        
                    # Case if valid name is in quotes:
                    } else {
                        
                        # Store name with quotes as qualifier for first name:
                        binomial <- rbind(binomial, c("quotes", strsplit(gsub("\"", "", taxon.names[i]), " ")[[1]][1], "", strsplit(gsub("\"", "", taxon.names[i]), " ")[[1]][2]))
                        
                    }
                    
                # Case if formal name:
                } else {
                    
                    # Case if no modification needed:
                    if(length(strsplit(taxon.names[i], " ")[[1]]) <= 2) {
                        
                        # Add single name to binomial with indet. added:
                        if(length(strsplit(taxon.names[i], " ")[[1]]) == 1) binomial <- rbind(binomial, c("", strsplit(taxon.names[i], " ")[[1]][1], "", "indet."))
                        
                        # Add binomial:
                        if(length(strsplit(taxon.names[i], " ")[[1]]) == 2) binomial <- rbind(binomial, c("", strsplit(taxon.names[i], " ")[[1]][1], "", strsplit(taxon.names[i], " ")[[1]][2]))
                        
                    # Case if some kind of qualifier (e.g., cf., aff. etc.) is present:
                    } else {
                        
                        # Split name using spaces:
                        taxon.name.in.parts <- strsplit(taxon.names[i], " ")[[1]]
                        
                        # Is there a qualifier to the first name:
                        if(length(grep("cf.|aff.|\\?", taxon.name.in.parts[1])) > 0) {
                            
                            # Store qualifer:
                            qualifier.first.name <- taxon.name.in.parts[1]
                            
                            # Remove qualifier from name:
                            taxon.name.in.parts <- taxon.name.in.parts[-1]
                            
                        # If there is no qualifier to first name:
                        } else {
                            
                            # Store empty qualifier:
                            qualifier.first.name <- ""
                            
                        }
                        
                        # Is there a qualifier to the second name:
                        if(length(grep("cf.|aff.|\\?", taxon.name.in.parts[2])) > 0) {
                            
                            # Store qualifer:
                            qualifier.second.name <- taxon.name.in.parts[2]
                            
                            # Remove qualifier from name:
                            taxon.name.in.parts <- taxon.name.in.parts[-2]
                            
                        # If there is no qualifier to the second name:
                        } else {
                            
                            # Store empty qualifier:
                            qualifier.second.name <- ""
                            
                        }
                        
                        # Store binomial with qualifier(s):
                        binomial <- rbind(binomial, c(qualifier.first.name, taxon.name.in.parts[1], qualifier.second.name, taxon.name.in.parts[2]))
                        
                    }
                    
                }
                
            }
            
            # Find collection number line:
            collection.number.line <- html[grep("PaleoDB collection", html)]
            
            # Get colletion number:
            collection.number <- strsplit(strsplit(collection.number.line, "PaleoDB collection ")[[1]][2], ":")[[1]][1]
            
            # Compile output into a table:
            out <- cbind(rep(collection.number, nrow(binomial)), binomial, matrix(rep(max.min, nrow(binomial)), ncol=2, byrow=TRUE), matrix(rep(present.coords, nrow(binomial)), ncol=2, byrow=TRUE), matrix(rep(palaeo.coords, nrow(binomial)), ncol=2, byrow=TRUE))
            
            # Add column names to table:
            colnames(out) <- c("Collection.number", "First.name.qualifier", "First.name", "Second.name.qualifier", "Second.name", "Max.age", "Min.age", "Current.lat", "Current.lon", "Palaeo.lat", "Palaeo.lon")
            
        # If there are no taxa at the site:
        } else {
            
            # Return null:
            out <- NULL
            
        }
        
        # Return output:
        return(out)
        
    }
    
    # Case if not using taxon number:
    if(is.null(taxon.no)) {
        
        # Create hypothetical web address for taxon:
        if(db == "Paleobiodb") web.address <- paste("http://paleobiodb.org/cgi-bin/bridge.pl?a=basicTaxonInfo&taxon_name=", gsub(" ", "%20", taxon), collapse="", sep="")
        
        # Create hypothetical web address for taxon:
        if(db == "Fossilworks") web.address <- paste("http://fossilworks.org/cgi-bin/bridge.pl?a=basicTaxonInfo&taxon_name=", gsub(" ", "%20", taxon), collapse="", sep="")
        
    # Case if using taxon number:
    } else {
        
        # Create hypothetical web address for taxon number:
        if(db == "Paleobiodb") web.address <- paste("http://paleobiodb.org/cgi-bin/bridge.pl?a=basicTaxonInfo&taxon_no=", taxon.no, collapse="", sep="")
        
        # Create hypothetical web address for taxon number:
        if(db == "Fossilworks") web.address <- paste("http://fossilworks.org/cgi-bin/bridge.pl?a=basicTaxonInfo&taxon_no=", taxon.no, collapse="", sep="")
        
    }	
    
    # Grab web page:
    x <- PaleoDBScraper(web.address)
    
    # Set variables as NULL to begin:
    occurrences.table <- naming.year <- taxon.number <- parent.name <- parent.number <- NULL
    
    # If database draws a blank:
    if(length(grep("There is no taxonomic or distributional information about", x)) > 0) {
        
        # Store name as "Taxon missing from database":
        taxon.name <- "Taxon missing from database"
        
    # If database has a hit:
    } else {
        
        # Case if mutiples "hits" in database:
        if(length(grep("Please select a taxonomic name", x)) == 1) {
            
            # Special case of single-letter genus name (e.g., T. rex):
            if(length(grep("[A-Z]{1}. ", taxon)) > 0) {
                
                # Get web addresses for each taxon:
                web.addresses <- GetPaleoDBTaxonAddresses(x)
                
                # Output warning with web addresses:
                cat(paste("Multiple equally valid hits found. Check the following addresses manually and resolve to full name:", paste(web.addresses, collapse="\n"), sep="\n"))
                
                # Store name as "Multiple valid hits in database":
                taxon.name <- "Multiple valid hits in database"
                
            # Case of valid name having multiple hits:
            } else {
                
                # Get web addresses for each taxon:
                web.addresses <- GetPaleoDBTaxonAddresses(x)
                
                # Grab web pages for each hit:
                x <- lapply(as.list(web.addresses), PaleoDBScraper)
                
                # Get parent names for each taxon hit:
                parent.names <- lapply(x, GetParentName)
                
                # Get subtaxa for each taxon hit:
                subtaxa.names <- lapply(x, GetSubtaxa)
                
                # Multiple taxa are the same:
                if(all(unique(rle(sort(unlist(subtaxa.names)))$lengths) == length(x)) && length(unique(unlist(parent.names))) == 1) {
                    
                    # Set x as first hit:
                    x <- x[[1]]
                    
                    # Get taxon name:
                    taxon.name <- GetTaxonName(x)$taxon.name
                    
                # Multiple taxa are different:
                } else {
                    
                    # Output warning with web addresses:
                    cat(paste("Multiple equally valid hits found. Check the following addresses manually and resolve to full name:", paste(web.addresses, collapse="\n"), sep="\n"))
                    
                    # Store name as "Multiple valid hits in database":
                    taxon.name <- "Multiple valid hits in database"
                    
                }
                
            }
            
        # Case if single hit in database:
        } else {
            
            # Get taxon name:
            taxon.name <- GetTaxonName(x)$taxon.name
            
            # Store taxon number:
            taxon.number <- GetTaxonNumber(x)
            
        }
        
    }
    
    # Only continue if there is a single valid taxon:
    if(taxon.name != "Multiple valid hits in database" && taxon.name != "Taxon missing from database") {
        
        # Get naming year:
        naming.year <- GetTaxonName(x)$year.named
        
        # Get parent name:
        parent.name <- GetParentName(x)
        
        # Get parent name:
        parent.number <- GetParentNumber(x)
        
        # If looking for occurrences too:
        if(occurrences) {
            
            # If not a species:
            if(length(grep(" ", taxon.name)) == 0) {
                
                # Get subtaxa numbers:
                subtaxa.numbers <- GetAllSubtaxa(x, db)
                
            # If a species:
            } else {
                
                # Create empty vector for subtaxa numbers:
                subtaxa.numbers <- vector(mode="numeric")
                
            }
            
            # Remove cases of no subtaxa:
            subtaxa.numbers <- subtaxa.numbers[-grep(TRUE, subtaxa.numbers == "None")]
            
            # Get all taxon numbers (of active taxon and any subtaxa):
            taxon.numbers <- c(taxon.number, subtaxa.numbers)
            
            # Make web addresses for each taxon:
            if(db == "Paleobiodb") taxon.pages <- paste("http://paleobiodb.org/cgi-bin/bridge.pl?a=basicTaxonInfo&taxon_no=", taxon.numbers, sep="")
            
            # Make web addresses for each taxon:
            if(db == "Fossilworks") taxon.pages <- paste("http://fossilworks.org/cgi-bin/bridge.pl?a=basicTaxonInfo&taxon_no=", taxon.numbers, sep="")
            
            # Get html for each taxon:
            taxon.pages <- lapply(as.list(taxon.pages), PaleoDBScraper)
            
            # Grab collection numbers for all subtaxa:
            collection.numbers <- lapply(taxon.pages, GetCollectionNumbers, db=db)
            
            # Isolate collection numbers as convert to numeric:
            collection.numbers <- sort(as.numeric(setdiff(unique(unlist(collection.numbers)), "There are no listed occurrences for this taxon")))
            
            # Only continue if there are collections:
            if(length(collection.numbers) > 0) {
                
                # Make web addresses for each collection:
                if(db == "Paleobiodb") collection.pages <- paste("http://paleobiodb.org/cgi-bin/bridge.pl?a=basicCollectionSearch&collection_no=", collection.numbers, sep="")
                
                # Make web addresses for each collection:
                if(db == "Fossilworks") collection.pages <- paste("http://fossilworks.org/cgi-bin/bridge.pl?a=basicCollectionSearch&collection_no=", collection.numbers, sep="")
                
                # Get html for each collection:
                collection.pages <- lapply(as.list(collection.pages), PaleoDBScraper)
                
                # Get valid name for each taxon:
                taxon.names <- matrix(unlist(lapply(taxon.pages, GetTaxonName)), ncol=2, byrow=TRUE)[, 1]
                
                # Get collection output:
                collection.output <- lapply(collection.pages, GetCollectionData)
                
                # Create empty occurrences table:
                occurrences.table <- vector(mode="character")
                
                # Combine collection output into single occurrences table:
                for(i in 1:length(collection.output)) occurrences.table <- rbind(occurrences.table, collection.output[[i]])
                
                # Create empty vector to store row numbers to keep:
                rows.to.keep <- vector(mode="numeric")
                
                # For each taxon name:
                for(i in 1:length(taxon.names)) {
                    
                    # Case if taxon name is a species
                    if(length(strsplit(taxon.names[i], " ")[[1]]) == 2) {
                        
                        # Add to rows to keep:
                        rows.to.keep <- sort(unique(c(rows.to.keep, intersect(grep(strsplit(taxon.names[i], " ")[[1]][1], occurrences.table[, "First.name"]), grep(strsplit(taxon.names[i], " ")[[1]][2], occurrences.table[, "Second.name"])))))
                        
                    # Case if taxon name is supraspecific:
                    } else {
                        
                        # Add to rows to keep:
                        rows.to.keep <- sort(unique(c(rows.to.keep, grep(taxon.names[i], occurrences.table[, "First.name"]))))
                        
                    }
                    
                }
                
                # Collapse table to just the rows for the taxon:
                occurrences.table <- occurrences.table[rows.to.keep, ]
                
            # If there are no occurrences of the taxon:
            } else {
                
                # Store warning for output:
                occurrences.table <- "There are no occurrences in the database."
                
            }
            
        }
        
    }
    
    # Reset input taxon name if using taxon number:
    if(!is.null(taxon.no)) taxon <- taxon.name
    
    # Compile output:
    out <- list(taxon, taxon.name, taxon.number, naming.year, parent.name, parent.number, occurrences.table)
    
    # Add names to output:
    names(out) <- c("input.name", "taxon.name", "taxon.number", "year.named", "parent.name", "parent.number", "occurrences.table")
    
    # Return output:
    return(out)
    
}
