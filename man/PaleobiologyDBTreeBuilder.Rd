% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PaleobiologyDBTreeBuilder.R
\name{PaleobiologyDBTreeBuilder}
\alias{PaleobiologyDBTreeBuilder}
\title{Palaeobiology Database Tree Builder}
\usage{
PaleobiologyDBTreeBuilder(
  taxon_nos = NULL,
  taxon_names = NULL,
  original = FALSE,
  interval = NULL,
  extant = "include",
  stopfororphans = TRUE,
  validonly = TRUE,
  returnrank = "3",
  breaker = 100,
  plot.tree = FALSE,
  TimeScale = FALSE
)
}
\arguments{
\item{taxon_nos}{Either a vector of Paleobiology Database taxon numbers that will serve as tips, or a single number that will define the clade requested. Number(s) must match Paleobiology Database records.}

\item{taxon_names}{Either a vector of Paleobiology Database taxon names that will serve as tips, or a single name that will define the clade requested. Name(s) must match Paleobiology Database records.}

\item{original}{Option to be passed to \link{PaleobiologyDBDescendantFinder} or \link{PaleobiologyDBTaxaQuerier}.}

\item{interval}{Option to be passed to \link{PaleobiologyDBDescendantFinder} or \link{PaleobiologyDBTaxaQuerier}.}

\item{extant}{Option to be passed to \link{PaleobiologyDBDescendantFinder} or \link{PaleobiologyDBTaxaQuerier}.}

\item{stopfororphans}{Option to be passed to \link{PaleobiologyDBTaxaQuerier}.}

\item{validonly}{Option to be passed to \link{PaleobiologyDBDescendantFinder}.}

\item{returnrank}{Option to be passed to \link{PaleobiologyDBDescendantFinder}. Default is "3" (species level).}

\item{breaker}{Option to be passed to \link{PaleobiologyDBDescendantFinder} or \link{PaleobiologyDBTaxaQuerier}.}

\item{plot.tree}{Logical whether or not to produce a plot of the resulting tree alongside the output (default is FALSE).}

\item{TimeScale}{Logical indicating whether or not to timescale the tree. (NOT OPERATIONAL YET!)}
}
\value{
A phylo object in ape format with node labels. Note that sometimes multiple labels may be valid for a node and if so these are separated by "_et_".
}
\description{
Using the Paleobiology Database taxonomy builds a phylogenetic tree.
}
\details{
Taxonomies such as those in the Paleobiology Database (\code{paleobiodb.org}; queriable via the API, Peters and McLennen 2016) can also be represented as phylogenetic trees. These can then be used either to aid in reconstructing meta-analytical phylogenies (e.g., Lloyd et al. 2016), or as surrogates for phylogenies where none exist (Soul and Friedman 2015).

The function presented here can either be handed a single higher taxon (e.g., Dinosauria) or a series of (presumably) lower level taxa (genera, species) for which a taxonomic "phylogeny" is desired. Note that all names must appear in the Paleobiology Database, and that it is recommended you actual use the taxonomic number desired to avoid any potential homonym issues, either between animals and plants or lower- and higher-level taxa.

Internally the function will call \link{PaleobiologyDBDescendantFinder} and/or \link{PaleobiologyDBTaxaQuerier} and most options refer to these functions. Note that this means the function will typically take several seconds to run so do not expect an immediate results, even if the desired clade is very small.

The resulting tree is in ape format with node labels, but currently no ages or branch lengths (these are planned future additions).

Note also that the resulting tree should be treated with caution as it will rely solely on what is currently accessible in the Paleobiology Database (\code{paleobiodb.org}). As this database is incomplete it will not necessarily include all valid species (i.e., it logically cannot include taxa not entered in the database). Similarly, the resulting taxonomy is a synthesis of the set of opinions currently entered into the database for the taxa concerned.

Finally, note that typically the resulting tree will contain multiple multifurcations (polytomies) and hence may not be usable in a lot of phylogenetic applications without first resolving these to form a bifurcating tree. (Interested users should consult the functions \code{paleotree::resolveTreeChar} or \code{paleotree::timeLadderTree}.) Caution should also be applied in how this is done, e.g., see Bell and Lloyd (2015).
}
\examples{

# Build a phylogenetic tree out of some famous dinosaurs:
PaleobiologyDBTreeBuilder(taxon_names =
  c("Tyrannosaurus rex", "Triceratops horridus",
  "Brontosaurus excelsus"), plot.tree = TRUE)

}
\references{
Bell, M. A. and Lloyd, G. T., 2015. strap: an R package for plotting phylogenies against stratigraphy and assessing their stratigraphic congruence. \emph{Palaeontology}, \bold{58}, 379-389.

Lloyd, G. T., Bapst, D. W., Friedman, M. and Davis, K. E., 2016. Probabilistic divergence time estimation without branch lengths: dating the origins of dinosaurs, avian flight, and crown birds. \emph{Biology Letters}, \bold{12}, 20160609.

Peters, S. E. and McClennen, M., 2016. The Paleobiology Database application programming interface. \emph{Paleobiology}, \bold{42}, 1-7.

Soul, L. C. and Friedman, M., 2015, Taxonomy and phylogeny can yield comparable results in comparative paleontological analyses. \emph{Systematic Biology}, \bold{64}, 608-620.
}
\seealso{
See also the \code{makePBDBtaxonTree} function in the \code{paleotree} package.
}
\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}
