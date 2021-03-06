% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Trees2MPTsAndStrict.R
\name{Trees2MPTsAndStrict}
\alias{Trees2MPTsAndStrict}
\title{Converts TNT tree output to ape format}
\usage{
Trees2MPTsAndStrict(filename)
}
\arguments{
\item{filename}{The filename of the TNT output file to read in.}
}
\value{
\item{mpts}{A phylo or multi.phylo object containing the most parsimoniou tree(s).}
\item{strict}{A phylo object containing the strict component consensus tree.}
}
\description{
Converts combined most parsimonious and strict component consensus trees from TNT output to ape format.
}
\details{
This is a function I use to process TNT (Goloboff and Catalano 2016) tree output from re-analysing source matrices, i.e., the MPT(s) and SC links you will see on my \href{http://www.graemetlloyd.com/matr.html}{site}. It assumes that something like the following block of code was used in TNT:

\preformatted{bbreak=tbr;
nelsen*;
export -MYFILENAMEtntmpts_plus_strict.nex;}

In other words, after a final round of tree bisection reconnection (bbreak=tbr) the full set of most parsimonious trees (mpts) will be in memory and to these a strict component consensus will be appended (nelsen*) prior to writing out (export) just the trees (-) in TNT's custom format. NB: This is not quite Newick or NEXUS tree format and hence this function exists to read the TNT format and convert this into ape (Paradis and Schliep 2019) format, whilst also separating out the mpts from the strict component consensus ready for further processing.
}
\examples{

# An example line may look something like this:
#x <- Trees2MPTsAndStrict("tntoutput.tre")

# The mpts can be accessed with:
#x$mpts

# And the strict component consensus with:
#x$strict

}
\references{
Goloboff, P. A. and Catalano, S. A., 2016. TNT version 1.5, including a full implementation of phylogenetic morphometrics. \emph{Cladistics}, \bold{32}, 221-238.

Lloyd, G. T., Bapst, D. W., Friedman, M. and Davis, K. E., 2016. Probabilistic divergence time estimation without branch lengths: dating the origins of dinosaurs, avian flight, and crown birds. \emph{Biology Letters}, \bold{12}, 20160609.

Paradis, E. and Schliep, K., 2019. ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. \emph{Bioinformatics}, \bold{35}, 526–528.
}
\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}
