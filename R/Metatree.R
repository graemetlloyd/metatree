#' Builds metatree from source data
#'
#' @description
#'
#' Builds a metatree data set from a set of source data.
#'
#' @param MRPDirectory The directory in which the MRP files are to be read from. See details.
#' @param XMLDirectory The directory in which the XML files are to be read from. See details.
#' @param InclusiveDataList A vector of the data sets to include in the metatree. Can be left empty to just read all files in \code{MRPDirectory} and \code{XMLDirectory}.
#' @param ExclusiveDataList A vector of any data sets to exclude from the metatree. Can be left empty if all data sets in \code{MRPDirectory} and \code{XMLDirectory} are valid. (Intended to exclude things like oogenera or footprint analyses, other supertree data sets etc.)
#' @param TargetClade The name of the target clade of the metatree (e.g., "Dinosauria"). OTUs outside of this clade will be pruned.
#' @param HigherTaxaToCollapse Vector of any higher taxa to collapse (e.g., if you are focused on relationships in a stem-group). NB: It is very important that these are safely monophyletic or the results will be confounded.
#' @param SpeciesToExclude Vector of any individual species to be excluded from the final metatree. Intended to deal with problematic taxa, for example, the dinosaurs Eshanosaurus and Ricardoestesia.
#' @param MissingSpecies What to do with species assigned to the target clade, but not present in the source data. Options are: "exclude" (excludes these missing species; the default and safest option), "genus" (include those species in a genus-level polytomy if the genus is sampled in the source data), and "all" (every species assigned to the target clade will be included). Note that neither "genus" or "all" should be used without careful checks of the taxonomy.
#' @param Interval If restricting the sample to a specific interval of geologic time then use this option (passed to \link{PaleobiologyDBChildFinder} which should be consulted for formatting). Default is NULL (no restriction on ages of tips to be included).
#' @param VeilLine A logical indicating whether to remove older data sets that do not increase taxonomic coverage (TRUE; the default and recommended) or not (FALSE). See Lloyd et al. (2016) and the details section below for more information.
#' @param IncludeSpecimenLevelOTUs A logical indicating whether specimen-level OTUs should (TRUE; the default) or should not (FALSE) be included in the metatree. See details.
#' @param BackboneConstraint The file name of one of the source data sets to be used as a backbone constraint (will enforce topology in final metatree but allows taxa not in topology to fall out inside the constraint). This is not required and the default (NULL) will mean no constraint is applied. See details for more information.
#' @param MonophylyConstraint The file name of one of the source data sets to be used as a monophyly constraint (will enforce topology in final metatree and forces taxa not in topology to fall outside the constraint). This is not required and the default (NULL) will mean no constraint is applied. See details for more information.
#' @param RelativeWeights A numeric vector of four values (default \code{c(1, 1, 1, 1)}) giving the respective weights to use for: 1) the input weights (the weights read in from the source MRP files), 2) the publication year weights (from equation 1 in the supplement of Lloyd et al. 2016), 3) the data set dependency weights (1 / the number of "sibling" data sets; see Lloyd et al. 2016), and 4) the within-matrix weights of individual clades (1 / number of conflicitng clades). Zeroes exclude particular weighting types. E.g., to only use input weights use \code{c(1, 0, 0, 0)}.
#' @param WeightCombination How to combine the weights above. Must be one of either "product" or "sum". Note product will exclude zero weight values to avoid zero weight output. E.g., if only using input weights the result of combining weights will not be all zeroes simply because the other types of weight are set at zero.
#' @param ReportContradictionsToScreen Logical indicating whether or not to print any taxonomy-phylogeny contradictions found to the screen. These can aid checking for congruence betwen taxonomy and phylogeny, i.e., they inform the user on whether either the Paleobiology Database or the metadata might need amending.
#'
#' @details
#'
#' \bold{Introduction}
#'
#' Broadly speaking this function is an implementation and extension of the approach to generating composite phylogeneic trees laid out in Lloyd et al. (2016), which itself builds on Lloyd et al. (2008), namely the "metatree" approach. Metatrees are most comparable to formal supertrees but differ in that instead of published trees (figures in source publications) the input data are the original character-taxon matrices and/or sequence alignments. Thus metatrees can be considered superior to formal supertrees if you want to: 1) standardise the way input data are analysed, 2) choose the optimality criterion applied to inference instead of being forced to use whatever the original study used, 3) include non-focal species that may have been removed from published figures (improving taxon overlap), and 4) more properly incorporate phylogenetic uncertainty rather than being restricted to the use of consensus topologies.
#'
#' \bold{Input values}
#'
#' \emph{Formatting input data}
#'
#' The implementation here assumes this data set reanalysis has already been performed, and the results have been encoded in NEXUS (Maddison et al. 1997) format using Matrix Representation with Parsimony (MRP; Baum 1992; Ragan 1992). Lloyd et al. (2016) further suggested that such MRP should represent every biparition present across a sample of trees (all most parsimonious trees, all trees in a posterior sample). An example of how this should be formatted is shown below:
#'
#' \preformatted{#NEXUS
#'
#' BEGIN DATA;
#'   DIMENSIONS  NTAX=4 NCHAR=3;
#'   FORMAT SYMBOLS= " 0 1" MISSING=? GAP=- ;
#' MATRIX
#'
#' Ancilla      000
#' Turrancilla  011
#' Ancillista   101
#' Amalda       111
#' ;
#' END;
#'
#' BEGIN ASSUMPTIONS;
#'   OPTIONS  DEFTYPE=unord PolyTcount=MINSTEPS ;
#' END;}
#'
#' Note that all characters must be either zero or one (i.e., the function cannot currently deal with Purvis MRP; Purvis 1995). Each MRP file should also have a corresponding metadata file in a specific format expressed as XML. Multiple examples of this format are available at \href{http://www.graemetlloyd.com/matr.html}{graemetlloyd.com} and a much more detailed description of the XML structure can be found in the help file for the \link{ReadMetatreeXML} function, but a simple version is also shown below. Note that this corresponds to the MRP file above.
#'
#' \preformatted{<?xml version="1.0" standalone="yes"?>
#' <SourceTree>
#'   <Source>
#'     <Author>
#'       <List>Michaux, B.</List>
#'     </Author>
#'     <Year>1989</Year>
#'     <Title>Cladograms can reconstruct phylogenies: an example from the fossil record</Title>
#'     <Journal>Alcheringa</Journal>
#'     <Volume>13</Volume>
#'     <Pages>21-36</Pages>
#'     <Booktitle/>
#'     <Publisher/>
#'     <City/>
#'     <Editor/>
#'   </Source>
#'   <Taxa number="4">
#'     <List recon_name="Ancilla" recon_no="10760">Ancilla</List>
#'     <List recon_name="DELETE" recon_no="-1">Turrancilla</List>
#'     <List recon_name="Ancillista" recon_no="10763">Ancillista</List>
#'     <List recon_name="Amalda" recon_no="10743">Amalda</List>
#'   </Taxa>
#'   <Characters>
#'     <Molecular/>
#'     <Morphological number="11">
#'       <Type>Osteology</Type>
#'     </Morphological>
#'     <Behavioural/>
#'     <Other/>
#'   </Characters>
#'   <Analysis>
#'     <Type>Maximum Parsimony</Type>
#'   </Analysis>
#'   <Notes>Based on reanalysis of the original matrix.</Notes>
#'   <Filename>Michaux_1989aa</Filename>
#'   <Parent/>
#'   <Sibling/>
#' </SourceTree>}
#'
#' This format is itself a modification of that used by the Supertree Toolkit (Hill and Davis 2014) and contains two pieces of critical information used by the function: 1) how the Operational Taxonomic Units (OTUs) of the source data set should be reconciled to real taxa (the \code{<Taxa>} tag), and 2) how the source data are related to each other (the \code{<Parent>} and \code{<Sibling>} tags). At present the former operates exclusively using the \href{https://paleobiodb.org}{Paleobiology Database}, which is considered as the taxonomic authority. (NB: This means currently the metatree approach is better suited to extinct-only or total-group inference, but extant taxa are also included there.) For each OTU both the formal name (\code{recon_name}) as spelled in the database, spaces should be replaced with underscores) and its' associated database number should be provided (\code{recon_no}). For example, the OTU name "Trex" would be reconciled to "Tyrannosaurus_rex" and the number "54833". These numbers can be found by looking up individual taxa in the database and referring to the URL (e.g., \href{https://paleobiodb.org/classic/basicTaxonInfo?taxon_no=54833}{54833}).
#'
#' \emph{Input directories}
#'
#' Once these files are prepared for use they should be placed in \emph{separate} directories (one for MRP files and one for XML files) and the paths to these folders are the primary inputs to the function (as \code{MRPDirectory} and \code{XMLDirectory}, respectively). Note that aside from the file extension (\code{.nex}, \code{.xml}) filenames must match perfectly between corresponding MRP and XML files. For example, \code{Eppes_etal_2005amrp.nex} in the MRP directory must have a corresponding \code{Eppes_etal_2005a.xml} file in the XML directory. (Note the "mrp.nex" ending.)
#'
#' \emph{Data set inclusion and exclusion}
#'
#' By default the function will use all MRP and XML files found in the two directories, but the user may also wish to supply a limited list of data sets instead. For only data to include the user should use \code{InclusiveDataList} and a vector of strings corresponding to the file names (without extension). For example, \code{InclusiveDataList = c("Eppes_etal_2005a", "White_et_Pinkman_2008a")}. \code{ExclusiveDataList} works the same way but in reverse and might be preferable if the data sets to exclude are fewer than those to include.
#'
#' \emph{Target clade}
#'
#' As generally speaking not all taxa found in the source trees will want to be included in the final analysis it is necessary to supply a target clade for the function to use. Again, this depends on the Paleobiology Database and so this must be a valid taxon there and spelled correctly (e.g., \code{TargetClade = "Dinosauria"}). In practice any species not assigned to this taxon will be excluded from the final tree, so in practice you may wish to use a slightly more inclusive taxon (e.g., "Dinosauromorpha" rather than "Dinosauria"). In theory if you wanted to include all taxa you could supply the top of the taxonomic hierarchy, i.e., \code{TargetClade = "Life"}.
#'
#' \emph{Taxonomic options}
#'
#' It may be that a more specific complement of taxa is desired and there are multiple ways to achieve this. By default the function will only include species present in the sample (species-level reconciliations amongst the source data sets; \code{MissingSpecies = "species"}). This can be reduced or expanding depending on what is required.
#'
#' More species can be excluded by changing the \code{MissingSpecies} option to something else. Additional taxa can be included if you use "genus" (all species in the Paleobiology Database assigned to a sampled genus are also included, emerging in a polytomy from that genus) or "all" (all species assigned to the target clade are included, again emerging from polytomies for each corresponding supraspecific taxon nested inside the target clade). These options represent a trade-off between taxonomic coverage and precision. If in doubt it is recommended that the user stick to the default here as more inclusive options are also likely to include problematic or obsolete taxa excluded from source data sets for good reasons, or inaccurate placements because the current taxonomic synthesis needs updating.
#'
#' Species can also be individually excluded and in multiple ways. A simple means is by using the \code{SpeciesToExclude} option where the user may supply a vector of names to exclude from the final analysis. For example, \code{SpeciesToExclude = c("Tyrannosaurus_rex", "Triceratops_horridus")}. Again, underscores should separate genus and species and names must match the Paleobiology Database version. However, in practice it may be that the desired species to exclude really represent a clade (supraspecific taxon) and hence it would make more sense to provide a single name instead of many. This can be done with the \code{HigherTaxaToCollapse} option. Note that this will not remove the taxa completely but replace them with a single OTU that will appear in ALL CAPS in the final tree. For example, if the desired target clade is Dinosauria exclusive of crown-birds you could use \code{HigherTaxaToCollapse = c("Neornithes")}, which will replace all crown-bird OTUs with the single taxon "NEORNITHES". Note that this option should only be applied if the user is certain that the higher taxon is monophyletic in the Paleobiology Database, otherwise the results may be compromised.
#'
#' \emph{Specifying a sampling interval}
#'
#' Another way of more specifically sampling taxa might be temporal. For example, if the desired tree is Cretaceous dinosaurs only. This can be specified using the \code{Interval} option. This works by using the \link{PaleobiologyDBChildFinder} function to identify valid species assigned to the specified interval. Here both a highest and a lowest interval must be specified, so for our Cretaceous example these would be the same, i.e., \code{Interval = c("Cretaceous", "Cretaceous")}. Note that currently the function only accepts geologic periods and not finer subdivisions. This option should also be used with caution as not all taxa in the database have temporal information (meaning Cretaceous species could be excluded unintentionally simply because their fossil occurrence(s) have not yet been entered into the database). Additionally, specimen-level OTUs (see below) cannot be excluded this way.
#'
#' \emph{Use of a veil line}
#'
#' A common criticism of formal supertrees is that they treat all source data sets equally, regardless of age or quality. However, this is actually an issue of poor implementation rather than a limitation of the approach - good meta-analysis can (and should) differentially weight source data (see discussion on weights below). Another way to deal with this issue is to apply some form of "veil line" where a specified year is used and data sets older than this year excluded from the analysis. However, there are good reasons to not do this a priori. For example, choosing a year would ideally be based on some explicit quantitative criteri(a) that is actually tested. Additionally, excluding older data sets a priori might lead to critical dependency information also being excluded, compromising both data set weighting and pruning (see below).
#'
#' Here the veil line criterion of Lloyd et al. (2016) is applied by default (\code{VeilLine = TRUE}). This assumes that the desired optimality criterion is taxonomic coverage and hence searches for the year corresponding to the last year between then and the present that all possible taxa can be included. To put this another way, for each valid OTU it finds the most recent source data set in which it was included. The oldest of these (the taxon least recently included in a phylogenetic analysis) will set the veil line.
#'
#' Any data set older than this will not appear in the final analysis. However, any information included in such data sets pertaining to non-dependence is retained. In practice the major advantage here is to reduce the amount of data in the final output, but also to exclude older, potentially less accurate data.
#'
#' \emph{Specimen-level OTUs}
#'
#' By default the function assumes the desired taxonomic-level of OTUs is the species level. Aside from not being implemented here, use of higher-levels such as the genus are strongly cautioned against as they can incorporate all kinds of problems into the analysis (e.g., many genera, extinct and extant, are para- or polyphyletic). However, in practice many source data sets of fossil taxa will include OTUs without a valid species name. For example, the dinosaur fossil nicknamed "Dave" known from specimen NGMC91 (\href{https://en.wikipedia.org/wiki/Sinornithosaurus}{Wikipedia page}), which has been included in numerous phylogenetic analyses of theropod dinosaurs. Such specimen-level OTUs can still be included in the analysis by assigning them to the lowest-level taxon to which they can be safely assigned (e.g., genus, family etc.) and then including their specimen number in the name to ensure they are unique. For example, Lloyd et al. (2016) included two specimen-level OTUs assigned to the family Alvarezsauridae and given the names Alvarezsauridae_indet_MPC_100_99and120 and Alvarezsauridae_indet_YPM_1049 to designate them as separate OTUs.
#'
#' In general such specimen-level OTUs are considered to be valuable to many forms of analysis that might be applied to the resulting metatree. For example, time-scaling (they may be the oldest or youngest members of their larger clade) or biogeography (they may represent some form of spatial extreme, e.g., the only member of the clade known from a specifc region). However, they can complicate issues as well. For example, they are not recorded as "proper" taxa in the Paleobiolgy Database hence their validity must be determined by the user. Especially vexatious here is that the fate of many of these specimens is to be given valid names and hence without vigilant checking by the user the same OTU may be inadvertently included twice, first as a specimen-level OTU before it was named and then as its' proper valid name. Thus the option to exclude such taxa (\code{IncludeSpecimenLevelOTUs = FALSE} is offered here. Note that if these do not exist amongst the source data the issue is irrelevant.
#'
#' \emph{Applying constraint(s)}
#'
#' As with any form of phylogenetic inference it may be desirable to specify some form of constraint to restrict the resulting topolog(ies). For example, a molecular scaffold may be desired or multiple data sets might be generated that reflect competing hypotheses of relationships.
#'
#' In most phylogenetic software a constraint tree is specified as a single tree, but here the constraint must be included in the source data. In other words, it must be expressed as an MRP file and XML file just like any other data set. (To build an MRP file from a tree the user should consult the \link{Tree2MRP} function.) This is to limit the many problems that can come from separately specifying a tree, such as ensuring the taxa match up properly. However, it offers additional benefits too. In particular, the MRP encoding means the constraint can represent not just a single tree but a set of trees (e.g., a posterior sample from a Bayesian analysis of molecular data). Thus the result can be limited to a specific set of biparitions without having to specify a single (consensus) tree that would unintentionally allow relationships not found in the original sample.
#'
#' Only a single data set can be specified as a constraint, but two options may be used, one of either: \code{BackboneConstraint} or \code{MonophylyConstraint}. These represent a constraint tree that either allows taxa not included in the constraint to fit anywhere else (backbone) or forcing them to fall outside the constraint (monophyly). Note that if the constraint tree includes all the sampled OTUs then this option is irrelevant - either of \code{BackboneConstraint} or \code{MonophylyConstraint} could be used and the result would be identical.
#'
#' Note that the constraint works by simply upweighting the constraint data set's MRP block over the remaining characters.
#'
#' \emph{Weighting source data}
#'
#' A major failing of many formal supertree approaches is to weight all input data equally by default without any consideration of variation in information quality. This is problematic for multiple reasons, but it's also very easily dealt with as most parsimony-based inference software (the intended inference software for the final MRP matrix) can take weighting information as input and apply it accordingly. Here the function considers four possible types of weighting:
#'
#' \enumerate{
#'   \item Weighting by input values
#'   \item Weighting by year of publication
#'   \item Weighting by non-independence
#'   \item Weighting by data set size
#' }
#'
#' Weighting by input values means applying the weights read in from the original MRP files. If none are found these are considered to be one for every character. In practice such weights could represent, for example, the posterior probability (relative frequency amongst a posterior sample) of each biparition (MRP character). Alternatively they could simply represent some custom weighting scheme favoured by the user. A word of caution with using such weights, however. As the function works taxa can be altered (duplicated, pruned etc.) by the taxonomic reconciliation steps and hence characters can be deleted and thus the weighting potentially (partially) corrupted. The exact effect of this will depend on specific circumstances, but is something the user should be aware of.
#'
#' Weighting by year of publication seems logical on a basic level: we assume that scientific knowledge improves over time. In other words more recent data sets should be upweighted over older data sets. Lloyd et al. (2016) proposed the following equation (slightly modified here) to capture this:
#'
#' \deqn{W_i = 2^(0.5(x_i - t_0))}
#'
#' Where W_i is the weight of the ith data set, \eqn{x_i} is the year of publication of the ith dataset, and \eqn{t_0} is the year of publication of the oldest included data set. In practice this describes a scenario where the weight assigned doubles in value every two years, e.g., a data set from 2000 would be weighted 1/32 the value of a dataset from 2010.
#'
#' Weighting by non-independence is a means of acknowledging the fact that morphological data sets in particular tend not to be independent hypotheses of phylogeny, but rather frequently reuse previously published data sets often with little modification. The function used here can also exclude data sets based on non-independence. For example, if data set A is reused by two later authors (data sets B and C), then typically data set A can be excluded (equivalent of weighting zero). However, data sets B and C have equal claim to data set A's "weight" and hence should be weighted half each. Note that this process is more complex as branching trees of successive reuse emerge, making multiple data sets redundant and sharing the overall initial data set weight over sometimes a very large number of "descendant" data sets. This complexity is automatically handled by the function.
#'
#' Weighting by data set size is complex. Generally speaking, source data sets with more OTUs will lead to more MRP characters and hence greater influence over the resulting super- or metatree. This in of itself is not strictly a problem that requires intervention as more comprehensive data sets are generally of higher value, increasing both taxonomic overlap and total information. However, the metatree approach complicates this issue as it can take \emph{samples} of trees as input and hence MRP characters can grow not just with the number of OTUs but with the amount of phylogenetic uncertainty in the source data set. Here intervention is required as otherwise uncertain data sets will have undue influence, especially as variation \emph{within} subclades will effectively reinforce the existence of that subclade in the first place (the concern that prompted Purvis' MRP modification; Purvis 1995). This issue is dealt with here by searching for characters within data sets that conflict and downweighting these such that they sum to one (the weight assigned to characters that do not conflict). If there is no conflict (all characters are congruent) then they will all be weighted one.
#'
#' Each of the four different weightings will initially be weighted on a zero to one scale, although in practice nothing will be weighted zero as this effectively excludes it from the analysis (except for the data sets made redundant through non-independence). However, this does not mean all different types of weighting will, or should, be applied in the final analysis. For example, if a priori weights are preferred then the other three types will presumably want to be ignored. Similarly, even if including two or more types of weighting it may be that the user will prefer to rely on one more than others. Both what to include and what to favour can be captured with the \code{RelativeWeights} option, which is simply a vector of four numbers corresponding to the four types of weight listed above and in that order. Thus to only include a priori weights everything else would be set to zero (\code{RelativeWeights = c(1, 0, 0, 0)}). Note that here the default (\code{RelativeWeights = c(1, 1, 1, 1)}) is not the recommended option. However, if wanting to account for both non-independence and year of publication but favour the year of publication a larger weight could be used for the latter (e.g., \code{RelativeWeights = c(0, 10, 1, 0)}).
#'
#' If using two or more weights another consideration is how these should be combined in order to produce a single weight for each character. Lloyd et al. (2016) used the product, but here the sum is also offered as an option. This choice is made using the \code{WeightCombination} variable, with the sum as the default. Importantly, if using the product and setting some relative weights to zero does not mean the result will be a zero (instead the zero values are excluded).
#'
#' Final weights may still differ from what the user expects, however, and this is because of two factors: 1) consideration of taxonomic information in inferring relationships and 2) weight limits set by TNT (Goloboff and Catalano 2016), the assumed inference software.
#'
#' Following Lloyd et al. (2016) the metatree approach includes an additional "hidden" source tree based on the current Paleobiology Database taxonomic hierarchy. This is to provide some minimal information on where taxa go in the resulting tree, apply the missing species option (described above) and potentially break ties where no other information is available. However, unlike other approaches this is \emph{not} a constraint and more generally phylogenetic source data is allowed to overrule taxonomy. This is achieved by setting the weight of all taxonomic MRP characters to one and \emph{starting} all phylogenetic characters at a weight of 10 (i.e., minimally an order of magnitude higher than the taxonomic characters).
#'
#' As the resulting MRP matrix is meant for parsimony analysis the most obvious destination software is TNT, which is optimised for fast searching of large data sets, but sets restrictions on the range of weights. Specifically, TNT can only accept input weights between 0.5 and 1000 and at a maximum precision of two decimal places (0.50-1000.00 total range). Thus weights between 10 and whatever the maximum value is are rescaled to fall on to a 10.00 to 1000.00 scale (maintaining the lower weights for taxonomy). (Lloyd et al. 2016 explored a means to extend this by using more states for each MRP character, effectively increasing the weight of an individual character from 0.50 to 31000.00. However, in practice using more character states causes a dramatic slowdown in performance in TNT and so that option is not implemented here.)
#'
#' \emph{Taxonomy-phylogeny contradictions and "chunking"}
#'
#' In an ideal scenario taxonomic and phylogenetic hierarchies will be perfectly congruent (i.e., all taxa will be monophyletic), making the generation of metatrees a simple and fast affair. However, contradictions will arise when input phylogenies incorporate (for example) paraphyletic groups as supraspecific OTUs or the currently favoured Paleobiology Database opinion for a supraspecific taxon is non-monophyletic.
#'
#' These conflicts may be correctable, either in the metadata (editing the XML file(s)) or the database (entering new data into the Paleobiology Database). However, before either can occur the user must be aware of the issues and this is not easily done through manual inspection of the data. For this purpose a series of checks are made to find any contradictions. This can be ignored by the user (\code{ReportContradictionsToScreen = FALSE}) or output to the screen for manual checking later (\code{ReportContradictionsToScreen = TRUE}).
#'
#' Ultimately this information will also be used for another purpose: "chunking" the analysis. As metatrees get larger in size, moving from hundreds to thousands of tips, inference becomes computationally and temporally more expensive. However, given the reliance on taxonomy as an input tree if some clades are found to be monophyletic (no contradictions between taxonomy and phylogeny) then they can safely be broken out into separate analyses, or "chunks". Critically, this can dramatically speed up inference as it is much easier to infer optimal topologies for, say, ten analyses of 100 tips, than one analysis of 1000 tips.
#'
#' Although this chunking isn't yet implemented automatically it can be done manually by the user using the \code{MonophyleticTaxa} output (see below) and the \code{TargetClade} and \code{HigherTaxaToCollapse} input options. As an example, lets imagine we are interested in generating a metatree of Pseudosuchia (crocodile-line archosaurs). We can initially set this as our target clade and run the metatree function. Once this is complete we can check the monophyletic taxa and might find a large subclade is monophyletic (e.g., Eusuchia, the crocodile crown-group). We can save some time in inference if we then use this information to generate two separate metatrees: 1) a Eusuchia tree (using this as the target clade) and, 2) a Pseudosuchian tree with the Eusuchian portion collapsed to a single tip.
#'
#' Note that in practice large subclades may not be monophyletic and although a similar manual "chunking" could still be forced it is not recommended as it would compromise the results.
#'
#' \bold{Operational steps}
#'
#' In practice the function operates in a linear fashion, taking the data through a multi-step process and associated checks, and informing the user with messages as it goes. It can be slow (minutes to hours) depending on the data set size and options applied, but primarily because it makes multiple calls to the Paleobiology Database API as it goes. These are essential to the operation of the function and are not easy to to speed up so the user should factor in run time to their pipelines.
#'
#' Users will also typically find the first time they run it on their data multiple errors will be found. However, these are typically caught by informative messages and easily fixed by modifying the input data. In other words, human error is pretty much inevitable and typographic errors in reconciled names (for example) will usually be found somewhere. However, these are critical to fix to avoid confounded output data. (It is always true that the worst errors are generated by code that doesn't produce any.)
#'
#' [BELOW IS UNFINISHED]
#'
#' Need to add:
#'
#' - Linear list of operations of function
#'
#' - Output explanations
#'
#' @return
#'
#' \item{FullMRPMatrix}{The full MRP matrix in the same format imported by \code{\link[Claddis]{ReadMorphNexus}}. This can be written out to a file in either NEXUS or TNT format using \code{\link[Claddis]{WriteMorphNexus}} or \code{\link[Claddis]{WriteMorphTNT}}.}
#' \item{STRMRPMatrix}{The safe taxonomic reduction MRP matrix in the same format imported by \code{\link[Claddis]{ReadMorphNexus}}. This can be written out to a file in either NEXUS or TNT format using \code{\link[Claddis]{WriteMorphNexus}} or \code{\link[Claddis]{WriteMorphTNT}}. This is the matrix recommended for analysis as it will be smaller and faster than the full version and taxa can still be reinserted later using \code{\link[Claddis]{SafeTaxonomicReinsertion}}.}
#' \item{TaxonomyTree}{The taxonomic hierarchy of the included taxa presented as an ape "phylo" object, with supraspeciifc taxa as node labels.}
#' \item{MonophyleticTaxa}{A vector of taxa which can be considered monophyletic (no phylogenetic data in the sample contradicts the existence of these clades). The intended use of this is to identify smaller subsets of the data that can be analysed separately to "chunk" the metatree process into smaller, faster parts.}
#' \item{SafelyRemovedTaxa}{The results of the safe taxonomic reduction. This is the \code{$str.list} part of the output of \code{\link[Claddis]{SafeTaxonomicReduction}} and can be used to reinsert taxa later with the \code{\link[Claddis]{SafeTaxonomicReinsertion}} function.}
#' \item{RemovedSourceData}{A vector of source data removed throughout the Metatree function. Note that currently the function does not distinguish between the reasons for this (e.g., too many invalid taxa, too few taxa, redundant through non-independence, removed through the veil year process etc.). Importantly it is not therefore safe to remove these data sets from the input as they may still be contributing to the non-independence information.}
#' \item{VeilYear}{The veil year applied (i.e., only data sets this age or younger are included in the output).}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Baum, B. R., 1992. Combining trees as a way of combining data sets for phylogenetic inference, and the desirability of combining gene trees. \emph{Taxon}, \bold{41}, 3-10.
#'
#' Goloboff, P. A. and Catalano, S. A., 2016. TNT version 1.5, including a full implementation of phylogenetic morphometrics. \emph{Cladistics}, \bold{32}, 221-238.
#'
#' Hill, J. and Davis, K. E., 2014. The Supertree Toolkit 2: a new and improved software package with a Graphical User Interface for supertree construction. \emph{Biodiversity Data Journal}, \bold{2}, e1053.
#'
#' Lloyd, G. T., Davis, K. E., Pisani, D., Tarver, J. E., Ruta, M., Sakamoto, M., Hone, D. W. E., Jennings, R. & Benton, M. J., 2008. Dinosaurs and the Cretaceous Terrestrial Revolution. \emph{Proceedings of the Royal Society B}, \bold{275}, 2483-2490.
#'
#' Lloyd, G. T., Bapst, D. W., Friedman, M. and Davis, K. E., 2016. Probabilistic divergence time estimation without branch lengths: dating the origins of dinosaurs, avian flight, and crown birds. \emph{Biology Letters}, \bold{12}, 20160609.
#'
#' Maddison, D. R., Swofford, D. L. and Maddison, W. P., 1997. NEXUS: an extensible file format for systematic information. \emph{Systematic Biology}, \bold{46}, 590-621.
#'
#' Purvis, A., 1995. A modification to Baum and Ragan's method for combining phylogenetic trees. \emph{Systematic Biology}, \bold{44}, 251-255.
#'
#' Ragan, M., 1992. Phylogenetic inference based on matrix representation of trees. \emph{Molecular Phylogenetics and Evolution}, \bold{1}, 113-126.
#'
#' @examples
#'
#' # Local test for Ichthyopterygia:
#' #Metatree(MRPDirectory = "/Users/eargtl/Documents/Homepage/www.graemetlloyd.com/mrp",
#' #  XMLDirectory = "/Users/eargtl/Documents/Homepage/www.graemetlloyd.com/xml",
#' #  TargetClade = "Ichthyopterygia", InclusiveDataList = sort(c(GetFilesForClade("matricht.html"),
#' #  "Bickelmann_etal_2009a", "Caldwell_1996a", "Chen_etal_2014ba", "Chen_etal_2014bb",
#' #  "deBraga_et_Rieppel_1997a", "Gauthier_etal_1988b", "Laurin_et_Reisz_1995a", "Muller_2004a",
#' #  "Reisz_etal_2011a", "Rieppel_et_Reisz_1999a", "Rieppel_et_deBraga_1996a", "Young_2003a")),
#' #  ExclusiveDataList = c("Averianov_inpressa", "Bravo_et_Gaete_2015a", "Brocklehurst_etal_2013a",
#' #  "Brocklehurst_etal_2015aa", "Brocklehurst_etal_2015ab", "Brocklehurst_etal_2015ac",
#' #  "Brocklehurst_etal_2015ad", "Brocklehurst_etal_2015ae", "Brocklehurst_etal_2015af",
#' #  "Bronzati_etal_2012a", "Bronzati_etal_2015ab", "Brusatte_etal_2009ba", "Campbell_etal_2016ab",
#' #  "Carr_et_Williamson_2004a", "Carr_etal_2017ab", "Frederickson_et_Tumarkin-Deratzian_2014aa",
#' #  "Frederickson_et_Tumarkin-Deratzian_2014ab", "Frederickson_et_Tumarkin-Deratzian_2014ac",
#' #  "Frederickson_et_Tumarkin-Deratzian_2014ad", "Garcia_etal_2006a", "Gatesy_etal_2004ab",
#' #  "Grellet-Tinner_2006a", "Grellet-Tinner_et_Chiappe_2004a", "Grellet-Tinner_et_Makovicky_2006a",
#' #  "Knoll_2008a", "Kurochkin_1996a", "Lopez-Martinez_et_Vicens_2012a", "Lu_etal_2014aa",
#' #  "Norden_etal_inpressa", "Pisani_etal_2002a", "Ruiz-Omenaca_etal_1997a", "Ruta_etal_2003ba",
#' #  "Ruta_etal_2003bb", "Ruta_etal_2007a", "Selles_et_Galobart_2016a", "Sereno_1993a", "Sidor_2001a",
#' #  "Skutschas_etal_inpressa", "Tanaka_etal_2011a", "Toljagic_et_Butler_2013a",
#' #  "Tsuihiji_etal_2011aa", "Varricchio_et_Jackson_2004a", "Vila_etal_2017a", "Wilson_2005aa",
#' #  "Wilson_2005ab", "Zelenitsky_et_Therrien_2008a"), MissingSpecies = "exclude",
#' #  BackboneConstraint = "Moon_inpressa", RelativeWeights = c(0, 100, 10, 1))
#'
#' @export Metatree
Metatree <- function(MRPDirectory, XMLDirectory, InclusiveDataList = c(), ExclusiveDataList = c(), TargetClade = "", HigherTaxaToCollapse = c(), SpeciesToExclude = c(), MissingSpecies = "exclude", Interval = NULL, VeilLine = TRUE, IncludeSpecimenLevelOTUs = TRUE, BackboneConstraint = NULL, MonophylyConstraint = NULL, RelativeWeights = c(1, 1, 1, 1), WeightCombination = "sum", ReportContradictionsToScreen = FALSE) {
  
  # DOUBLE CHECK PARENT REPLACEMENT LINE NOW MULTIPLE PARENTS EXIST - SEEMS TO WORK BUT MIGHT NOT.
  # Weights are also super slow (IntraMatrixWeights really?). Can this be sped up somehow? E.g., way STR is.
  
  # DEFO NEED TO IMPROVE NAME CHECKER AFTER RECON COS THAT IS WHERE MOST OF THE LATER ISSUES ARE (E.G. GENUS NAME NUMBER USED FOR SPECIES)
  # FOR HIGHER TAXA TO COLLAPSE HAVE TO ALSO EDIT CONSTRAINT TREES (AND CHECK THEY CAN EVEN MESH!)
  # CHECK FOR SPECIES THAT BELONG TO A GENUS DIFFERENT TO THE ONE IN THEIR NAME!
  # NEED TO CATCH ISSUE WHERE GENUS NUMBER IS USED FOR A SPECIES (HARD TO CHECK SO FAR DUE TO INDETERMINATES CONTINGENCY); I.E., A SAFER CHECK THAT RECON_NAME MATCHES DATABASE NAME FOR RECON_NO
  # NEED SOME TEST THAT HELPS DETERMINE IF MULTIPLE OCCURRENCES OF SAME TAXON AFTER RECONCILIATION IS CORRECT OR AN ERROR
  # CHECK THERE ARE MULTIPLE TAXA PRE-RECONCILIATION
  # CHECK INDETS DO NOT GIVE MULTIPLE MATCHES
  # ADD BLOCK NAMES TO OUTPUT (I.E. SOURCE DATA); COULD BUILD MRP MATRIX THIS WAY? (ADDITIONAL OUTPUT TYPE PROBABLY BEST WAY TO DO IT)
  # THINK ABOUT BETTER WASY TO HANDLE WEIGHTS AS CHARACTERS GET PRUNED/MERGED/DUPLICATED ETC.
  # CHECK INTERACTION BETWEEN COLLAPSED HIGHER TAXA AND EMPTY HIGHER TAXA THAT CURENTLY GET DELETED
  # ADD INCLUDED DATA SETS TO OUTPUT ALONGSIDE REMOVED!
  # EXTANT TAXA INCLUDE/EXCLUDE OPTION?
  # OTU LEVEL, E.G. A GENUS OPTION? FOR OUTPUT THAT IS. MESSY THOUGH. WHAT ABOUT SPECIMENS?
  # MAKE TAXONOTREE OPTIONAL? I.E., AS AN INCLUDED PART OF THE DATA.
  # RECORD WHY EXCLUDED DATA SETS WERE EXCLUDED (CAN ALLOW FASTER WORKING IN FUTURE BY ELIMINATING TOTALLY REDUNDANT DATA SETS
  
  # OPTIONS TO ADD IN FUTURE:
  #
  # 1. Way to do historical metatrees (should be easy as just prune out starting data sets then run as normal - except setting current year to historical year).
  # 2. Allow Purvis coding instead of Baum and Ragan. (Involves dealing with NA issues that would be caused currently.)
  # 3. Proper chunking and maybe even terminal/TNT calls.
  # 4. Species to include option as alternative to species to exclude.
  # 5. Make Safe Taxonomic Reduction optional.
  # 6. Make operable as web site calls as well as/instead of local directory calls. (Needs proper XMLs to be available.)
  
  # Subfunction that gives just MRPs where matrix is still intact (has rows and columns):
  ActiveMRP <- function(MRPList) unname(which(unlist(lapply(MRPList, function(x) prod(dim(x$Matrix)))) > 0))
  
  # Subfunction to make multi-taxon reconciliations unique OTUs:
  SeparateMultiTaxonReconciliations <- function(ListBlock) {
    
    # Find comma rows (multiple taxa in initial reconciliation):
    commarows <- grep(",", rownames(ListBlock$Matrix))
    
    # If there is at least one multiple-taxon reconciliation:
    if(length(commarows) > 0) {
      
      # For each multiple-taxon reconciliation in reverse order (to avoid later rows not matching):
      for(j in rev(commarows)) {
        
        # Get multiple names of reconciliation:
        multiplenames <- strsplit(rownames(ListBlock$Matrix)[j], "%%%%")[[1]]
        
        # Get multiple-taxon numbers:
        multitaxonnumbers <- strsplit(multiplenames[1], ";")[[1]]
        
        # Get multiple-taxon names:
        multitaxonnames <- strsplit(multiplenames[2], ",")[[1]]
        
        # Check data integrity with respect to multiple-taxon values:
        if(length(multitaxonnumbers) != length(multitaxonnames)) stop(paste("Problem with multiple-taxon reconciliation(s) in ", ListBlock$FileName, " (check commas and semi-colons are correct; i.e., of same length).", sep = ""))
        
        # Add new rows at base of matrix:
        ListBlock$Matrix <- rbind(ListBlock$Matrix, matrix(rep(ListBlock$Matrix[j, ], length(multitaxonnumbers)), nrow = length(multitaxonnumbers), byrow = TRUE, dimnames = list(paste(multitaxonnumbers, multitaxonnames, sep = "%%%%"), c())))
        
        # Remove now redundant row from matrix:
        ListBlock$Matrix <- ListBlock$Matrix[-j, , drop = FALSE]
        
      }
      
    }
    
    # Return updated list block:
    return(ListBlock)
    
  }
  
  # Subfunction to find contradicting MRP characters between a string (single character) and a matrix (multiple characters):
  MRPCharacterContradiction <- function(MRPCharacterString, MRPCharacterMatrix) {
    
    # Check MRP string has names and stop and warn user if not:
    if(is.null(names(MRPCharacterString))) stop("MRPCharacterString must have names. Add and try again.")
    
    # Check MRP matrix has names and stop and warn user if not:
    if(is.null(rownames(MRPCharacterMatrix))) stop("MRPCharacterMatrix must have row names. Add and try again.")
    
    # Check MRP string and matrix match in size and stop and warn user if not:
    if(length(MRPCharacterString) != nrow(MRPCharacterMatrix)) stop("MRPCharacterString must have the same length as the number of rows of MRPCharacterMatrix. Check data and try again.")
    
    # Check names of MRP string and matrix match and stop and warn user if not:
    if(!all(sort(names(MRPCharacterString)) == sort(rownames(MRPCharacterMatrix)))) stop("MRPCharacterString names must match row names of MRPCharacterMatrix. Check names and try again.")
    
    # Check only zeroes and ones are coded and stop and warn user if not:
    if(length(setdiff(unique(c(MRPCharacterString, MRPCharacterMatrix)), c("0", "1"))) > 0) stop("Both MRPCharacterString and MRPCharacterMatrix must consist exclusively of the characters \"0\" and \"1\". Check data and try again.")
    
    # Check MRP string contains both zeroe and ones and stop and warn user if not:
    if(length(unique(MRPCharacterString)) < 2) stop("MRPCharacterString must contain both zeroes and ones.")
    
    # Check every MRP matrix column contains both zeroe and ones and stop and warn user if not:
    if(length(unique(as.vector(MRPCharacterMatrix))) < 2) stop("MRPCharacterMatrix must contain both zeroes and ones.")
    
    # Find names corresponding to scores of "0" in the MRP string:
    StringZeroNames <- names(which(MRPCharacterString == "0"))
    
    # Find names corresponding to scores of "1" in the MRP string:
    StringOneNames <- names(which(MRPCharacterString == "1"))
    
    # Find any matrix columns that contradict the string character (both "0" and "1" coded for zero and one matches in MRP string):
    ContradictionColumns <- which(apply(MRPCharacterMatrix, 2, function(x) length(unique(x[StringZeroNames])) == 2 && length(unique(x[StringOneNames])) == 2))
    
    # Return contradicting columns:
    return(ContradictionColumns)
    
  }
  
  # Subfunction to produce intra matrix weights:
  MRPIntraMatrixWeights <- function(MRPMatrix) {
    
    # Get a list of characters that contradict with each character in turn:
    ContradictionList <- lapply(apply(MRPMatrix, 2, as.list), function(x) MRPCharacterContradiction(unlist(x), MRPMatrix))
    
    # Add every character to its' own list:
    ContradictionList <- mapply(function(x, y) c(x, y), x = as.list(1:ncol(MRPMatrix)), y = ContradictionList)
    
    # Build vector of intra matrix weights:
    IntraMatrixWeights <- 1 / unlist(lapply(ContradictionList, function(x) length(unique(unlist(ContradictionList[x])))))
    
    # Return intra matrix weights:
    return(IntraMatrixWeights)
    
  }
  
  # Subfunction to collapse vector of string to single formatted string:
  WriteListAsString <- function(ListOfItems, OxfordComma = TRUE) {
    
    # If not using Oxford comma set final bridge as not having one:
    if(!OxfordComma) FinalBridge <- " and "
    
    # If using Oxford comma set final bridge as having one:
    if(OxfordComma) FinalBridge <- ", and "
    
    # If list is a single item make that the output:
    if(length(ListOfItems) == 1) Output <- ListOfItems
    
    # If list is two items join with a simple and:
    if(length(ListOfItems) == 2) Output <- paste(ListOfItems, collapse = " and ")
    
    # If three or more items format as list with Oxford Comma option included:
    if(length(ListOfItems) > 2) Output <- paste(paste(ListOfItems[1:(length(ListOfItems) - 2)], collapse = ", "), paste(ListOfItems[(length(ListOfItems) - 1):length(ListOfItems)], collapse = FinalBridge), sep = ", ")
    
    # Return output:
    return(Output)
    
  }
  
  # Subfunction to find taxonomy-phylogeny contradictions and turn them into warning messages:
  ListContradictions <- function(TaxonomyMRP, MRPMatrix, ContradictionTaxa, DataSetName) {
    
    # Subfunction to build warning messages:
    BuildWarningMessages <- function(BoundMatrix, ContradictionTaxon, DataSetName) {
      
      # Create four blocks for matching:
      ZeroZero <- ZeroOne <- OneZero <- OneOne <- BoundMatrix
      
      # Fill blocks with zeroes:
      ZeroOne[, 1] <- OneZero[, 2] <- ZeroZero[1:length(ZeroZero)] <- "0"
      
      # Fill blocks with ones:
      ZeroOne[, 2] <- OneZero[, 1] <- OneOne[1:length(OneOne)] <- "1"
      
      # Build empty list of names:
      NamesList <- list()
      
      # Store names where scores are zero and zero:
      NamesList[["ZeroZeroNames"]] <- names(which(apply(BoundMatrix == ZeroZero, 1, all)))
      
      # Store names where scores are zero and one:
      NamesList[["ZeroOneNames"]] <- names(which(apply(BoundMatrix == ZeroOne, 1, all)))
      
      # Store names where scores are one and zero:
      NamesList[["OneZeroNames"]] <- names(which(apply(BoundMatrix == OneZero, 1, all)))
      
      # Store names where scores are one and one:
      NamesList[["OneOneNames"]] <- names(which(apply(BoundMatrix == OneOne, 1, all)))
      
      # Prune down to just the minimum list lengths (the likley problem candidate(s)):
      NamesList <- NamesList[unlist(lapply(NamesList, length)) == min(unlist(lapply(NamesList, length)))]
      
      # Collapse each item of the list to a singel string:
      NamesList <- lapply(NamesList, WriteListAsString)
      
      # Find datasets with second case (apparent taxa outside clade that should be inside):
      CaseOne <- match(c("ZeroZeroNames", "OneZeroNames"), names(NamesList))
      
      # Find datasets with second case (apparent taxa inside clade that should be outside):
      CaseTwo <- match(c("ZeroOneNames", "OneOneNames"), names(NamesList))
      
      # Remove NAs from first case matches:
      CaseOne <- CaseOne[!is.na(CaseOne)]
      
      # Remove NAs from second case matches:
      CaseTwo <- CaseTwo[!is.na(CaseTwo)]
      
      # If first case exists then reformat string with message to user:
      if(length(CaseOne) > 0) NamesList[CaseOne] <- lapply(NamesList[CaseOne], function(x) paste("In ", DataSetName, " the following taxa were found outside ", ContradictionTaxon, " when taxonomy suggests they should be inside: ", x, ". Check data set and/or taxonomy that this is correct.\n", sep = ""))
      
      # If second case exists then reformat string with message to user:
      if(length(CaseTwo) > 0) NamesList[CaseTwo] <- lapply(NamesList[CaseTwo], function(x) paste("In ", DataSetName, " the following taxa were found inside ", ContradictionTaxon, " when taxonomy suggests they should be outside: ", x, ". Check data set and/or taxonomy that this is correct.\n", sep = ""))
      
      # Reformat names list as a vector for output:
      Output <- unname(unlist(NamesList))
      
      # If multiple values then reformat into a single value:
      if(length(Output) > 1) Output <- paste("In ", DataSetName, " one of the following is true:\n", paste(paste(paste(1:length(Output), ". T", sep = ""), gsub(paste("In ", DataSetName, " t| Check data set and/or taxonomy that this is correct.\n", sep = ""), "", Output), rep("\n", length(Output)), sep = ""), collapse = ""), "Check data set and/or taxonomy that this is correct.\n", sep = "")
      
      # Return output:
      return(Output)
      
    }
    
    # Build output into unique values (as can get duplicates):
    Output <- unique(unlist(lapply(as.list(ContradictionTaxa), function(x) lapply(as.list(MRPCharacterContradiction(TaxonomyMRP[, x], MRPMatrix)), function(y) {BoundMatrix <- cbind(TaxonomyMRP[rownames(MRPMatrix), x], MRPMatrix[, y]); BuildWarningMessages(BoundMatrix, x, DataSetName)}))))
    
    # Return output:
    return(Output)
    
  }
  
  # Subfunction to "stretch" a vector of numbers so they fall on the Min-Max range:
  NumberStretcher <- function(X, Min, Max) {
    
    # Check minimum is smaller than maximum value:
    if(Min >= Max) stop("Minimum values exceeds maximum value.")
    
    # Get maximum value:
    MaxVal <- max(X)
    
    # Get minimum value:
    MinVal <- min(X)
    
    # Build multiplication factor:
    MultiplicationFactor <- 1 / ((MaxVal - MinVal) / (Max - Min))
    
    # Build addition factor:
    AdditionFactor <- 10 - (MinVal * MultiplicationFactor)
    
    # Return values streched over the min-max range:
    return((X * MultiplicationFactor) + AdditionFactor)
    
  }

  # Check MRPDirectory is formatted correctly adn stop and warn user if not:
  if(!all(is.character(MRPDirectory)) || length(MRPDirectory) != 1) stop("MRPDirectory must be a single character string indicating the path to the folder containing the MRP files.")
  
  # Check XMLDirectory is formatted correctly adn stop and warn user if not:
  if(!all(is.character(XMLDirectory)) || length(XMLDirectory) != 1) stop("XMLDirectory must be a single character string indicating the path to the folder containing the XML files.")
  
  # Check TargetClade is formatted correctly adn stop and warn user if not:
  if(!all(is.character(TargetClade)) || length(TargetClade) != 1) stop("TargetClade must be a single character string indicating the desired clade the metatree will represent.")
  
  # Check MissingSpecies respresents a valid option:
  if(length(setdiff(MissingSpecies, c("all", "exclude", "genus"))) > 0) stop("MissingSpecies must be one of \"all\", \"exclude\", or \"genus\".")
  
  # Check VeilLine is a logical and stop and warn user if not:
  if(!is.logical(VeilLine)) stop("VeilLine must be a logical (TRUE or FALSE).")
  
  # Check IncludeSpecimenLevelOTUs is a logical and stop and warn user if not:
  if(!is.logical(IncludeSpecimenLevelOTUs)) stop("IncludeSpecimenLevelOTUs must be a logical (TRUE or FALSE).")
  
  # Set defualt of constraint in use to FALSE:
  ConstraintInUse <- FALSE
  
  # Check that there is a maximum of one constraint tree being used and stop and warn user if so.
  # Technically it ought to be possible to do this, but it leaves open some potentially horrendous disasters best avoided for now:
  if(!is.null(BackboneConstraint) && !is.null(MonophylyConstraint)) stop("Cannot currently apply a backbone constraint and a monophyly constraint simultaneously.")
  
  # If backbone constraint is set:
  if(!is.null(BackboneConstraint)) {
    
    # Set constraint in use to TRUE:
    ConstraintInUse <- TRUE
    
    # Check is only a single value and stop and warn user if not:
    if(length(BackboneConstraint) > 1) stop("BackboneConstraint must be a single value. (Cannot apply two backbone constraints simultaneously.)")
    
    # Check is a string and stop and warn user if not:
    if(!is.character(BackboneConstraint)) stop("BackboneConstraint must be a text string. Reformat and try again.")
    
    # Set constraint data set to backbone constraint:
    ConstraintDataSet <- BackboneConstraint
    
    # Set constraint type to backbone:
    ConstraintType <- "backbone"

  }
  
  # If monophyly onstraint is set:
  if(!is.null(MonophylyConstraint)) {
    
    # Set constraint in use to TRUE:
    ConstraintInUse <- TRUE

    # Check is only a single value and stop and warn user if not:
    if(length(MonophylyConstraint) > 1) stop("MonophylyConstraint must be a single value. (Cannot apply two monophyly constraints simultaneously.)")
    
    # Check is a string and stop and warn user if not:
    if(!is.character(MonophylyConstraint)) stop("MonophylyConstraint must be a text string. Reformat and try again.")
    
    # Set constraint data set to monophyly constraint:
    ConstraintDataSet <- MonophylyConstraint
    
    # Set constraint type to monophyly:
    ConstraintType <- "monophyly"
    
  }
  
  # Check VeilLine is a logical and stop and warn user if not:
  if(!is.logical(ReportContradictionsToScreen)) stop("ReportContradictionsToScreen must be a logical (TRUE or FALSE).")
  
  # List of types of resolution that require finding a senior synonym:
  synonyms <- c("corrected to", "misspelling of", "objective synonym of", "obsolete variant of", "recombined as", "replaced by", "subjective synonym of")
  
  # List of types of resolution that require changing reconciliation to DELETE:
  DeletionCategories <- c("nomen dubium", "nomen vanum", "nomen nudum", "nomen oblitum", "invalid subgroup of")
  
  # Print current processing status:
  cat("Reading MRP data...")
  
  # Set working directory as MRP directory:
  setwd(MRPDirectory)
  
  # List MRP files (or just use inclusivedatalist if set):
  MRPFileList <- strsplit(ifelse(length(InclusiveDataList) > 0, paste(setdiff(sort(unique(InclusiveDataList)), sort(unique(ExclusiveDataList))), "mrp.nex", sep = "", collapse = "%%"), paste(setdiff(gsub("mrp\\.nex", "", list.files()), sort(unique(ExclusiveDataList))), "mrp.nex", sep = "", collapse = "%%")), "%%")[[1]]
  
  # If a backbone constraint is used and is a file check it is present in the source data and stop and warn user if not:
  if(!is.null(BackboneConstraint)) if(is.na(match(paste(BackboneConstraint, "mrp.nex", sep = ""), MRPFileList))) stop("Backbone constraint file not found in data. Check name and try again.")
  
  # If a monophyly constraint is used and is a file check it is present in the source data and stop and warn user if not:
  if(!is.null(MonophylyConstraint)) if(is.na(match(paste(MonophylyConstraint, "mrp.nex", sep = ""), MRPFileList))) stop("Monophyly constraint file not found in data. Check name and try again.")
  
  # Check there are four relative weights values and stop and warn user if not:
  if(length(RelativeWeights) != 4) stop("RelativeWeights must consist of exactly four values. Fix and try again.")
  
  # Check relative weights are numeric and stop and warn user if not:
  if(!is.numeric(RelativeWeights)) stop("RelativeWeights must consist of numeric values. Fix and try again.")
  
  # Check there are no negative weights and at least one positive weight is being used and stop and wanr user if not:
  if(sum(RelativeWeights) <= 0 || any(RelativeWeights < 0)) stop("RelativeWeights must include at least one positive value and negative values are not allowed. Fix and try again.")
  
  # Check only a single weight combination value is being used and stop and warn user if not:
  if(length(WeightCombination) != 1) stop("WeightCombination must consist of a single value.")
  
  # Check weight combintion is a valid option and stop and warn if not:
  if(length(setdiff(WeightCombination, c("product", "sum"))) > 0) stop("WeightCombination must be one of \"product\" or \"sum\" only. Check spelling and try again.")
  
  # Read in all MRP files and store in a list (include duplicate headers to store parent sibling info later):
  MRPList <- lapply(lapply(as.list(MRPFileList), Claddis::ReadMorphNexus), function(x) {y <- list(x$Matrix_1$Matrix, x$Matrix_1$Weights, "", "", ""); names(y) <- c("Matrix", "Weights", "FileName", "Parent", "Sibling"); y})
  
  # Set names of MRP files:
  names(MRPList) <- gsub("mrp.nex", "", MRPFileList)
  
  # Find maximum input weight:
  MaximumInputWeight <- max(unlist(lapply(MRPList, function(x) x$Weights)))
  
  # Resacle all input weights zero to one by dividing through by maximum input weight:
  MRPList <- lapply(MRPList, function(x) {x$Weights <- x$Weights / MaximumInputWeight; x})
  
  # Print current processing status:
  cat("Done\nReading XML data...")
  
  # Set working directory as XML (i.e., metadata) directory:
  setwd(XMLDirectory)
  
  # List MRP files (or just use inslusivedatalist if set):
  XMLFileList <- strsplit(ifelse(length(InclusiveDataList) > 0, paste(setdiff(sort(unique(InclusiveDataList)), sort(unique(ExclusiveDataList))), ".xml", sep = "", collapse = "%%"), paste(setdiff(gsub("\\.xml", "", list.files()), sort(unique(ExclusiveDataList))), ".xml", sep = "", collapse = "%%")), "%%")[[1]]
  
  # Check there are no MRPs not listed as XMLs and vice versa (should return empty vector):
  MRPXMLunion <- c(setdiff(gsub("\\.xml", "", XMLFileList), gsub("mrp\\.nex", "", MRPFileList)), setdiff(gsub("mrp\\.nex", "", MRPFileList), gsub("\\.xml", "", XMLFileList)))
  
  # Stop if MRP datasets not listed as XMLs and vice versa:
  if(length(MRPXMLunion) > 0) stop(paste("Datasets do not match (MRP and XML)!:", MRPXMLunion, collapse = " "))
  
  # Read in all XML files and store in a list (reformatting subgenera as GenusSubgenus along the way):
  XMLList <- lapply(as.list(XMLFileList), function(x) {y <- ReadMetatreeXML(x); y$SourceTree$Taxa$TagContents[, "recon_name"] <- gsub("_\\(|\\)", "", y$SourceTree$Taxa$TagContents[, "recon_name"]); y})
  
  # Add names to XML list:
  names(XMLList) <- gsub(".xml", "", XMLFileList)
  
  # Collapse to just pertinent information:
  XMLList <- lapply(XMLList, function(x) {y <- list(); y[["TaxonMatrix"]] <- x$SourceTree$Taxa$TagContents; y[["FileName"]] <- unname(unlist(x$SourceTree$Filename)); y[["Parent"]] <- unname(unlist(x$SourceTree$Parent)); y[["Sibling"]] <- unname(unlist(x$SourceTree$Sibling)); y})
  
  # Find any files that contain duplicated taxon names:
  FilesWithDuplicatedTaxonNames <- names(XMLList)[which(unlist(lapply(XMLList, function(x) any(duplicated(x$TaxonMatrix[, "ListValue"])))))]
  
  # If duplicate names were found stop and warn user:
  if(length(FilesWithDuplicatedTaxonNames) > 0) stop(paste("The following files contain duplicate taxon names: ", paste(FilesWithDuplicatedTaxonNames, collapse = ", "), ". Ensure all taxon names are unique and try again.", sep = ""))
  
  # Find any taxon names that do not match between MRP and XML:
  TaxonMismatches <- mapply(function(x, y) {MRPNames <- rownames(x$Matrix); XMLNames <- y$TaxonMatrix[, "ListValue"]; c(setdiff(MRPNames, XMLNames), setdiff(XMLNames, MRPNames))}, x = MRPList[names(MRPList)], y = XMLList[names(MRPList)])
  
  # Find any files with mismatching taxon names between MRP and XML:
  FilesWithTaxonMismatches <- names(TaxonMismatches)[which(unlist(lapply(TaxonMismatches, function(x) length(x))) > 0)]
  
  # If such files are found then stop and warn user:
  if(length(FilesWithTaxonMismatches) > 0) stop(paste("The following files contain mismatching taxon names between the MRP and XML versions: ", paste(FilesWithTaxonMismatches, collapse = ", "), ". Ensure all taxon names match and try again.", sep = ""))
  
  # Compile any name issues:
  NameIssues <- lapply(XMLList, function(x) {TaxonMatrix <- x$TaxonMatrix; SpacesFound <- c(grep(" ", TaxonMatrix[, "recon_name"]), grep(" ", TaxonMatrix[, "recon_no"]), grep(" ", TaxonMatrix[, "ListValue"])); EmptyValuesFound <- c(which(TaxonMatrix[, "recon_name"] == ""), which(TaxonMatrix[, "recon_no"] == ""), which(TaxonMatrix[, "ListValue"] == "")); RogueNumberCharacters <- setdiff(unique(unlist(strsplit(TaxonMatrix[, "recon_no"], ""))), c(0:9, ";", "-")); RogueNameCharacters <- setdiff(unique(c(unlist(strsplit(TaxonMatrix[, "recon_name"], "")), unlist(strsplit(TaxonMatrix[, "ListValue"], "")))), c(LETTERS, letters, 0:9, "_", ",")); y <- list(SpacesFound, EmptyValuesFound, RogueNumberCharacters, RogueNameCharacters); names(y) <- c("SpacesFound", "EmptyValuesFound", "RogueNumberCharacters", "RogueNameCharacters"); y})
  
  # Find any files with spaces in taxon names:
  FilesWithSpaces <- names(NameIssues)[unlist(lapply(NameIssues, function(x) length(x$SpacesFound))) > 0]
  
  # Find any values with empty values for taxon names:
  FilesWithEmptyValues <- names(NameIssues)[unlist(lapply(NameIssues, function(x) length(x$EmptyValuesFound))) > 0]
  
  # Files with rogue values in the recon number field:
  FilesWithRogueTaxonNumbers <- names(NameIssues)[unlist(lapply(NameIssues, function(x) length(x$RogueNumberCharacters))) > 0]
  
  # Files with rogue values in the name fields:
  FilesWithRogueTaxonNames <- names(NameIssues)[unlist(lapply(NameIssues, function(x) length(x$RogueNameCharacters))) > 0]
  
  # If issues with spaces in names stop and warn user:
  if(length(FilesWithSpaces) > 0) stop(paste("The following files contain spaces in the taxonomic reconciliation (names or numbers): ", paste(FilesWithSpaces, collapse = ", "), ". Remove spaces and try again.", sep = ""))
  
  # If issues with empty names stop and warn user:
  if(length(FilesWithEmptyValues) > 0) stop(paste("The following files contain empty values in the taxonomic reconciliation (names or numbers): ", paste(FilesWithEmptyValues, collapse = ", "), ". Ensure all values are filled and try again.", sep = ""))
  
  # If issues with rogue characters in number field stop and warn user:
  if(length(FilesWithRogueTaxonNumbers) > 0) stop(paste("The following files contain rogue values in the taxonomic reconciliation (numbers): ", paste(FilesWithRogueTaxonNumbers, collapse = ", "), ". Ensure all taxon numbers only include semicolon(s) (the separating character) or dashes (for negative values) and try again.", sep = ""))
  
  # If issues with rogue characters in name field stop and warn user:
  if(length(FilesWithRogueTaxonNames) > 0) stop(paste("The following files contain rogue values in the taxonomic reconciliation (names): ", paste(FilesWithRogueTaxonNames, collapse = ", "), ". Ensure all taxon names are formed from alphanumerics, commas (the separating character) or underscores and try again.", sep = ""))
  
  # Reconcile OTU names with XML version:
  MRPList <- mapply(function(x, y) {rownames(x$Matrix)[unlist(lapply(as.list(y$TaxonMatrix[, "ListValue"]), function(z) which(rownames(x$Matrix) == z)))] <- paste(y$TaxonMatrix[, "recon_no"], y$TaxonMatrix[, "recon_name"], sep = "%%%%"); x$FileName <- y$FileName; if(!is.null(y$Parent)) x$Parent <- y$Parent; if(!is.null(y$Sibling)) x$Sibling <- y$Sibling; x}, x = MRPList[names(MRPList)], y = XMLList[names(MRPList)], SIMPLIFY = FALSE)

  # Print current processing status:
  cat("Done\nChecking for unsampled parents and siblings...")
  
  # Extract parent and sibling names:
  ParentAndSiblingNames <- sort(unlist(lapply(as.list(unique(unname(unlist(lapply(MRPList, '[', c("Parent", "Sibling")))))), function(x) x[nchar(x) > 0])))
  
  # Warn user about any unsampled parents and/or siblings:
  if(length(setdiff(ParentAndSiblingNames, names(MRPList))) > 0) print(paste("The following parents and siblings are not in the sample (check they are correct or add them into the sample): ", paste(setdiff(ParentAndSiblingNames, names(MRPList)), collapse = ", "), sep = ""))
  
  # Print current processing status:
  cat("Done\nFinding initial multiple-taxon reconciliations...")
  
  # Separate out multi-taxon reconcilations:
  MRPList <- lapply(MRPList, SeparateMultiTaxonReconciliations)
  
  # If excluding specimen-level OTUs:
  if(!IncludeSpecimenLevelOTUs) {
    
    # Print current processing status:
    cat("Done\nRemoving specimen-level OTUs...")
    
    # Convert specimen-level OTUs to taxa to DELETE:
    MRPList <- lapply(MRPList, function(x) {RowNamesToDelete <- which(unlist(lapply(strsplit(rownames(x$Matrix), split = ""), function(y) sum(y == "_") > 2))); if(length(RowNamesToDelete) > 0) rownames(x$Matrix)[RowNamesToDelete] <- "0%%%%DELETE"; x})
    
  }
  
  # Print current processing status:
  cat("Done\nRemoving taxa with initial reconciliations of \"DELETE\"...")
  
  # Remove any taxa reconciled as DELETE:
  MRPList <- lapply(MRPList, function(x) {DeleteRows <- which(unlist(lapply(strsplit(rownames(x$Matrix), split = "%%%%"), function(y) y[2])) == "DELETE"); if(length(DeleteRows) > 0) x$Matrix <- x$Matrix[-DeleteRows, , drop = FALSE]; x})
  
  # Prune matrices following deletion:
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- PisaniMRPPrune(Claddis::MakeMorphMatrix(x$Matrix, Weights = x$Weights, ignore.duplicate.taxa = TRUE)); x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})

  # Print current processing status:
  cat("Done\nSearching for and collapsing pre-reconciliation duplicated taxa...")
  
  # Collapse any duplicate taxon names:
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- Claddis::MakeMorphMatrix(x$Matrix, Weights = x$Weights, ignore.duplicate.taxa = TRUE); if(any(duplicated(rownames(y$Matrix_1$Matrix)))) {DuplicateNames <- setdiff(unlist(lapply(strsplit(rownames(y$Matrix_1$Matrix)[duplicated(rownames(y$Matrix_1$Matrix))], split = "%%%%"), '[', 2)), "DELETE"); if(length(DuplicateNames) > 0) cat(paste("\nDuplicate resolved OTU name(s) found in ", x$FileName, ": ", paste(DuplicateNames, collapse = ", "), ". Check this is correct.", sep = "")); y <- CollapseDuplicateTaxonMRP(y)}; x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
  
  # Print current processing status:
  cat("Done\nBuilding initial taxonomy matrix...")
  
  # Create taxonomy matrix to store all taxon resolution data:
  TaxonomyMatrix <- do.call(rbind, strsplit(unique(unname(unlist(lapply(MRPList[ActiveMRP(MRPList)], function(x) rownames(x$Matrix))))), split ="%%%%"))
  
  # Add column names:
  colnames(TaxonomyMatrix) <- c("TaxonNo", "TaxonName")
  
  # Print current processing status:
  cat("Done\nChecking for missing taxon numbers...")
  
  # If any "-1" taxa found stop and tell user:
  if(any(TaxonomyMatrix[, "TaxonNo"] == "-1")) stop(paste("The following taxa have the reconciliation number \"-1\": ", paste(TaxonomyMatrix[TaxonomyMatrix[, "TaxonNo"] == "-1", "TaxonName"], collapse = ", "), sep = ""))
  
  # Print current processing status:
  cat("Done\nBuilding initial Paleobiology Database reconciliation list...")
  
  # Create resolved taxon numbers matrix:
  ResolvedTaxonNumbers <- cbind(unique(TaxonomyMatrix[, "TaxonNo"]), PaleobiologyDBTaxaQuerier(taxon_nos = unique(TaxonomyMatrix[, "TaxonNo"]), interval = NULL))
  
  # Deal with subgenera:
  ResolvedTaxonNumbers[, "TaxonName"] <- gsub(" \\(|\\)", "", ResolvedTaxonNumbers[, "TaxonName"])
  
  # Add column names to first value (input number):
  colnames(ResolvedTaxonNumbers)[1] <- "InputNo"
  
  # If specifying an Interval:
  if(!all(is.null(Interval))) {
    
    # Do same query for just taxa in Interval:
    ResolvedTaxonNumbersInterval <- cbind(unique(TaxonomyMatrix[, "TaxonNo"]), PaleobiologyDBTaxaQuerier(taxon_nos = unique(TaxonomyMatrix[, "TaxonNo"]), interval = Interval))
    
    # Deal with subgenera:
    ResolvedTaxonNumbersInterval[, "TaxonName"] <- gsub(" \\(|\\)", "", ResolvedTaxonNumbersInterval[, "TaxonName"])
    
    # Invert variable so only includes taxa outside Interval:
    ResolvedTaxonNumbersInterval <- ResolvedTaxonNumbers[is.na(ResolvedTaxonNumbersInterval[, "TaxonName"]), ]
    
    # Find any nomen dubia etc. to delete:
    DeleteRows <- which(unlist(lapply(lapply(lapply(as.list(ResolvedTaxonNumbersInterval[, "TaxonValidity"]), match, x = DeletionCategories), sort), length)) > 0)
    
    # If there are deletes then remove them from the matrix:
    if(length(DeleteRows) > 0) ResolvedTaxonNumbersInterval <- ResolvedTaxonNumbersInterval[-DeleteRows, , drop = FALSE]
    
  }
  
  # Print current processing status:
  cat("Done\nChecking taxon names match with database version...")
  
  # Vector to store any failed matches:
  failedmatches <- vector(mode = "character")
  
  # For each initially reconciled name:
  for(i in 1:nrow(TaxonomyMatrix)) {
    
    # Get resolved name:
    resolvedname <- gsub(" ", "_", ResolvedTaxonNumbers[which(ResolvedTaxonNumbers[, "InputNo"] == TaxonomyMatrix[i, "TaxonNo"]), "TaxonName"])
    
    # Get input name:
    inputname <- TaxonomyMatrix[i, "TaxonName"]
    
    # Check names truly match (i.e., deals with case of indets where direct match not possible) and store if not:
    if(resolvedname != inputname && any(is.na(match(strsplit(resolvedname, "_")[[1]], strsplit(inputname, "_")[[1]])))) failedmatches <- c(failedmatches, paste("Input name ", inputname, " does not match database name for corresponding number (", resolvedname, ").", sep = ""))
    
  }
  
  # If there are failed matches:
  if(length(failedmatches) > 0) {
    
    # Provide list to user:
    cat(paste(failedmatches, collapse = "\n"))
    
    # Stop function (will break later otherwise):
    stop("")
    
  }
  
  # Print current processing status:
  cat("Done\nChecking taxon validities...")
  
  # Check for any new kind of resolution (should be empty vector):
  newresolutions <- setdiff(sort(unique(ResolvedTaxonNumbers[, "TaxonValidity"])), c(DeletionCategories, synonyms))
  
  # Stop if new resolutiosn found (need to add these to the resolution types above):
  if(length(newresolutions) > 0) stop(paste("New resolution type found!: ", newresolutions, sep = ""))
  
  # Print current processing status:
  cat("Done\nBuilding synonymy tables...")
  
  # Empty vector to store rows that correspond to some form of junior synonym:
  synonymrows <- c()
  
  # Find all junior synonym rows:
  for(i in synonyms) synonymrows <- sort(c(synonymrows, which(ResolvedTaxonNumbers[, "TaxonValidity"] == i)))
  
  # Set junior synonym matrix:
  JuniorSynonyms <- ResolvedTaxonNumbers[synonymrows, , drop = FALSE]
  
  # Create empty matrix to store senior synoyms:
  SeniorSynonyms <- matrix(nrow = 0, ncol = 8, dimnames = list(c(), c("OriginalTaxonNo", "ResolvedTaxonNo", "TaxonName", "TaxonRank", "ParentTaxonNo", "TaxonValidity", "AcceptedNumber", "AcceptedName")))
  
  # Reconcile senior synonym with database:
  currenttaxa <- PaleobiologyDBTaxaQuerier(gsub("txn:", "", ResolvedTaxonNumbers[synonymrows, "AcceptedNumber"]))
  
  # Deal with subgenera:
  currenttaxa[, "TaxonName"] <- gsub(" \\(|\\)", "", currenttaxa[, "TaxonName"])
  
  # While there are taxa with validity issues:
  while(any(!is.na(currenttaxa[, "TaxonValidity"]))) {
    
    # Isolate rows to check (i.e., just rows where validity isn't NA (i.e., resolved):
    rowstocheck <- which(!is.na(currenttaxa[, "TaxonValidity"]))
    
    # Check just those taxa:
    currenttaxa[rowstocheck, ] <- PaleobiologyDBTaxaQuerier(taxon_nos = gsub("txn:", "", currenttaxa[rowstocheck, "AcceptedNumber"]), taxon_names = currenttaxa[rowstocheck, "AcceptedName"])
    
    # Deal with subgenera:
    currenttaxa[rowstocheck, "TaxonName"] <- gsub(" \\(|\\)", "", currenttaxa[rowstocheck, "TaxonName"])
    
  }
  
  # Make current taxa into senior synonyms list:
  SeniorSynonyms <- currenttaxa
  
  # If using an Interval:
  if(!all(is.null(Interval))) {
    
    # Update resolved taxon numbers to valid taxa only:
    ResolvedTaxonNumbersInterval[which(!is.na(match(ResolvedTaxonNumbersInterval[, "InputNo"], JuniorSynonyms[, "InputNo"]))), c("OriginalTaxonNo", "ResolvedTaxonNo", "TaxonName", "TaxonRank", "ParentTaxonNo", "TaxonValidity", "AcceptedNumber", "AcceptedName")] <- SeniorSynonyms[match(ResolvedTaxonNumbersInterval[, "InputNo"], JuniorSynonyms[, "InputNo"])[!is.na(match(ResolvedTaxonNumbersInterval[, "InputNo"], JuniorSynonyms[, "InputNo"]))], c("OriginalTaxonNo", "ResolvedTaxonNo", "TaxonName", "TaxonRank", "ParentTaxonNo", "TaxonValidity", "AcceptedNumber", "AcceptedName")]
    
  }
  
  # Only complete this step if including specimen-level OTUs (there will not be any at this stage anyway if set as FALSE):
  if(IncludeSpecimenLevelOTUs) {
    
    # Print current processing status:
    cat("Done\nChecking validity of indeterminate taxon reconciliations...")
    
    # Get list of indeterminates:
    indeterminates <- TaxonomyMatrix[which((unlist(lapply(lapply(strsplit(TaxonomyMatrix[, "TaxonName"], "_"), '==', "aff"), sum)) + unlist(lapply(lapply(strsplit(TaxonomyMatrix[, "TaxonName"], "_"), '==', "cf"), sum)) + unlist(lapply(lapply(strsplit(TaxonomyMatrix[, "TaxonName"], "_"), '==', "indet"), sum)) + unlist(lapply(lapply(strsplit(TaxonomyMatrix[, "TaxonName"], "_"), '==', "sp"), sum))) > 0), "TaxonName"]
    
    # For each indeterminate:
    for(i in indeterminates) {
      
      # Get resolved row number:
      resolvedrownumber <- which(ResolvedTaxonNumbers[, "InputNo"] == TaxonomyMatrix[which(TaxonomyMatrix[, "TaxonName"] == i), "TaxonNo"])
      
      # If a possible invalid taxon (validity is not blank):
      if(!is.na(ResolvedTaxonNumbers[resolvedrownumber, "TaxonValidity"])) {
        
        # Get accepted number of taxon (may be NA):
        AcceptedNumber <- gsub("txn:", "", ResolvedTaxonNumbers[resolvedrownumber, "AcceptedNumber"])
        
        # If accepted number is blank (NA) stop adn warn user taxon is invalid:
        if(is.na(AcceptedNumber)) stop(paste(i, " assigned to a taxon that is invalid, consider renaming.", sep = ""))
        
        # If accepted numebr is different to input number stop and warn user taxon is synonymised:
        if(AcceptedNumber != ResolvedTaxonNumbers[resolvedrownumber, "InputNo"]) stop(paste(i, " assigned to a taxon that is invalid, consider renaming.", sep = ""))
        
      }
      
    }
    
  }

  # Print current processing status:
  cat("Done\nDeleting taxa resolved as nomen dubium and the like...")
  
  # Get input numbers that should be deleted:
  NumbersToDelete <- ResolvedTaxonNumbers[unlist(lapply(as.list(DeletionCategories), function(x) which(ResolvedTaxonNumbers[, "TaxonValidity"] == x))), "InputNo"]
  
  # As long as there are numbers to delete:
  if(length(NumbersToDelete) > 0) {
    
    # Remove any taxa to delete:
    MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {TaxonNumbers <- do.call(rbind, strsplit(rownames(x$Matrix), split = "%%%%"))[, 1]; DeleteRows <- sort(match(NumbersToDelete, TaxonNumbers)); if(length(DeleteRows) > 0) x$Matrix <- x$Matrix[-DeleteRows, , drop = FALSE]; x})
    
    # Prune matrices following deletion:
    MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- PisaniMRPPrune(Claddis::MakeMorphMatrix(x$Matrix, Weights = x$Weights, ignore.duplicate.taxa = TRUE)); x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
    
  }
  
  # Print current processing status:
  cat("Done\nReplacing junior synonyms with senior synonyms...")
  
  # Rebuild junior synonyms into a vector of names:
  JuniorSynonymsVector <- paste(JuniorSynonyms[, "InputNo"], gsub(" ", "_", JuniorSynonyms[, "TaxonName"]), sep = "%%%%")
  
  # Rebuild senior synonyms into a vector of names:
  SeniorSynonymsVector <- paste(unlist(lapply(apply(SeniorSynonyms[, c("OriginalTaxonNo", "ResolvedTaxonNo"), drop = FALSE], 1, as.list), function(x) unname(gsub("txn:|var:", "", unlist(x)[!is.na(unlist(x))][1])))), gsub(" ", "_", SeniorSynonyms[, "TaxonName"]), sep = "%%%%")
  
  # Build synonym matrix:
  SynonymyMatrix <- cbind(JuniorSynonymsVector, SeniorSynonymsVector)
  
  # Replace junior synonyms with senior synonyms:
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {SynonymyRows <- sort(match(SynonymyMatrix[, 1], rownames(x$Matrix))); if(length(SynonymyRows) > 0) rownames(x$Matrix)[SynonymyRows] <- SynonymyMatrix[match(rownames(x$Matrix)[SynonymyRows], SynonymyMatrix[, 1]), 2]; x})

  # Collapse any duplicate taxa created by this substitution:
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- Claddis::MakeMorphMatrix(x$Matrix, Weights = x$Weights, ignore.duplicate.taxa = TRUE); if(any(duplicated(rownames(y$Matrix_1$Matrix)))) {DuplicateNames <- rownames(y$Matrix_1$Matrix)[duplicated(rownames(y$Matrix_1$Matrix))]; y <- CollapseDuplicateTaxonMRP(y)}; x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
  
  # Prune characters made redundant by these collapses:
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- PisaniMRPPrune(Claddis::MakeMorphMatrix(x$Matrix, Weights = x$Weights, ignore.duplicate.taxa = TRUE)); x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
  
  # GOT TO HERE WITH REFACTOR (BUT HAVE JUMPED AROUND A BUNCH, SO...)

  # Print current processing status:
  cat("Done\nBuilding taxonomy...")
  
  # Get a list of the valid OTU names (may be pruned down later to just those in target clade):
  ValidOTUNames <- unique(unlist(lapply(lapply(MRPList, '[[', "Matrix"), rownames)))[grep("_", unique(unlist(lapply(lapply(MRPList, '[[', "Matrix"), rownames))))]
  
  # Replace junior with senior synonyms in resolved names matrix:
  ResolvedTaxonNumbers[synonymrows, c("OriginalTaxonNo", "ResolvedTaxonNo", "TaxonName", "TaxonRank", "ParentTaxonNo", "TaxonValidity", "AcceptedNumber", "AcceptedName")] <- SeniorSynonyms
  
  # Overwrite resolved number with resolved taxon number:
  ResolvedTaxonNumbers[, "ResolvedTaxonNo"] <- gsub("txn:|var:", "", unlist(lapply(lapply(apply(ResolvedTaxonNumbers[, c("OriginalTaxonNo", "ResolvedTaxonNo")], 1, sort), rev), '[', 1)))
  
  # Remove original taxon number column:
  ResolvedTaxonNumbers <- ResolvedTaxonNumbers[, -which(colnames(ResolvedTaxonNumbers) == "OriginalTaxonNo")]
  
  # Remove deleted taxa from resolved names matrix (if there are any):
  if(length(which(!is.na(ResolvedTaxonNumbers[, "TaxonValidity"]))) > 0) ResolvedTaxonNumbers <- ResolvedTaxonNumbers[-which(!is.na(ResolvedTaxonNumbers[, "TaxonValidity"])), , drop = FALSE]
  
  # Reformat parent taxon numbers into just numbers:
  ResolvedTaxonNumbers[, "ParentTaxonNo"] <- gsub("txn:", "", ResolvedTaxonNumbers[, "ParentTaxonNo"])
  
  # Collapse resolved matrix to just field with values (i.e., drop valid and senior synonym columns):
  ResolvedTaxonNumbers <- ResolvedTaxonNumbers[, c("ResolvedTaxonNo", "TaxonName", "TaxonRank", "ParentTaxonNo")]
  
  # If doing something with missing species (i.e., those not currently included as OTUs, but existing in target clade):
  if(MissingSpecies != "exclude") {
    
    # Find all children of target clade:
    AllChildren <- PaleobiologyDBChildFinder(taxon_nos = "1", taxon_names = TargetClade, validonly = TRUE, returnrank = "3", interval = Interval)
    
    # Deal with subgenera:
    AllChildren[, "TaxonName"] <- gsub(" \\(|\\)", "", AllChildren[, "TaxonName"])
    
    # If inserting all missing species get all possible species parent numbers:
    if(MissingSpecies == "all") CurrentSpeciesParentNumbers <- unique(c(gsub("txn:", "", AllChildren[, "ParentTaxonNo"]), ResolvedTaxonNumbers[which(ResolvedTaxonNumbers[, "TaxonRank"] == 3), "ParentTaxonNo"]))
    
    # If only inserting missing species at genus-level find parent numbers of all current species (i.e., potential genera to add):
    if(MissingSpecies == "genus") CurrentSpeciesParentNumbers <- unique(ResolvedTaxonNumbers[which(ResolvedTaxonNumbers[, "TaxonRank"] == 3), "ParentTaxonNo"])
    
    # Find any parents not already present in resolved numbers matrix:
    AsYetUnsampledSpeciesParents <- setdiff(CurrentSpeciesParentNumbers, ResolvedTaxonNumbers[, "ResolvedTaxonNo"])
    
    # If such parents exist:
    if(length(AsYetUnsampledSpeciesParents) > 0) {
      
      # Find unsampled species parents:
      CurrentSpeciesParents <- PaleobiologyDBTaxaQuerier(taxon_nos = AsYetUnsampledSpeciesParents)
      
      # Deal with subgenera:
      CurrentSpeciesParents[, "TaxonName"] <- gsub(" \\(|\\)", "", CurrentSpeciesParents[, "TaxonName"])
      
      # Find rows corresponding to valid genera:
      ValidGenusRows <- intersect(which(is.na(CurrentSpeciesParents[, "TaxonValidity"])), which(CurrentSpeciesParents[, "TaxonRank"] == "5"))
      
      # If there are valid genera then add these to resolved taxon numbers:
      if(length(ValidGenusRows) > 0) ResolvedTaxonNumbers <- rbind(ResolvedTaxonNumbers, cbind(unname(gsub("txn:|var:", "", unlist(lapply(lapply(lapply(apply(CurrentSpeciesParents[ValidGenusRows, c("OriginalTaxonNo", "ResolvedTaxonNo"), drop = FALSE], 1, list), unlist), sort, decreasing = TRUE), '[', 1)))), CurrentSpeciesParents[ValidGenusRows, c("TaxonName", "TaxonRank")] , gsub("txn:", "", CurrentSpeciesParents[ValidGenusRows, "ParentTaxonNo"])))
      
    }
    
    # If including all species:
    if(MissingSpecies == "all") {
      
      # Update Valid OTUs accordingly
      ValidOTUNames <- unique(c(ValidOTUNames, gsub(" ", "_", paste(unlist(lapply(strsplit(gsub("NA", "", paste(AllChildren[, "OriginalTaxonNo"], AllChildren[, "ResolvedTaxonNo"], sep = "")), split = "var:|txn:"), '[[', 2)), AllChildren[, "TaxonName"], sep = "%%%%"))))
      
      # Set new children to add to resolved taxon numbers later:
      NewChildren <- matrix(c(gsub("txn:|var:", "", unlist(lapply(apply(AllChildren[, c("OriginalTaxonNo", "ResolvedTaxonNo")], 1, sort, decreasing = TRUE), '[', 1))), AllChildren[, c("TaxonName", "TaxonRank")], gsub("txn:", "", AllChildren[, "ParentTaxonNo"])), ncol = 4)
      
    }
    
    # If including only species assigned to genus-level OTUs:
    if(MissingSpecies == "genus") {
      
      # Get current genus numbers (to check what has already been included):
      CurrentGenusNumbers <- ResolvedTaxonNumbers[which(ResolvedTaxonNumbers[, "TaxonRank"] == 5), "ResolvedTaxonNo"]
      
      # Get children of sampled genera:
      GeneraChildren <- PaleobiologyDBChildFinder(taxon_nos = CurrentGenusNumbers, validonly = TRUE, returnrank = "3")
      
      # Deal with subgenera:
      GeneraChildren[, "TaxonName"] <- gsub(" \\(|\\)", "", GeneraChildren[, "TaxonName"])
      
      # Update valid OTUs with children of all sampled genera:
      ValidOTUNames <- unique(c(ValidOTUNames, gsub(" ", "_", paste(unlist(lapply(strsplit(gsub("NA", "", paste(GeneraChildren[, "OriginalTaxonNo"], GeneraChildren[, "ResolvedTaxonNo"], sep = "")), split = "var:|txn:"), '[[', 2)), GeneraChildren[, "TaxonName"], sep = "%%%%"))))
      
      # Set new children to add to resolved taxon numbers later:
      NewChildren <- matrix(c(gsub("txn:|var:", "", unlist(lapply(apply(GeneraChildren[, c("OriginalTaxonNo", "ResolvedTaxonNo")], 1, sort, decreasing = TRUE), '[', 1))), GeneraChildren[, c("TaxonName", "TaxonRank")], gsub("txn:", "", GeneraChildren[, "ParentTaxonNo"])), ncol = 4)
      
    }
    
    # Find any new children not already included in resolved taxon numbers list:
    ChildrenToAdd <- setdiff(NewChildren[, 1], ResolvedTaxonNumbers[, "ResolvedTaxonNo"])
    
    # If there are children to add then add them to resolved taxon numbers:
    if(length(ChildrenToAdd) > 0) ResolvedTaxonNumbers <- rbind(ResolvedTaxonNumbers, NewChildren[match(ChildrenToAdd, NewChildren[, 1]), ])
    
  }
  
  # Find parents of all resolved taxon numbers (need to make sure these are valid or will propogate errors):
  ParentMatrix <- PaleobiologyDBTaxaQuerier(ResolvedTaxonNumbers[, "ParentTaxonNo"], original = TRUE)
  
  # If any are invalid then update to valid version:
  if(any(!is.na(ParentMatrix[, "AcceptedNumber"]))) ParentMatrix[!is.na(ParentMatrix[, "AcceptedNumber"]), ] <- do.call(rbind, lapply(as.list(gsub("txn:", "", ParentMatrix[!is.na(ParentMatrix[, "AcceptedNumber"]), "AcceptedNumber"])), function(x) PaleobiologyDBTaxaQuerier(x, original = FALSE)))
  
  # Update parent numbers to valid versions only:
  ResolvedTaxonNumbers[, "ParentTaxonNo"] <- unname(unlist(lapply(apply(ParentMatrix[, c("OriginalTaxonNo", "ResolvedTaxonNo")], 1, list), function(x) {x <- unlist(x); gsub("txn:|var:", "", x[!is.na(x)][1])})))
  
  # Get initial parent child relationships based on OTUs:
  ParentChildRelationships <- paste(unlist(lapply(strsplit(ValidOTUNames, "%%%%"), '[', 1)), ResolvedTaxonNumbers[match(unlist(lapply(strsplit(ValidOTUNames, "%%%%"), '[', 1)), ResolvedTaxonNumbers[, "ResolvedTaxonNo"]), "ParentTaxonNo"], sep = " belongs to ")
  
  # If including specimen level OTUs:
  if(IncludeSpecimenLevelOTUs) {
    
    # Find which rows correspond to indeterminate and sp taxa (i.e., those where parent should be initial reconciliation):
    indetsandsps <- sort(c(which(unlist(lapply(lapply(lapply(lapply(lapply(strsplit(ValidOTUNames, "%%%%"), '[', 2), strsplit, split = "_"), unlist), '==', "indet"), any))), which(unlist(lapply(lapply(lapply(lapply(lapply(strsplit(ValidOTUNames, "%%%%"), '[', 2), strsplit, split = "_"), unlist), '==', "sp"), any)))))
    
    # If such taxa exist then update parent child relationships accordingly:
    if(length(indetsandsps) > 0) ParentChildRelationships[indetsandsps] <- paste(unlist(lapply(strsplit(ValidOTUNames[indetsandsps], "%%%%"), '[', 1)), unlist(lapply(strsplit(ValidOTUNames[indetsandsps], "%%%%"), '[', 1)), sep = " belongs to ")
    
  }
  
  # Get list of new children (for which parents are needed) - excludes "Life" which has no parent:
  newchildren <- setdiff(unlist(lapply(strsplit(ParentChildRelationships, " belongs to "), '[', 2)), "28595")
  
  # As long as there are still children in need of parents:
  while(length(newchildren) > 0) {
    
    # Find any numbers missing for the taxonomy name resolution matrix:
    missingfromresolutions <- newchildren[which(is.na(match(newchildren, ResolvedTaxonNumbers[, "ResolvedTaxonNo"])))]
    
    # If there are such numbers:
    if(length(missingfromresolutions) > 0) {
      
      # Get raw query data for new names
      rawquery <- PaleobiologyDBTaxaQuerier(taxon_nos = missingfromresolutions)
      
      # Look for synonymy rows (taxon validity is a synonym type):
      NewSynonymRows <- sort(unique(unname(unlist(lapply(as.list(synonyms), function(x) which(rawquery[, "TaxonValidity"] == x))))))
      
      # If synonyms were found:
      if(length(NewSynonymRows) > 0) {
        
        # Fix any synonyms (replace junior with senior):
        rawquery[NewSynonymRows, c("ResolvedTaxonNo", "TaxonName")] <- rawquery[NewSynonymRows, c("AcceptedNumber", "AcceptedName"), drop = FALSE]
        
        # Remove now obsolete validity data:
        rawquery[NewSynonymRows, c("TaxonValidity", "AcceptedNumber", "AcceptedName")] <- NA
        
      }
      
      # Deal with subgenera:
      rawquery[, "TaxonName"] <- gsub(" \\(|\\)", "", rawquery[, "TaxonName"])
      
      # Add formatted results of query to resolved names matrix:
      ResolvedTaxonNumbers <- rbind(ResolvedTaxonNumbers, unname(cbind(gsub("txn:|var:", "", unname(unlist(lapply(lapply(lapply(apply(rawquery[, c("OriginalTaxonNo", "ResolvedTaxonNo"), drop = FALSE], 1, list), unlist), sort, decreasing = TRUE), '[', 1)))), rawquery[, c("TaxonName", "TaxonRank"), drop = FALSE], gsub("txn:", "", rawquery[, "ParentTaxonNo"]))))
      
    }
    
    # Add new parent child relationships to list:
    ParentChildRelationships <- c(ParentChildRelationships, paste(ResolvedTaxonNumbers[match(newchildren, ResolvedTaxonNumbers[, "ResolvedTaxonNo"]), "ResolvedTaxonNo"], ResolvedTaxonNumbers[match(newchildren, ResolvedTaxonNumbers[, "ResolvedTaxonNo"]), "ParentTaxonNo"], sep = " belongs to "))
    
    # Update new children:
    newchildren <- setdiff(ResolvedTaxonNumbers[match(newchildren, ResolvedTaxonNumbers[, "ResolvedTaxonNo"]), "ParentTaxonNo"], "28595")
    
  }
  
  # If Life is missing then add it at bottom:
  if(all(!ResolvedTaxonNumbers[, "ResolvedTaxonNo"] == "28595")) ResolvedTaxonNumbers <- rbind(ResolvedTaxonNumbers, c("28595", "Life", "25", NA))
  
  # Convert parent-child relationships into a matrix (columns for child and parent):
  parentchildmatrix <- matrix(unlist(strsplit(ParentChildRelationships, split = " belongs to ")), ncol = 2, byrow = TRUE, dimnames = list(c(), c("Child", "Parent")))
  
  # Update parent-child matrix with child names:
  parentchildmatrix[, "Child"] <- ResolvedTaxonNumbers[match(parentchildmatrix[, "Child"], ResolvedTaxonNumbers[, "ResolvedTaxonNo"]), "TaxonName"]
  
  # Update parent-child matrix with parent names:
  parentchildmatrix[, "Parent"] <- ResolvedTaxonNumbers[match(parentchildmatrix[, "Parent"], ResolvedTaxonNumbers[, "ResolvedTaxonNo"]), "TaxonName"]
  
  # Add valid OTU names into parent-child matrix:
  parentchildmatrix[c(1:length(ValidOTUNames)), "Child"] <- unlist(lapply(strsplit(ValidOTUNames, "%%%%"), '[', 2))
  
  # Add missing taxon ("Life") to parent-child matrix:
  parentchildmatrix[which(is.na(parentchildmatrix[, "Parent"])), "Parent"] <- "Life"
  
  # Create empty taxonomy MRP matrix:
  TaxonomyMRP <- matrix(0, nrow = length(ValidOTUNames), ncol = length(sort(unique(parentchildmatrix[, "Parent"]))), dimnames = list(sort(unlist(lapply(strsplit(ValidOTUNames, "%%%%"), '[', 2))), sort(unique(parentchildmatrix[, "Parent"]))))
  
  # Remove duplicates:
  parentchildmatrix <- matrix(unlist(strsplit(unique(paste(parentchildmatrix[, "Child"], parentchildmatrix[, "Parent"], sep = "%%%%")), "%%%%")), ncol = 2, byrow = TRUE, dimnames = list(c(), c("Child", "Parent")))
  
  # Check for duplicate names:
  if(any(duplicated(sort(unlist(lapply(strsplit(ValidOTUNames, "%%%%"), '[[', 2)))))) stop(paste("The following OTU names are duplicated in the database (check and correct): ", paste(sort(unlist(lapply(strsplit(ValidOTUNames, "%%%%"), '[[', 2)))[duplicated(sort(unlist(lapply(strsplit(ValidOTUNames, "%%%%"), '[[', 2))))], sep = ", "), sep = ""))
  
  # For each OTU (traces up through hierarchy until its presence in every higher taxon to which it belongs is assigned)::
  for(i in 1:length(ValidOTUNames)) {
    
    # Set starting current child taxon:
    currentchild <- parentchildmatrix[i, "Child"]
    
    # Set starting current parent taxon:
    currentparent <- setdiff(parentchildmatrix[which(parentchildmatrix[, "Child"] == currentchild), "Parent"], currentchild)
    
    # Record presence of child in parent in taxonomy matrix:
    TaxonomyMRP[parentchildmatrix[i, "Child"], currentparent] <- 1
    
    # As long as the parent is not "Life" (top of taxonomic hierarchy not reached):
    while(currentparent != "Life") {
      if(length(currentparent) > 1) stop("")
      # Update child with previous parent:
      currentchild <- currentparent
      
      # Update parent with new parent:
      currentparent <- setdiff(parentchildmatrix[which(parentchildmatrix[, "Child"] == currentchild), "Parent"], currentchild)
      
      # Record presence of child in parent:
      TaxonomyMRP[parentchildmatrix[i, "Child"], currentparent] <- 1
      
    }
    
  }
  
  # Print current processing status:
  cat("Done\nDealing with subgenera...")
  
  # Find any subgenera names (as supraspecific column names only):
  subgenerarows <- grep("\\(", colnames(TaxonomyMRP))
  
  # If subgenera found make these single names (i.e., removes parentheses that will screw up Newick trees later):
  if(length(subgenerarows) > 0) colnames(TaxonomyMRP)[subgenerarows] <- gsub("\\(|\\)| ", "", colnames(TaxonomyMRP)[subgenerarows])
  
  # Correct subgenera in MRP list too:
  MRPList <- lapply(MRPList, function(x) {NamesToCheck <- rownames(x$Matrix); NamesToCheck <- gsub("_\\(|\\)", "", NamesToCheck); rownames(x$Matrix) <- NamesToCheck; x})
  
  # Print current processing status:
  cat("Done\nTidying up taxonomy...")
  
  # Check target clade is actually found:
  if(length(which(colnames(TaxonomyMRP) == TargetClade)) == 0) stop("Target clade not found in taxonomy. Check spelling/Paleobiology database validity.")
  
  # Work out which taxa are actually valid OTUs (belong to target clade):
  NewValidOTUs <- names(which(TaxonomyMRP[, TargetClade] == 1))
  
  # Find any subspecies names:
  SubspeciesNames <- NewValidOTUs[unlist(lapply(strsplit(NewValidOTUs, split = ""), function(x) all(c(sum(x == "_") == 2, length(grep("[:A-Z:]", x)) == 1))))]
  
  # Only continue if subspecies were found:
  if(length(SubspeciesNames) > 0) {
    
    # Find any subspecues where species is not in sample:
    SubspeciesWhereSpeciesIsNotFound <- SubspeciesNames[!unlist(lapply(as.list(SubspeciesNames), function(x) any(NewValidOTUs == paste(strsplit(x, split = "_")[[1]][1:2], collapse = "_"))))]
    
    # If subspecies without sampled species :
    if(length(SubspeciesWhereSpeciesIsNotFound) > 0) {
      
      # Rename these with species names:
      NewValidOTUs[match(SubspeciesWhereSpeciesIsNotFound, NewValidOTUs)] <- unlist(lapply(strsplit(SubspeciesWhereSpeciesIsNotFound, split = "_"), function(x) paste(x[1:2], collapse = "_")))
      
      # Rename taxonomy MRP rownames with species names too:
      rownames(TaxonomyMRP)[match(SubspeciesWhereSpeciesIsNotFound, rownames(TaxonomyMRP))] <- unlist(lapply(strsplit(SubspeciesWhereSpeciesIsNotFound, split = "_"), function(x) paste(x[1:2], collapse = "_")))
      
    }
    
    # Collapse subspecies back to just the names where the species is already sampled:
    SubspeciesNames <- setdiff(SubspeciesNames, SubspeciesWhereSpeciesIsNotFound)
    
    # If there are subspecies where species already exists then prune these from the taxonomy MRP:
    if(length(SubspeciesNames) > 0) NewValidOTUs <- NewValidOTUs[-match(SubspeciesNames, NewValidOTUs)]
    
  }
  
  # Modify this if using intervals:
  if(!all(is.null(Interval))) NewValidOTUs <- setdiff(NewValidOTUs, gsub(" ", "_", ResolvedTaxonNumbersInterval[, "TaxonName"]))
  
  # Can now strip out numbers from taxon names:
  for(i in 1:length(MRPList)) if(!is.null(rownames(MRPList[[i]]$Matrix))) rownames(MRPList[[i]]$Matrix) <- unlist(lapply(strsplit(rownames(MRPList[[i]]$Matrix), "%%%%"), '[', 2))
  
  # Collapse taxonomy MRP to just new valid taxa:
  TaxonomyMRP <- TaxonomyMRP[NewValidOTUs, ]
  
  # Make taxonomy MRP into list:
  TaxonomyMRPlist <- split(TaxonomyMRP, rep(1:ncol(TaxonomyMRP), each = nrow(TaxonomyMRP)))
  
  # Add column names to list:
  names(TaxonomyMRPlist) <- colnames(TaxonomyMRP)
  
  # Find higher taxa for which every taxon is present:
  redundanthighertaxa <- colnames(TaxonomyMRP)[intersect(which(unlist(lapply(lapply(TaxonomyMRPlist, unique), length)) == 1), which(unlist(lapply(lapply(TaxonomyMRPlist, unique), '[', 1)) == 1))]
  
  # Empty higher taxa:
  emptyhighertaxa <- colnames(TaxonomyMRP)[intersect(which(unlist(lapply(lapply(TaxonomyMRPlist, unique), length)) == 1), which(unlist(lapply(lapply(TaxonomyMRPlist, unique), '[', 1)) == 0))]
  
  # Find taxonomic autapomorphies (those with just one OTU and hence redundant):
  taxonomicautapomorphies <- names(which(unlist(lapply(TaxonomyMRPlist, sum)) == 1))
  
  # Collapse taxonomy MRP by removing constant characters (i.e., most of the subgroups just established - not autapomorphies as they can be substitutes later!):
  TaxonomyMRP <- TaxonomyMRP[, -match(c(redundanthighertaxa, emptyhighertaxa), colnames(TaxonomyMRP)), drop = FALSE]
  
  # Print current processing status:
  cat("Done\nSubstituting valid OTUs for supraspecific taxa...")
  
  # Find datasets with (valid) surpaspecific OTUs:
  datasetswithsupraspecificOTUs <- which(unlist(lapply(lapply(lapply(lapply(MRPList, '[[', "Matrix"), rownames), intersect, y = colnames(TaxonomyMRP)), length)) > 0)
  
  # If such data sets exist:
  if(length(datasetswithsupraspecificOTUs) > 0) {
    
    # For each such data set:
    for(i in datasetswithsupraspecificOTUs) {
      
      # Find higher taxa that will need to be replaced:
      highertaxatoreplace <- intersect(rownames(MRPList[[i]]$Matrix), colnames(TaxonomyMRP))
      
      # For each higher taxon to replace:
      for(j in highertaxatoreplace) {
        
        # Find substitue names from taxonomy:
        substitutenames <- names(which(TaxonomyMRP[, j] == 1))
        
        # Add these to end of matrix using coding for higher taxon:
        MRPList[[i]]$Matrix <- rbind(MRPList[[i]]$Matrix, matrix(rep(MRPList[[i]]$Matrix[j, ], length(substitutenames)), nrow = length(substitutenames), byrow = TRUE, dimnames = list(substitutenames, c())))
        
        # Remove now replaced higher taxon from matrix:
        MRPList[[i]]$Matrix <- MRPList[[i]]$Matrix[-which(rownames(MRPList[[i]]$Matrix) == j), , drop = FALSE]
        
      }
      
    }
    
  }
  
  # Print current processing status:
  cat("Done\nRetracting subspecies into species...")
  
  # Replace all subspecies with their species name:
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {UnderscoreAndCapitalCounts <- matrix(unlist(lapply(strsplit(rownames(x$Matrix), split = ""), function(y) c(sum(y == "_"), length(grep("[:A-Z:]", y))))), ncol = 2, byrow = TRUE, dimnames = list(c(), c("Underscores", "Capitals"))); SubspeciesRows <- intersect(which(UnderscoreAndCapitalCounts[, "Underscores"] == 2), which(UnderscoreAndCapitalCounts[, "Capitals"] == 1)); if(length(SubspeciesRows) > 0) rownames(x$Matrix)[SubspeciesRows] <- unlist(lapply(strsplit(rownames(x$Matrix)[SubspeciesRows], split = "_"), function(z) paste(z[1:2], collapse = "_"))); x})
  
  # Collapse any duplicate taxon names:
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- Claddis::MakeMorphMatrix(x$Matrix, Weights = x$Weights, ignore.duplicate.taxa = TRUE); if(any(duplicated(rownames(y$Matrix_1$Matrix)))) {DuplicateNames <- rownames(y$Matrix_1$Matrix)[duplicated(rownames(y$Matrix_1$Matrix))]; if(length(DuplicateNames) > 0) cat(paste("\nDuplicate resolved OTU name(s) found post higher-taxon substitution in ", x$FileName, ": ", paste(DuplicateNames, collapse = ", "), ". Check this is correct.", sep = "")); y <- CollapseDuplicateTaxonMRP(y)}; x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
  
  # Print current processing status:
  cat("Done\nFurther tidying of taxonomy...")
  
  # If applying an Interval:
  if(!all(is.null(Interval))) {
    
    # Find any rows to delete (because they represent taxa outside the Interval):
    RowsToDelete <- sort(match(gsub(" ", "_", ResolvedTaxonNumbersInterval[ResolvedTaxonNumbersInterval[, "TaxonRank"] == "3", "TaxonName"]), rownames(TaxonomyMRP)))
    
    # If found then remove these from the taxonomy MRP:
    if(length(RowsToDelete) > 0) TaxonomyMRP <- TaxonomyMRP[-RowsToDelete, , drop = FALSE]
    
  }
  
  # Collapse taxonomy MRP by removing autapomorphic characters (if any):
  if(length(taxonomicautapomorphies) > 0) TaxonomyMRP <- TaxonomyMRP[, -match(taxonomicautapomorphies, colnames(TaxonomyMRP)), drop = FALSE]
  
  # Overwrite taxonomy MRP list with strings for each column:
  TaxonomyMRPlist <- apply(TaxonomyMRP, 2, paste, collapse = "")
  
  # Find any duplicated MRP strings:
  duplicatedMRPstrings <- rle(sort(TaxonomyMRPlist))$values[which(rle(sort(TaxonomyMRPlist))$lengths > 1)]
  
  # If there are duplicated columns (i.e., redundant MRP characters in the taxonomy):
  if(length(duplicatedMRPstrings) > 0) {
    
    # For each duplicated character:
    for(i in duplicatedMRPstrings) {
      
      # Get duplicated columns for current MRP string:
      duplicatedcolumns <- which(TaxonomyMRPlist == i)
      
      # Form new column name by collapsing higher taxa that are duplicated to a single string:
      newcolumnname <- paste(names(duplicatedcolumns), collapse = "_et_")
      
      # Overwrite first duplicated column name with new collapsed name:
      colnames(TaxonomyMRP)[duplicatedcolumns[1]] <- newcolumnname
      
      # Remove redundant columns from matrix:
      TaxonomyMRP <- TaxonomyMRP[, -duplicatedcolumns[2:length(duplicatedcolumns)], drop = FALSE]
      
      # Remove redundant columns from list:
      TaxonomyMRPlist <- TaxonomyMRPlist[-duplicatedcolumns[2:length(duplicatedcolumns)]]
      
    }
    
  }
  
  # Print current processing status:
  cat("Done\nRemoving outgroups, empty supraspecific taxa and those outside the sampling interval...")
  
  # Find any remaining taxa that now need to be deleted (outgroups to target clade and empty higher taxa):
  TaxaToDelete <- setdiff(unlist(lapply(lapply(MRPList, '[[', "Matrix"), rownames)), NewValidOTUs)
  
  # If applying an Interval then add taxa outside of it to the deletes list:
  if(!all(is.null(Interval))) TaxaToDelete <- unique(c(TaxaToDelete, gsub(" ", "_", ResolvedTaxonNumbersInterval[ResolvedTaxonNumbersInterval[, "TaxonRank"] == "3", "TaxonName"])))
  
  # If there are species to exclude:
  if(length(SpeciesToExclude) > 0) {
    
    # Build vector of all current OTU names:
    OTUNames <- unique(unlist(lapply(MRPList, function(x) rownames(x$Matrix))))
    
    # Find any missing names (in exclude list but not in tree):
    MissingNames <- setdiff(SpeciesToExclude, OTUNames)
    
    # If any are found stop and warn user:
    if(length(MissingNames) > 0) stop(paste("The following SpeciesToExclude were not actually found in the data: ", paste(MissingNames, collapse = ", "), ". Check they are spelled correctly and try again.", sep = ""))
    
    # Add species to exclude to taxa to delete:
    TaxaToDelete <- unique(c(TaxaToDelete, SpeciesToExclude))
    
    # Remove species to exclude from Taxonomy MRP (as lomg as they are still there):
    if(length(intersect(SpeciesToExclude, rownames(TaxonomyMRP))) > 0) TaxonomyMRP <- TaxonomyMRP[-match(intersect(SpeciesToExclude, rownames(TaxonomyMRP)), rownames(TaxonomyMRP)), , drop = FALSE]
    
    # Find any columns to delete (duplicated, autapomorphic or constant):
    ColumnsToDalete <- unique(c(which(duplicated(apply(TaxonomyMRP, 2, paste, collapse = ""))), unname(which(apply(TaxonomyMRP, 2, sum) < 2))))
    
    # If columns are to be deleted then delete them:
    if(length(ColumnsToDalete) > 0) TaxonomyMRP <- TaxonomyMRP[, -ColumnsToDalete]
  
  }
  
  # Delete taxa from every matrix they occur in:
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {DeleteRows <- match(intersect(TaxaToDelete, rownames(x$Matrix)), rownames(x$Matrix)); if(length(DeleteRows) > 0) x$Matrix <- x$Matrix[-DeleteRows, , drop = FALSE]; x})
  
  # Prune redundant characters from matrices following taxon deletion(s):
  MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- PisaniMRPPrune(Claddis::MakeMorphMatrix(x$Matrix, Weights = x$Weights, ignore.duplicate.taxa = TRUE)); x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
  
  # Print current processing status:
  cat("Done\nProducing taxonomy tree...")
  
  # Duplicated Taxonomy MRP to create a collapsable version for generating taxonomy Newick string:
  TaxonomyMRPNewick <- TaxonomyMRP
  
  # Get order of columns to collapse to form MRP
  columncollapseorder <- order(apply(TaxonomyMRPNewick, 2, sum))
  
  # For each column ("clade") in order from smallest to largest:
  for(i in columncollapseorder) {
    
    # Get taxa present in current clade:
    taxonrows <- which(TaxonomyMRPNewick[, i] == 1)
    
    # Create new partial Newick string for current clade (node):
    newNewickstring <- paste("(", paste(names(taxonrows), collapse = ","), ")", colnames(TaxonomyMRPNewick)[i], sep = "")
    
    # Replace row name with new Newick string:
    rownames(TaxonomyMRPNewick)[taxonrows[1]] <- newNewickstring
    
    # Remove now redundant taxa from Newick matrix:
    TaxonomyMRPNewick <- TaxonomyMRPNewick[-taxonrows[2:length(taxonrows)], , drop = FALSE]
    
  }
  
  # Complete Newick string and over write Taxonomy Newick matrix:
  TaxonomyMRPNewick <- paste("(", paste(rownames(TaxonomyMRPNewick), collapse = ","), ")", TargetClade, ";", sep = "")
  
  # Ladderize taxonomy tree for neatness!:
  TaxonomyMRPTree <- ladderize(read.tree(text = TaxonomyMRPNewick))

  # If there are higher taxa to collapse:
  if(length(HigherTaxaToCollapse) > 0) {
    
    # Print current processing status:
    cat("Done\nCollapsing higher taxa...")
    
    # Find any higher taxa actually present in target clade:
    HigherTaxaInTargetClade <- colnames(TaxonomyMRP)
    
    # Find any missing names (in collapse list but not in taxonomy):
    MissingHigherTaxa <- setdiff(HigherTaxaToCollapse, HigherTaxaInTargetClade)
    
    # If any are found stop and warn user:
    if(length(MissingHigherTaxa) > 0) stop(paste("The following HigherTaxaToCollapse were not actually found in the data: ", paste(MissingHigherTaxa, collapse = ", "), ". Check they are spelled correctly and are valid (according to the Paleobiology Database) and try again.", sep = ""))
    
    # Check the clades are all unique (not internested) and if not then stop and warn user:
    if(any(duplicated(unlist(lapply(as.list(HigherTaxaToCollapse), function(x) rownames(TaxonomyMRP[TaxonomyMRP[, colnames(TaxonomyMRP) == x] == 1, ])))))) stop("HigherTaxaToCollapse contains clades that are internested. Remove the internested clades and try again.")
    
    # Build taxa to rename matrix:
    TaxaToRenameMatrix <- do.call(rbind, lapply(as.list(HigherTaxaToCollapse), function(x) unname(cbind(x, rownames(TaxonomyMRP[TaxonomyMRP[, colnames(TaxonomyMRP) == x] == 1, ])))))
    
    # Build list of each clade's species compliment:
    CladeContentsList <- lapply(as.list(TaxonomyMRPTree$node.label), function(x) {NodeNumber <- which(TaxonomyMRPTree$node.label == x) + Ntip(TaxonomyMRPTree); TaxonomyMRPTree$tip.label[FindDescendants(n = NodeNumber, tree = TaxonomyMRPTree)]})
    
    # Add names of clades:
    names(CladeContentsList) <- TaxonomyMRPTree$node.label
    
    # Find any subsumed clades (to be removed from taxonomy MRP):
    SubsumedClades <- unlist(lapply(as.list(HigherTaxaToCollapse), function(x) {CurrentClade <- which(names(CladeContentsList) == x); TempCladeContents <- CladeContentsList[-CurrentClade]; names(which(unlist(lapply(TempCladeContents, function(x) length(setdiff(x, CladeContentsList[[CurrentClade]])))) == 0))}))
    
    # Build block to add to taxonomy MRP out of first taxon inside each clade to collapse:
    BlockToAddToTaxonomyMRP <- do.call(rbind, lapply(as.list(HigherTaxaToCollapse), function(x) TaxonomyMRP[TaxaToRenameMatrix[which(TaxaToRenameMatrix[, 1] == x)[1], 2], ]))
    
    # Add uppercase rownames to new block:
    rownames(BlockToAddToTaxonomyMRP) <- toupper(HigherTaxaToCollapse)
    
    # Add to taxonomy MRP:
    TaxonomyMRP <- rbind(TaxonomyMRP, BlockToAddToTaxonomyMRP)
    
    # Collapse taxonomy MRP down by removing clades and the species from the collapsed clades:
    TaxonomyMRP <- TaxonomyMRP[-unlist(lapply(as.list(TaxaToRenameMatrix[, 2]), function(x) which(rownames(TaxonomyMRP) == x))), -unlist(lapply(as.list(c(HigherTaxaToCollapse, SubsumedClades)), function(x) which(colnames(TaxonomyMRP) == x))), drop = FALSE]
    
    # Remove all but one collapsed taxa from each clade from the tree:
    TaxonomyMRPTree <- drop.tip(TaxonomyMRPTree, TaxaToRenameMatrix[-unlist(lapply(as.list(HigherTaxaToCollapse), function(x) which(TaxaToRenameMatrix[, 1] == x)[1])), 2])
    
    # Ladderize taxonomy tree for neatness!:
    TaxonomyMRPTree <- ladderize(TaxonomyMRPTree)
    
    # Build tips to replace matrix:
    TipsToReplaceMatrix <- do.call(rbind, lapply(as.list(HigherTaxaToCollapse), function(x) TaxaToRenameMatrix[which(TaxaToRenameMatrix[, 1] == x)[1], ]))
    
    # Replace tip names in tree:
    TaxonomyMRPTree$tip.label[unlist(lapply(apply(TipsToReplaceMatrix, 1, as.list), function(x) {x <- unlist(x); which(TaxonomyMRPTree$tip.label == x[2])}))] <- toupper(TipsToReplaceMatrix[, 1])
    
    # Replace names in MRP matrices:
    MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {CurrentRownames <- rownames(x$Matrix); NamesToReplace <- intersect(CurrentRownames, TaxaToRenameMatrix[, 2]); if(length(NamesToReplace) > 0) rownames(x$Matrix)[match(NamesToReplace, rownames(x$Matrix))] <- toupper(TaxaToRenameMatrix[match(NamesToReplace, TaxaToRenameMatrix[, 2]), 1]); x})
    
    # Collapse any duplicate taxa created by this substitution (very likely!):
    MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- Claddis::MakeMorphMatrix(x$Matrix, Weights = x$Weights, ignore.duplicate.taxa = TRUE); if(any(duplicated(rownames(y$Matrix_1$Matrix)))) {DuplicateNames <- rownames(y$Matrix_1$Matrix)[duplicated(rownames(y$Matrix_1$Matrix))]; y <- CollapseDuplicateTaxonMRP(y)}; x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
    
    # Prune constant characters and collapse duplicated characters:
    MRPList[ActiveMRP(MRPList)] <- lapply(MRPList[ActiveMRP(MRPList)], function(x) {y <- PisaniMRPPrune(Claddis::MakeMorphMatrix(x$Matrix, Weights = x$Weights, ignore.duplicate.taxa = TRUE)); x$Matrix <- y$Matrix_1$Matrix; x$Weights <- y$Matrix_1$Weights; x})
    
    # Update new valid OTUs:
    NewValidOTUs <- sort(rownames(TaxonomyMRP))
    
  }
  
  # If specimen-level OTUs are included:
  if(IncludeSpecimenLevelOTUs) {
    
    # Print current processing status:
    cat("Done\nAdding NAs for indet. and sp. subclades to taxonomy MRP...")
    
    # Find indeterminates and sps:
    indetsandsps <- NewValidOTUs[which((unlist(lapply(lapply(strsplit(NewValidOTUs, "_"), '==', "indet"), sum)) + unlist(lapply(lapply(strsplit(NewValidOTUs, "_"), '==', "sp"), sum))) > 0)]
    
    # If there are indeterminates and/or sps:
    if(length(indetsandsps) > 0) {
      
      # For each such taxon:
      for(i in indetsandsps) {
        
        # Find higher taxon to which it belongs:
        highertaxon <- colnames(TaxonomyMRP)[which(unlist(lapply(lapply(strsplit(colnames(TaxonomyMRP), "_et_"), '==', strsplit(i, "_")[[1]][1]), sum)) == 1)]
        
        # Find any sub (suprapseicifc taxa) for that higher taxon:
        subtaxa <- colnames(TaxonomyMRP)[which(unlist(lapply(lapply(lapply(split(TaxonomyMRP[names(which(TaxonomyMRP[, highertaxon] == 1)), , drop = FALSE], rep(1:ncol(TaxonomyMRP[names(which(TaxonomyMRP[, highertaxon] == 1)), , drop = FALSE]), each = nrow(TaxonomyMRP[names(which(TaxonomyMRP[, highertaxon] == 1)), , drop = FALSE]))), sort), unique), length)) == 2)]
        
        # If these exist then set ith taxon as being NA with respect to belonging to the subtax(a):
        if(length(subtaxa) > 1) TaxonomyMRP[i, subtaxa] <- NA
        
      }
      
    }
    
  }
  
  # For empty data sets make sure matrix is zero-by-zero and weights have no lengths:
  MRPList <- lapply(MRPList, function(x) {MatrixSize <- nrow(x$Matrix) * ncol(x$Matrix); if(MatrixSize == 0) {x$Matrix <- matrix(nrow = 0, ncol = 0); x$Weights <- vector(mode = "numeric")}; x})
  
  # Print current processing status:
  cat("Done\nGetting weighting data (publication year and dependencies)...")
  
  # Set current year (used multiple times later):
  CurrentYear <- strsplit(as.character(Sys.Date()), "-")[[1]][1]
  
  # Add publication year to each data set (in presses are ascribed the current year):
  MRPList <- lapply(MRPList, function(x) {x$PublicationYear <- gsub("[:A-Z:a-z:]|_|-", "", gsub("inpress", CurrentYear, x$FileName)); x})
  
  # Find any missing parents:
  MissingParents <- setdiff(unique(unname(unlist(lapply(MRPList, function(x) x$Parent[nchar(x$Parent) > 0])))), names(MRPList))
  
  # Find all parent data set names:
  ParentDataSets <- sort(unique(unname(unlist(lapply(MRPList, function(x) x$Parent[nchar(x$Parent) > 0])))))
  
  # Get child data sets for each parent:
  ChildDataSets <- lapply(as.list(ParentDataSets), function(x) unname(unlist(mapply(function(x, y) y$FileName[y$Parent == x], x = x, y = MRPList))))
  
  # Add names to child data sets:
  names(ChildDataSets) <- ParentDataSets
  
  # Now include any grandchildren, greatgrandchildren etc.:
  ChildDataSets <- lapply(ChildDataSets, function(x) sort(unique(c(x, unname(unlist(ChildDataSets[intersect(x, names(ChildDataSets))]))))))
  
  # Add sibling relationships to data sets with shared parents and update parents field with grandparents, greatgrandparents etc.:
  MRPList <- lapply(MRPList, function(x) {SiblingVector <- c(x$Sibling, setdiff(ChildDataSets[[match(x$Parent, ParentDataSets)]], x$FileName)); if(any(nchar(SiblingVector)) > 0) SiblingVector <- SiblingVector[nchar(SiblingVector) > 0]; x$Sibling <- unique(SiblingVector); if(!is.null(MRPList[[x$Parent[1]]]$Parent)) while(nchar(MRPList[[x$Parent[1]]]$Parent) > 0) x$Parent <- c(MRPList[[x$Parent[1]]]$Parent, x$Parent); x})
  
  # OLD LINE FOR ABOVE THAT I AM PRETTY SURE IS BROKEN BUT AM LEAVING HERE FOR NOW IN CASE IT AIN'T
  # Add sibling relationships to data sets with shared parents and update parents field with grandparents, greatgrandparents etc.:
  #MRPList <- lapply(MRPList, function(x) {SiblingVector <- c(x$Sibling, setdiff(ChildDataSets[[match(x$Parent, ParentDataSets)]], x$FileName)); if(any(nchar(SiblingVector)) > 0) SiblingVector <- SiblingVector[nchar(SiblingVector) > 0]; x$Sibling <- unique(SiblingVector); x$Parent <- names(which(unlist(lapply(ChildDataSets, function(y) length(intersect(y, x$FileName)))) > 0)); x})
  
  # Build an empty redundant parent list for later use if no parents exist:
  RedundantParents <- vector(mode = "character")
  
  # If parents exist build a vector of those that are redundant (at least one child data set contains their full taxonomic complement):
  if(length(unlist(lapply(MRPList[ActiveMRP(MRPList)], function(x) x$Parent[nchar(x$Parent) > 0]))) > 0) RedundantParents <- unlist(lapply(as.list(unique(unlist(lapply(MRPList[ActiveMRP(MRPList)], function(x) x$Parent[nchar(x$Parent) > 0])))), function(x) {Children <- unname(unlist(lapply(MRPList[ActiveMRP(MRPList)], function(y) if(any(y$Parent == x)) y$FileName))); if(any(unlist(lapply(as.list(Children), function(y) length(setdiff(rownames(MRPList[[x]]$Matrix), rownames(MRPList[[y]]$Matrix))))) == 0)) x}))
  
  # If redundant parents were found then collapse these data sets back to empty matrix and weights::
  if(length(RedundantParents) > 0) MRPList[RedundantParents] <- lapply(MRPList[RedundantParents], function(x) {x$Matrix <- matrix(ncol = 0, nrow = 0); x$Weights <- vector(mode = "numeric"); x})
  
  # Find any remaining active parent data sets:
  ActiveParents <- names(which(unlist(lapply(MRPList[setdiff(ParentDataSets, MissingParents)], function(x) nrow(x$Matrix) * ncol(x$Matrix))) > 0))
  
  # If active parents remain:
  if(length(ActiveParents) > 0) {
    
    # Add children of active parent to siblings:
    MRPList[ActiveParents] <- lapply(MRPList[ActiveParents], function(x) {SiblingVector <- unique(c(ChildDataSets[[x$FileName]], x$Sibling)); x$Sibling <- unique(SiblingVector[nchar(SiblingVector) > 0]); x})
    
    # Add parent as sibling of offspring:
    MRPList <- lapply(MRPList, function(x) {SiblingVector <- c(x$Sibling, intersect(x$Parent, ActiveParents)); if(any(nchar(SiblingVector) > 0)) SiblingVector <- SiblingVector[nchar(SiblingVector) > 0]; x$Sibling <- sort(unique(SiblingVector)); x})
    
  }
  
  # Find any empty data sets to remove:
  RemovedSourceData <- sort(names(which(unlist(lapply(MRPList, function(x) nrow(x$Matrix) * ncol(x$Matrix))) == 0)))
  
  # Remove data sets from MRPList:
  MRPList[RemovedSourceData] <- NULL
  
  # Remove any dead siblings:
  MRPList <- lapply(MRPList, function(x) {x$Sibling <- intersect(x$Sibling, names(MRPList)); x})
  
  # If usinga  veil line:
  if(VeilLine) {
    
    # Print current processing status:
    cat("Done\nApplying veil line...")
    
    # Start with current year as veil year:
    CurrentVeilYear <- as.numeric(CurrentYear)
    
    # Set current taxa included as being from current veil year (to present):
    CurrentTaxaIncluded <- unique(unlist(lapply(MRPList[as.numeric(unlist(lapply(MRPList, function(x) x$PublicationYear))) >= CurrentVeilYear], function(x) rownames(x$Matrix))))
    
    # Make a stop point (where all taxa are sampled):
    StopPoint <- length(unique(unlist(lapply(MRPList, function(x) rownames(x$Matrix)))))
    
    # While not all taxa are included in current sample:
    while(length(CurrentTaxaIncluded) < StopPoint) {
      
      # Increment one year back in time:
      CurrentVeilYear <- CurrentVeilYear - 1

      # Update current taxa included:
      CurrentTaxaIncluded <- unique(unlist(lapply(MRPList[as.numeric(unlist(lapply(MRPList, function(x) x$PublicationYear))) >= CurrentVeilYear], function(x) rownames(x$Matrix))))
      
    }
    
    # Find any data sets to remove (older than veil year):
    DataSetsToRemove <- which(as.numeric(unlist(lapply(MRPList, function(x) x$PublicationYear))) < CurrentVeilYear)
    
    # If data sets to remove:
    if(length(DataSetsToRemove) > 0) {
      
      # Add to removed source data vector:
      RemovedSourceData <- sort(c(RemovedSourceData, names(MRPList)[DataSetsToRemove]))
      
      # Remove from MRP list:
      MRPList <- MRPList[-DataSetsToRemove]
      
      # Remove any new dead siblings:
      MRPList <- lapply(MRPList, function(x) {x$Sibling <- intersect(x$Sibling, names(MRPList)); x})
    
    }

  # If not using veil line:
  } else {
    
    # Set current veil year as oldest data set:
    CurrentVeilYear <- min(as.numeric(unlist(lapply(MRPList, function(x) x$PublicationYear))))
    
  }
  
  # Print current processing status:
  cat("Done\nCalculating weights...")
  
  # Update current year to youngest data set (in case this is not the actual current year) as will screw uo weights otherwise:
  CurrentYear <- as.character(max(as.numeric(unlist(lapply(MRPList, function(x) x$PublicationYear)))))
  
  # Reformat weights to be input weights, publication year weights (equation 1 in Supplementary Information of Lloyd et al. 2016),
  # data set dependence weights (1 / (N siblings + 1)), and clade contradiction weights (1 / frequency of contradictory clades).
  # All weights are set on a zero to one scale initially and then multiplied by RelativeWeights:
  MRPList <- lapply(MRPList, function(x) {InputWeights <- x$Weights; x$Weights <- NULL; x$InputWeights <- RelativeWeights[1] * InputWeights; x$PublicationYearWeights <- RelativeWeights[2] * (rep(((2 ^ (0.5 * (as.numeric(x$PublicationYear) - CurrentVeilYear + 1))) - 1) / ((2 ^ (0.5 * (as.numeric(CurrentYear) - CurrentVeilYear + 1))) - 1), length(InputWeights))); x$DataSetDependenceWeights <- RelativeWeights[3] * rep(1 / (length(x$Sibling) + 1), length(InputWeights)); x$CladeContradictionWeights <- RelativeWeights[4] * MRPIntraMatrixWeights(x$Matrix); x})
  
  # If using sum to combine weights collapse weights to just their sum:
  if(WeightCombination == "sum") MRPList <- lapply(MRPList, function(x) {x$Weights <- apply(rbind(x$InputWeights, x$PublicationYearWeights, x$DataSetDependenceWeights, x$CladeContradictionWeights), 2, sum); x$InputWeights <- NULL; x$PublicationYearWeights <- NULL; x$DataSetDependenceWeights <- NULL; x$CladeContradictionWeights <- NULL; x})
  
  # If using product to combine weights collapse weights to just their product (excluding zeroes) and remove other weights from list:
  if(WeightCombination == "product") MRPList <- lapply(MRPList, function(x) {x$Weights <- apply(rbind(x$InputWeights, x$PublicationYearWeights, x$DataSetDependenceWeights, x$CladeContradictionWeights), 2, function(y) {z <- as.character(y); z[z == "0"] <- NA; prod(as.numeric(z), na.rm = TRUE)}); x$InputWeights <- NULL; x$PublicationYearWeights <- NULL; x$DataSetDependenceWeights <- NULL; x$CladeContradictionWeights <- NULL; x})
  
  # Get current maximum weight:
  MaximumWeight <- max(unlist(lapply(MRPList, function(x) x$Weights)))
  
  # Get current minimum weight:
  MinimumWeight <- min(unlist(lapply(MRPList, function(x) x$Weights)))
  
  # Calculate the multiplication factor for weight rescaling (10 to 1000), but use minimum weight if there is no variance:
  MultiplicationFactor <- ifelse(1 / ((MaximumWeight - MinimumWeight) / 990) == Inf, MinimumWeight, 1 / ((MaximumWeight - MinimumWeight) / 990))
  
  # Calculate the addition factor for weight rescaling (10 to 1000):
  AdditionFactor <- 10 - (MinimumWeight * MultiplicationFactor)
  
  # Rescale weights (10 to 1000) and round results to two decimal places (best TNT can cope with):
  MRPList <- lapply(MRPList, function(x) {x$Weights <- round((x$Weights * MultiplicationFactor) + AdditionFactor, 2); x})
  
  # Print current processing status:
  cat("Done\nChecking for phylogeny-taxonomy contradictions...")
  
  # Do first pass to find any contradictions between taxonomy MRP and the fully reconciled source matrix:
  MRPList <- lapply(MRPList, function(x) {TaxonomyMRPStrings <- TaxonomyMRP[rownames(x$Matrix), ]; TaxonomyMRPStrings[is.na(TaxonomyMRPStrings)] <- "0"; TaxonomyMRPStrings <- TaxonomyMRPStrings[, apply(TaxonomyMRPStrings, 2, function(y) length(unique(y))) == 2, drop = FALSE]; x$TaxonomyContradictions <- names(which(unlist(lapply(apply(TaxonomyMRPStrings, 2, list), function(z) length(MRPCharacterContradiction(unlist(z), x$Matrix))) > 0))); x$TaxonomyContradictionProportion <- length(x$TaxonomyContradictions) / ncol(TaxonomyMRPStrings); x})
  
  # Store monophyletic taxa (those not contradicted by any phylogenetic characters - useful for chunking larger data if found):
  MonophyleticTaxa <- setdiff(colnames(TaxonomyMRP), unique(unlist(lapply(MRPList, function(x) x$TaxonomyContradictions))))
  
  # If reporting contradiction issues to the screen:
  if(ReportContradictionsToScreen) {
    
    # Find any data sets with taxonomy-phylogeny contradictions:
    ContradictionIssueDataSets <- names(unlist(lapply(MRPList, function(x) length(x$TaxonomyContradictions))) > 0)
    
    # If there are data sets with contradictions:
    if(length(ContradictionIssueDataSets) > 0) {
      
      # Build vector of contradiction warnings:
      ContradictionWarnings <- unname(unlist(lapply(MRPList[ContradictionIssueDataSets], function(x) {TaxonomyMRPSubset <- TaxonomyMRP[rownames(x$Matrix), x$TaxonomyContradictions, drop = FALSE]; if(any(is.na(TaxonomyMRPSubset))) TaxonomyMRPSubset[is.na(TaxonomyMRPSubset)] <- "0"; ListContradictions(TaxonomyMRP = TaxonomyMRPSubset, MRPMatrix = x$Matrix, ContradictionTaxa = x$TaxonomyContradictions, DataSetName = x$FileName)})))
      
      # Print warnings to screen:
      cat(ContradictionWarnings)
      
      # NEED TO BREAK THIS DOWN FURTHER AS CLEARLY SOME REDUNDANCY! (E.G. GROUPING HIGHER TAXA WITH SAME ISSUE, OR DATA SETS WITH SAME ISSUE)
      
    }
    
  }
  
  # If using a constraint:
  if(ConstraintInUse) {
    
    # Print current processing status:
    cat("Done\nApplying constraint tree...")
    
    # If a monophyly constraint add all other taxa outside the constraint (makes NAs zeroes):
    if(ConstraintType == "monophyly") MRPList[[grep(ConstraintDataSet, names(MRPList))]]$Matrix <- rbind(MRPList[[grep(ConstraintDataSet, names(MRPList))]]$Matrix, matrix("0", ncol = ncol(MRPList[[grep("Constraint", names(MRPList))]]$Matrix), nrow = length(setdiff(NewValidOTUs, rownames(MRPList[[grep(ConstraintDataSet, names(MRPList))]]$Matrix))), dimnames = list(setdiff(NewValidOTUs, rownames(MRPList[[grep(ConstraintDataSet, names(MRPList))]]$Matrix)), c())))
    
    # Build vector of atxa in constraint:
    TaxaInConstraint <- rownames(MRPList[[grep(ConstraintDataSet, names(MRPList))]]$Matrix)
    
    # Get combined weight of all non-constraint that contradicts constraint (need to know to correctly weight the constraint data):
    NonConstraintWeightsTotal <- sum(unlist(lapply(MRPList[-grep(ConstraintDataSet, names(MRPList))], function(x) {TaxaInBoth <- intersect(rownames(x$Matrix), TaxaInConstraint); if(length(TaxaInBoth) > 2) {ConstraintMRPStrings <- x$Matrix[TaxaInBoth, ]; ConstraintMatrix <- MRPList[[grep(ConstraintDataSet, names(MRPList))]]$Matrix[TaxaInBoth, ]; ConstraintMatrix[is.na(ConstraintMatrix)] <- "0"; ConstraintMatrix <- ConstraintMatrix[, apply(ConstraintMatrix, 2, function(y) length(unique(y))) == 2, drop = FALSE]; if(length(unique(as.vector(ConstraintMRPStrings))) > 1) {x$ConstraintContradictions <- x$Weights[unique(unlist(lapply(apply(ConstraintMatrix, 2, list), function(z) {MRPCharacterContradiction(unlist(z), ConstraintMRPStrings)})))]} else {x$ConstraintContradictions <- integer(0)}} else {x$ConstraintContradictions <- integer(0)}; x$ConstraintContradictions})))
    
    # If NonConstraintWeightsTotal is less than 10000 then set it at 10000 to ensure it is upweighted:
    if(NonConstraintWeightsTotal < 10000) NonConstraintWeightsTotal <- 10000
    
    # Update weights of constraint tree to maximum (allowing for weights to fall in the 999-1000 range if they represent conflicting clades):
    MRPList[[ConstraintDataSet]]$Weights <- round(MRPIntraMatrixWeights(MRPList[[ConstraintDataSet]]$Matrix) + 999, 2)
    
    # Get order of magnitude of current character count of constraint tree (need to check this won't be too big):
    OrderOfMagnitudeOfCurrentCharacterCountOfConstraint <- nchar(as.character(ceiling(NonConstraintWeightsTotal / 1000) * ncol(MRPList[[ConstraintDataSet]]$Matrix)))
    
    # If this order exceeds 10^6 (too big for memory):
    if(OrderOfMagnitudeOfCurrentCharacterCountOfConstraint > 6) {
      
      # Get order of magnitude to reduce weights by:
      OrderOfMagnitudeToReduceWeightsBy <- (OrderOfMagnitudeOfCurrentCharacterCountOfConstraint - 6) * 10
      
      # Calculate the multiplication factor for weight rescaling (10 to 1000):
      MultiplicationFactor <- 1 / (990 / ((1000 / OrderOfMagnitudeToReduceWeightsBy) - 10))
      
      # Calculate the addition factor for weight rescaling (10 to 1000):
      AdditionFactor <- 10 - (10 * MultiplicationFactor)
      
      # Rescale weights (10 to 1000) and round results to two decimal places (best TNT can cope with):
      MRPList[-grep(ConstraintDataSet, names(MRPList))] <- lapply(MRPList[-grep(ConstraintDataSet, names(MRPList))], function(x) {x$Weights <- round((x$Weights * MultiplicationFactor) + AdditionFactor, 2); x})
      
      # Update NonConstraintWeightsTotal:
      NonConstraintWeightsTotal <- sum(unname(unlist(lapply(MRPList[-grep(ConstraintDataSet, names(MRPList))], function(x) x$Weights))))
      
    }
    
    # Embiggen MRP matrix so that weights are high enough to ensure constraint gets implemented:
    MRPList[[ConstraintDataSet]]$Matrix <- metatree::EmbiggenMatrix(Claddis::MakeMorphMatrix(MRPList[[ConstraintDataSet]]$Matrix, Weights = MRPList[[ConstraintDataSet]]$Weights), N = ceiling(NonConstraintWeightsTotal / 1000))$Matrix_1$Matrix
    
    # Update weights by replicating N times as with matrix embiggining:
    MRPList[[ConstraintDataSet]]$Weights <- rep(MRPList[[ConstraintDataSet]]$Weights, ceiling(NonConstraintWeightsTotal / 1000))
    
  }

  # Print current processing status:
  cat("Done\nBuilding MRP matrix...")
  
  # Add in missing taxa as NAs to every taxon:
  MRPList <- lapply(MRPList, function(x) {MissingTaxa <- setdiff(rownames(TaxonomyMRP), rownames(x$Matrix)); if(length(MissingTaxa) > 0) x$Matrix <- rbind(x$Matrix, matrix(nrow = length(MissingTaxa), ncol = ncol(x$Matrix), dimnames = list(MissingTaxa, c()))); x$Matrix <- x$Matrix[rownames(TaxonomyMRP), , drop = FALSE]; x})
  
  # Build full MRP matrix:
  FullMRPMatrix <- MakeMorphMatrix(CharacterTaxonMatrix = cbind(do.call(cbind, lapply(MRPList, function(x) x$Matrix)), TaxonomyMRP), Weights = c(unname(unlist(lapply(MRPList, function(x) x$Weights))), rep(1, ncol(TaxonomyMRP))))
  
  # Add all zero outgroup to matrix:
  FullMRPMatrix$Matrix_1$Matrix <- rbind(matrix("0", nrow = 1, ncol = ncol(FullMRPMatrix$Matrix_1$Matrix), dimnames = list("allzero", c())), FullMRPMatrix$Matrix_1$Matrix)

  # Print current processing status:
  cat("Done\nPerforming Safe Taxonomic Reduction...")
  
  # Perform STR on full matrix:
  STRData <- SafeTaxonomicReduction(FullMRPMatrix)
  
  # Create additional STR matrix from full matrix:
  STRMRPMatrix <- STRData$reduced.matrix
  
  # Print current processing status:
  cat("Done\nCompiling and returning output...")
  
  # Compile output:
  Output <- list(FullMRPMatrix = FullMRPMatrix, STRMRPMatrix = STRMRPMatrix, TaxonomyTree = TaxonomyMRPTree, MonophyleticTaxa = MonophyleticTaxa, SafelyRemovedTaxa = STRData$str.list, RemovedSourceData = RemovedSourceData, VeilYear = CurrentVeilYear)
  
  # Print current processing status:
  cat("Done")
  
  # Return output:
  return(Output)
  
}

#MRPDirectory <- "/Users/eargtl/Documents/Homepage/www.graemetlloyd.com/mrp" # MRP file directory
#XMLDirectory <- "/Users/eargtl/Documents/Homepage/www.graemetlloyd.com/xml" # XML file directory
#TargetClade <- "Cetacea"
#InclusiveDataList <- c("Aguirre-Fernandez_etal_2009a", "Aguirre-Fernandez_et_Fordyce_2014a", "Albright_etal_2018a", "Arnold_etal_2005a", "Benoit_etal_2011a", "Bianucci_2013a", "Bianucci_et_Gingerich_2011a", "Bianucci_et_Landini_2006a", "Bianucci_etal_2007a", "Bianucci_etal_2010a", "Bianucci_etal_2013a", "Bianucci_etal_2016a", "Bianucci_etal_2018aa", "Bianucci_etal_2018ab", "Bisconti_2008a", "Bisconti_2010a", "Bisconti_2015a", "Bisconti_et_Bosselaers_2016a", "Bisconti_etal_2013a", "Bisconti_etal_2019a", "Bisconti_2015a", "Boersma_et_Pyenson_2015a", "Boersma_et_Pyenson_2016a", "Boersma_etal_2017a", "Boessenecker_et_Fordyce_2015a", "Boessenecker_et_Fordyce_2015b", "Boessenecker_et_Fordyce_2017a", "Boessenecker_et_Fordyce_2015b", "Boessenecker_etal_2017a", "Bosselaers_et_Post_2010a", "Bouetel_et_de_Muizon_2006a", "Buono_et_Cozzuol_2013a", "Buono_etal_2017a", "Churchill_etal_2012a", "Churchill_etal_2016a", "Collareta_etal_2017a", "Colpaert_etal_2015a", "Demere_etal_2008a", "Dooley_etal_2004a", "Ekdale_etal_2011a", "El_Adli_etal_2014a", "Fajardo-Mellor_etal_2006a", "Fitzgerald_2010a", "Fordyce_1994a", "Fordyce_et_Marx_2013a", "Fordyce_et_Marx_2016a", "Fordyce_et_Marx_2018a", "Geisler_2001a", "Geisler_et_Luo_1996a", "Geisler_et_Sanders_2003a", "Geisler_et_Theodor_2009a", "Geisler_et_Uhen_2003a", "Geisler_et_Uhen_2005a", "Geisler_etal_2005a", "Geisler_etal_2011a", "Geisler_etal_2012a", "Geisler_etal_2014a", "Geisler_etal_2017a", "Gibson_etal_2018a", "Godfrey_etal_2016a", "Godfrey_etal_2017a", "Goldin_2018a", "Goldin_et_Startsev_2014a", "Goldin_et_Startsev_2017a", "Goldin_et_Steeman_2015a", "Goldin_et_Zvonok_2013a", "Goldin_etal_2014a", "Goldin_etal_2014b", "Hernandez_Cisneros_2018a", "Heyning_1997a", "Kimura_et_Hasegawa_2010a", "Lambert_2008a", "Lambert_2008b", "Lambert_et_Louwye_2016a", "Lambert_etal_2008a", "Lambert_etal_2009a", "Lambert_etal_2010a", "Lambert_etal_2013a", "Lambert_etal_2014a", "Lambert_etal_2015a", "Lambert_etal_2017a", "Lambert_etal_2017b", "Lambert_etal_2018a", "Lambert_etal_2018ba", "Lambert_etal_2018bb", "Lambert_etal_2019a", "Lambert_etal_inpressa", "Luo_et_Marsh_1996a", "Martinez-Caceres_etal_2017a", "Marx_2011a", "Marx_et_Fordyce_2015a", "Marx_et_Fordyce_2016a", "Marx_etal_2015a", "Marx_etal_2016a", "Marx_etal_2017a", "McGowen_etal_2009a", "Messenger_et_McGuire_1998a", "Mijan_etal_inpressa", "Muizon_etal_2019a", "Murakami_etal_2012a", "Murakami_etal_2012b", "Murakami_etal_2014a", "Murakami_etal_2014b", "Murakami_etal_inpressa", "Nelson_et_Uhen_inpressa", "OLeary_et_Gatesy_2008a", "Paolucci_etal_inpressa", "Peredo_et_Pyenson_2018a", "Peredo_et_Uhen_2016a", "Peredo_etal_inpressa", "Post_etal_2017a", "Pyenson_etal_2015a", "Racicot_etal_2019aa", "Racicot_etal_2019ab", "Ramassamy_2016a", "Sanders_et_Geisler_inpressa", "Solis-Anorve_etal_2019a", "Spaulding_etal_2009a", "Steeman_2007a", "Tanaka_et_Fordyce_2014a", "Tanaka_et_Fordyce_2015a", "Tanaka_et_Fordyce_2016a", "Tanaka_et_Fordyce_inpressa", "Tanaka_et_Watanabe_inpressa", "Tanaka_etal_2018a", "Theodor_et_Foss_2005a", "Thewissen_etal_2001a", "Thewissen_etal_2007a", "Tsai_et_Fordyce_2018a", "Tsai_et_Fordyce_inpressa", "Tsai_et_Fordyce_inpressb", "Uhen_1999a", "Uhen_2004a", "Uhen_2014a", "Uhen_et_Gingerich_2001a", "Velez-Juarbe_etal_2015a", "Velez-Juarbe_inpressa", "Viglino_etal_2019a", "Viglino_etal_inpressa", "Wichura_etal_2015a")
#ExclusiveDataList <- c("Averianov_inpressa", "Bravo_et_Gaete_2015a", "Brocklehurst_etal_2013a", "Brocklehurst_etal_2015aa", "Brocklehurst_etal_2015ab", "Brocklehurst_etal_2015ac", "Brocklehurst_etal_2015ad", "Brocklehurst_etal_2015ae", "Brocklehurst_etal_2015af", "Bronzati_etal_2012a", "Bronzati_etal_2015ab", "Brusatte_etal_2009ba", "Campbell_etal_2016ab", "Carr_et_Williamson_2004a", "Carr_etal_2017ab", "Frederickson_et_Tumarkin-Deratzian_2014aa", "Frederickson_et_Tumarkin-Deratzian_2014ab", "Frederickson_et_Tumarkin-Deratzian_2014ac", "Frederickson_et_Tumarkin-Deratzian_2014ad", "Garcia_etal_2006a", "Gatesy_etal_2004ab", "Grellet-Tinner_2006a", "Grellet-Tinner_et_Chiappe_2004a", "Grellet-Tinner_et_Makovicky_2006a", "Knoll_2008a", "Kurochkin_1996a", "Lopez-Martinez_et_Vicens_2012a", "Lu_etal_2014aa", "Norden_etal_inpressa", "Pisani_etal_2002a", "Ruiz-Omenaca_etal_1997a", "Ruta_etal_2003ba", "Ruta_etal_2003bb", "Ruta_etal_2007a", "Selles_et_Galobart_2016a", "Sereno_1993a", "Sidor_2001a", "Skutschas_etal_inpressa", "Tanaka_etal_2011a", "Toljagic_et_Butler_2013a", "Tsuihiji_etal_2011aa", "Varricchio_et_Jackson_2004a", "Vila_etal_2017a", "Wilson_2005aa", "Wilson_2005ab", "Zelenitsky_et_Therrien_2008a")
#HigherTaxaToCollapse = c()
#SpeciesToExclude = c()
#MissingSpecies = "exclude"
#Interval = NULL
#VeilLine = TRUE
#IncludeSpecimenLevelOTUs = TRUE
#BackboneConstraint = "McGowen_etal_2009a"
#MonophylyConstraint = NULL
#RelativeWeights = c(0, 100, 10, 1)
#WeightCombination = "sum"
#ReportContradictionsToScreen = FALSE
