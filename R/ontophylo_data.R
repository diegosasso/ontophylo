
#' Hymenoptera dated tree
#'
#' A phylogenetic tree modified from Klopfstein et al. (2013).
#' The tree was dated using penalized likelihood as implemented in TreePL (Smith & O’Meara, 2012) and
#' then pruned to contain 20 representative taxa used for the package examples. 
#'
#' @docType data
#'
#' @format A phylogenetic tree of class "phylo".
#'
#' @references Klopfstein, S., Vilhelmsen, L., Heraty, J.M., Sharkey, M. & Ronquist, F. (2013) The hymenopteran tree of life: evidence from protein-coding genes and objectively aligned ribosomal data. PLoS One, 8, e69344.
#' (\doi{10.1371/journal.pone.0069344})
#'
#' @examples
#' data(hym_tree)
"hym_tree"



#' Hymenoptera character matrix
#'
#' The character matrix from the first example data set used in the showcase analyses in Porto et al. (2023).
#' The dataset comprises 30 simulated binary characters and the original matrix was reduced to contain only the 
#' 20 representative taxa used for the package examples.
#'
#' @docType data
#'
#' @format A data table with 20 rows and 30 columns; each row indicates a species and each column a character.
#'
#' @references Porto, D.S., Uyeda, J., Mikó, I. & Tarasov, S. (2023) Supporting Data: ontophylo: Reconstructing the evolutionary dynamics of phenomes using new ontology-informed phylogenetic methods.
#' (\doi{10.5281/zenodo.10285424})
#'
#' @examples
#' data(hym_matrix) 
"hym_matrix"



#' Hymenoptera character annotations
#' 
#' The character annotations from the first example data set used in the showcase analyses in Porto et al. (2023).
#' The annotations comprise the character ids, ontology term ids, and term labels used for each character.
#'
#' @docType data
#'
#' @format A data table with three columns and 30 rows; "char_id" contains character ids (e.g., "CH1", "CH2", "CH3");
#' "onto_id" contains ontology term ids from the HAO ontology (e.g., "HAO:0000234", "HAO:0000101", "HAO:0000639");
#' "label" contains the respective ontology term labels (e.g., "cranium", "antenna", "mouthparts").
#' Rows indicate characters. 
#'
#' @references Porto, D.S., Uyeda, J., Mikó, I. & Tarasov, S. (2023) Supporting Data: ontophylo: Reconstructing the evolutionary dynamics of phenomes using new ontology-informed phylogenetic methods.
#' (\doi{10.5281/zenodo.10285424})
#'
#' @examples
#' data(hym_annot) 
"hym_annot"



#' Hymenoptera Anatomy Ontology (HAO)
#'
#' The anatomy ontology of Hymenoptera (Yoder et al. 2010). This same ontology was also used in Tarasov et al. (2022).
#' 
#' @docType data
#'
#' @format List containing various ontological relationships and terms.
#'
#' @references Tarasov, S., Mikó, I. & Yoder, M.J. (2022) ontoFAST: an r package for interactive and semi-automatic annotation of characters with biological ontologies. Methods in Ecology and Evolution, 13, 324–329.
#' (\doi{10.1111/2041-210X.13753})
#'
#' Yoder MJ, Mikó I, Seltmann KC, Bertone MA, Deans AR. 2010. A Gross Anatomy Ontology for Hymenoptera. PLoS ONE 5 (12): e15991.
#' (\doi{10.1371/journal.pone.0015991})
#'
#' \href{http://portal.hymao.org/projects/32/public/ontology/}{Hymenoptera Anatomy Ontology Portal}
#'
#' @examples
#' data(HAO)
"HAO"



#' Hymenoptera stochastic character maps
#'
#' List of stochastic character maps obtained from all characters of 
#' the first example data set used in the showcase analyses in Porto et al. (2023). 
#' Only 50 samples per character were included to reduce file size (originally 100 samples).
#' 
#' @docType data
#'
#' @format List containing 30 objects of class "multiSimmap".
#'
#' @references Porto, D.S., Uyeda, J., Mikó, I. & Tarasov, S. (2023) Supporting Data: ontophylo: Reconstructing the evolutionary dynamics of phenomes using new ontology-informed phylogenetic methods.
#' (\doi{10.5281/zenodo.10285424})
#'
#' @examples
#' data(hym_stm)
"hym_stm"



#' Hymenoptera amalgamated stochastic character maps
#'
#' List of amalgamated stochastic character maps obtained from  
#' the first example data set used in the showcase analyses in Porto et al. (2023). 
#' Character were amalgamated into groups of 10 characters representing 
#' three anatomical regions of Hymenoptera: "head", "mesosoma", and "metasoma".
#' Only 50 samples per character were included to reduce file size (originally 100 samples).
#'
#' @docType data
#'
#' @format List containing three objects of class "multiPhylo".
#'
#' @references Porto, D.S., Uyeda, J., Mikó, I. & Tarasov, S. (2023) Supporting Data: ontophylo: Reconstructing the evolutionary dynamics of phenomes using new ontology-informed phylogenetic methods.
#' (\doi{10.5281/zenodo.10285424})
#'
#' @examples
#' data(hym_stm_amalg)
"hym_stm_amalg"



#' Hamming distances object (Hymenoptera example)
#'
#' Hamming distances data object from the amalgamated characters of "head"
#' obtained from the first example data set used in the showcase analyses in Porto et al. (2023).
#' The Hamming distances data object was obtained using the function "path_hamming_over_trees_KDE".
#' This data object is used to run the examples of the ontophyo package.
#'
#' @docType data
#'
#' @format A data table with 11 columns; "t.start" contains the starting times of mapped states; 
#' "t.end" contains the ending times of mapped states; "States" contains the amalgamated states; 
#' "Edge.id" contains the edge ids onto the tree; "delta.t" contains the duration of each mapped state;
#' "Ham.dist" contains the Hamming distances from the root state;
#' "Ham.dist.n" contains the normalized Hamming distances from the root state;
#' "Pois.count" contains the number of discrete changes;
#' "Focal.Edge.id", "tree.id", and "tree.tip.id" contain internal ids used by the package functions.
#'
#' @references Porto, D.S., Uyeda, J., Mikó, I. & Tarasov, S. (2023) Supporting Data: ontophylo: Reconstructing the evolutionary dynamics of phenomes using new ontology-informed phylogenetic methods.
#' (\doi{10.5281/zenodo.10285424})
#'
#' @examples
#' data(hym_hm)
"hym_hm"



#' Changing times object (Hymenoptera example)
#'
#' Changing times data object from the amalgamated characters of "head"
#' obtained from the first example data set used in the showcase analyses in Porto et al. (2023).
#' The changing times data object was obtained using the function "make_data_NHPP_KDE_Markov_kernel".
#' This data object is used to run the examples of the ontophyo package.
#'
#' @docType data
#'
#' @format List containing changing times between states for all edges of the tree sample.
#'
#' @references Porto, D.S., Uyeda, J., Mikó, I. & Tarasov, S. (2023) Supporting Data: ontophylo: Reconstructing the evolutionary dynamics of phenomes using new ontology-informed phylogenetic methods.
#' (\doi{10.5281/zenodo.10285424})
#'
#' @examples
#' data(hym_nhpp)
"hym_nhpp"



#' pNHPP rates object (Hymenoptera example)
#'
#' The data object contains pNHPP rates estimated for the amalgamated characters of "head"
#' obtained from the first example data set used in the showcase analyses in Porto et al. (2023).
#' The data object contains the estimated rates for all edges of the tree sample;
#' $Maps.mean and $Maps.mean.norm contain the raw and normalized rates estimated using the function "estimate_edge_KDE"; 
#' $Maps.mean.loess and $Maps.mean.loess.norm contain the smoothed raw and normalized rates estimated using the function "loess_smoothing_KDE";
#' $lambda.mean, $loess.lambda.mean, and $loess.lambda.mean.deriv contain the posteriors estimated for the raw and normalized rates, and its derivative.
#' This data object is used to run the examples of the ontophyo package.
#'
#' @docType data
#'
#' @format List containing pNHPP rates estimated for all edges of the tree sample.
#'
#' @references Porto, D.S., Uyeda, J., Mikó, I. & Tarasov, S. (2023) Supporting Data: ontophylo: Reconstructing the evolutionary dynamics of phenomes using new ontology-informed phylogenetic methods.
#' (\doi{10.5281/zenodo.10285424})
#'
#' @examples
#' data(hym_kde)
"hym_kde"



#' Hymenoptera amalgamated phenome 
#'
#' Example of a stochastic character map obtained from the amalgamation of 394 characters
#' modified from the matrix of Sharkey et al. (2011) and using the tree from Klopfstein et al. (2013).
#' The tree was dated using penalized likelihood as implemented in TreePL (Smith & O’Meara, 2012).
#' This data object is used to run the examples of the ontophyo morphospace application.
#'
#' @docType data
#'
#' @format An stochastic character map of class "simmap".
#'
#' @references Sharkey, M.J., et al. 2011. Phylogenetic relationships among superfamilies of Hymenoptera. Cladistics 28(1), 80-112.
#' (\doi{10.1111/j.1096-0031.2011.00366.x})
#'
#' Klopfstein, S., Vilhelmsen, L., Heraty, J.M., Sharkey, M. & Ronquist, F. (2013) The hymenopteran tree of life: evidence from protein-coding genes and objectively aligned ribosomal data. PLoS One, 8, e69344.
#' (\doi{10.1371/journal.pone.0069344})
#'
#' @examples
#' data(hym_stm_mds)
"hym_stm_mds"



#' Hymenoptera graphics information
#'
#' Example of annotations of anatomy terms from the Hymenoptera Anatomy Ontology (HAO)
#' to layers representing anatomical entities in a vector image of a hymenopteran wasp.
#' This data object is used to run the examples of the ontophyo package.
#' 
#' @docType data
#'
#' @format A data table with three columns; "Term" contains the ontology term labels (e.g., "cranium", "antenna");
#' "ID" contains the respective ontology term ids from the HAO ontology (e.g., "HAO:0000234", "HAO:0000101");
#' "pic_id" contains the layer ids of the corresponding anatomical entities in the vector image of a hymenopteran wasp.
#'
#' @examples
#' data(hym_graph)
"hym_graph"



#' Hymenoptera vector image
#'
#' Example of a vector image of a hymenopteran wasp.
#' The original vector image in PostScript format was converted to XML format using (\href{https://www.ghostscript.com/}{GhostScript})
#' and the function "PostScriptTrace" from `grImport` (Murrell, 2009).
#' Then the XML file was imported into R using the function "readPicture" also from `grImport`.
#' This data object is used to run the examples of the ontophyo package.
#'
#' @docType data
#'
#' @format A data object of class "grImport".
#'
#' @references Murrell, P. (2009). Importing vector graphics: The grimport package for r. Journal of Statistical Software, 30:1–37.
#' (\doi{10.18637/jss.v030.i04})
#'
#' @examples
#' data(hym_img)
"hym_img"
