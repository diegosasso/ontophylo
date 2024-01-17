

###########################################################
#### Organize vector ID information and ontology terms ####
###########################################################


#' @title Get vector layer IDs for single term
#'
#' @description Given an ontology_index object, ontology term label, and data.frame with picture information (see examples),
#' produces a named vector with layer IDs to be used in the 'make_pic' function.
#'
#' @param term character. An ontology term label to get the corresponding layer IDs in the Picture object.
#' @param ONT ontology_index object.
#' @param GR data.frame. A data.frame with the picture information. It contains the matches between all ontology term labels and layer IDs in the Picture object.
#' The first column corresponds to the ontology term labels, the second to the ontology IDs, and the third to the layer IDs in the Picture object.
#'
#' @return A named vector with the layer IDs corresponding to or descending from the ontology term label queried.
#'
#' @author Sergei Tarasov
#'
#' @importFrom tibble as_tibble
#' @importFrom ontologyIndex get_descendants
#'
#' @examples
#' data("HAO", "hym_graph")
#' # Get picture layers from head.
#' get_vector_ids_per_term(term = "HAO:0000397", ONT = HAO, GR = hym_graph)
#'
#' @export
get_vector_ids_per_term <- function(term = 'HAO:0000349', ONT, GR) {

  GR <- as_tibble(GR)

  GL <- vector("list", nrow(GR))
  names(GL) <- GR[[2]]

  for (i in 1:nrow(GR)){

    GL[[i]]$name <- GR[[1]][i]
    GL[[i]]$layer <- strsplit(GR[[3]][i], ", ")[[1]]

  }

  all.ids <- names(GL)

  des <- get_descendants(ONT, roots = term, exclude_roots = FALSE)
  pp <- all.ids[(all.ids %in% des)]
  selected.ids <- lapply(GL[pp], function(x) x$layer) %>% unlist %>% as.numeric()
  selected.ids <- unique(selected.ids[!is.na(selected.ids)])

  return(selected.ids)

}


#' @title Wrapper for getting vector layer IDs for multiple terms
#'
#' @description Given an ontology_index object, data.frame with ontology term labels, and data.frame with picture information (see examples),
#' produces a named vector with layer IDs to be used in the 'make_pic' function.
#'
#' @param terms_list list. A named list with ontology terms to get layer IDs for.
#' The first column corresponds to the ontology term labels, the second to the ontology IDs.
#' @param ONT ontology_index object.
#' @param GR data.frame. A data.frame with the picture information. It contains the matches between all ontology term labels and layer IDs in the Picture object.
#' The first column corresponds to the ontology term labels, the second to the ontology IDs, and the third to the layer IDs in the Picture object.
#'
#' @return A named vector with the layer IDs corresponding to or descending from the ontology term label queried.
#'
#' @author Diego S. Porto
#'
#' @importFrom stats setNames
#'
#' @examples
#' data("HAO", "hym_graph")
#' # Get picture layers from three anatomical regions.
#' terms_list <- as.list(c("HAO:0000397", "HAO:0000576", "HAO:0000626"))
#' terms_list <- setNames(terms_list, c("head", "mesosoma", "metasoma"))
#' get_vector_ids_list(terms = terms_list , ONT = HAO, GR = hym_graph)
#'
#' @export
get_vector_ids_list <- function(terms_list, ONT, GR) {
  
  layers <- lapply(terms_list, function(x) get_vector_ids_per_term(term = x, ONT, GR))
  
  names(layers) <- names(terms_list)
  
  layers <- setNames(unlist(layers, use.names = F), rep(names(layers), lengths(layers)))
  
  return(layers)
  
}


####################
#### Get colors ####
####################


#' @title Make color palette for image plotting
#'
#' @description Produces a color scale for a given statistic of evolutionary rate.
#'
#' @param Stat numeric. A named vector where values are the statistics, and names are ontology term labels.
#' @param palette A character vector or function defining a color palette.
#'
#' @return A character vector where elements are color IDs and names are the input ontology term labels.
#'
#' @author Sergei Tarasov
#'
#' @importFrom stats setNames
#'
#' @examples
#' stat <- setNames(runif(5, 0.1, 10), 
#' c("cranium", "fore_wing", "hind_wing", "pronotum", "propectus") )
#' hm.palette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"), space = "Lab")
#' cols.maps <- make_colors(stat, palette = hm.palette(100))
#' cols.maps
#'
#' @export
make_colors <- function(Stat, palette) {

  # Colors
  ncols <- length(palette)
  lims <- c(min(Stat), max(Stat))

  #cols <- rainbow(ncols, start = 0, end = 0.7) %>% rev
  cols <- palette
  names(cols) <- 1:ncols

  # remap, scale colors
  odr <- (1+((ncols - 1)/(lims[2] - lims[1]))*(Stat-lims[1]))
  odr <- round(odr,0)
  cols <- cols[odr]

  cols.maps <- setNames(cols, names(Stat))

  return(cols.maps)

}


#' @title Make color palette for image plotting with relative scale
#'
#' @description Produces a relative color scale for a given statistic of evolutionary rate.
#'
#' @param Stat numeric. A named vector where values are the statistics, and names are ontology term labels.
#' @param palette A character vector or function defining a color palette.
#' @param lims numeric. A pair of values defining the lower and upper limits of the scale.
#'
#' @return A character vector where elements are color IDs and names are the input ontology term labels.
#'
#' @author Sergei Tarasov
#'
#' @importFrom stats setNames
#'
#' @examples
#' stat <- setNames(runif(5, 0.1, 10), 
#' c("cranium", "fore_wing", "hind_wing", "pronotum", "propectus") )
#' hm.palette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"), space = "Lab")
#' cols.maps <- make_colors_relative_scale(stat, palette = hm.palette(100), 
#' lims = c(min(stat), max(stat)))
#' cols.maps
#'
#' @export
make_colors_relative_scale <- function(Stat, palette, lims) {

  # Colors
  ncols <- length(palette)
  #lims <- c(min(Stat), max(Stat))

  #cols <- rainbow(ncols, start = 0, end = 0.7) %>% rev
  cols <- palette
  names(cols) <- 1:ncols

  # remap, scale colors
  odr <- (1+((ncols - 1)/(lims[2] - lims[1]))*(Stat-lims[1]))
  odr <- round(odr,0)  
  cols <- cols[odr]

  cols.maps <- setNames(cols, names(Stat))

  return(cols.maps)

}


#' @title Color bar
#'
#' @description Function to plot the color scale bar.
#'
#' @param pal character. A vector with color IDs.
#' @param min numeric. Value for lower limit of the scale.
#' @param max numeric. Value for upper limit of the scale.
#' @param ticks integer. A vector of values for the scale.
#' @param nticks integer. Number of subdivisions of the scale.
#' @param title character. A legend for the scale bar.
#'
#' @return A plot of the color scale bar.
#'
#' @author Sergei Tarasov
#'
#' @importFrom graphics rect axis
#'
#' @examples
#' stat <- runif(10, 0.25, 1)
#' hm.palette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"), space = "Lab")
#' color.bar(hm.palette(100), min = min(stat), max = max(stat),
#'           ticks = round(c(min(stat), max(stat)/2, max(stat)), 2), title = "")
#'
#' @export
color.bar <- function(pal, min, max = -min, nticks = 11, ticks = seq(min, max, len = nticks), title = '') {

  scale = (length(pal)-1)/(max-min)
  #dev.new(width=1.75, height=5)
  #par(mar=c(bottom=3, left=3, top=19, right=35))
  #plot(0:4,  type="n", axes=FALSE, xlab = "", ylab = "", bty="n")
  plot(c(0,1), c(min,max),  type = "n",  bty = 'n', xlab = '',  ylab = '', main = title, xaxt = 'n',  yaxt = 'n')
  #plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las = 1)

  for (i in 1:(length(pal)-1)) {

    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col = pal[i], border = NA)

  }

}


#######################
#### Build picture ####
#######################


#' @title Assign colors to picture ID layers
#'
#' @description Assigns colors to picture ID layers (@paths) of an object of class 'Picture'.
#' The object should be a PS or ESP vector illustration imported using the grImport package.
#' Colors are taken from cols.maps argument were the palette indicates the scale of the desired statistics for the evolutionary rates.
#'
#' @param picture grImport object. A vector image imported in R using the 'readPicture' function from grImport.
#' @param layers numeric. A named vector where values indicate the layer IDs in the Picture object and names indicate the anatomy ontology term labels.
#' @param cols.maps character. A named vector where elements correspond to color IDs and names indicate the anatomy ontology term labels.
#'
#' @return An object of class 'Picture' with the assigned colors to different anatomical regions.
#'
#' @author Sergei Tarasov
#'
#' @examples
#' data("HAO", "hym_graph", "hym_img", "hym_kde")
#' # Get picture.
#' picture <- hym_img
#' # Get picture layers from three anatomical regions.
#' terms_list <- as.list(c("HAO:0000397", "HAO:0000576", "HAO:0000626"))
#' terms_list <- setNames(terms_list, c("head", "mesosoma", "metasoma"))
#' anat_layers <- get_vector_ids_list(terms = terms_list , ONT = HAO, GR = hym_graph)
#' # Get mean rates all branches for the three anatomical regions.
#' plot_stat <- lapply(hym_kde, function(x) unlist(lapply(x$loess.lambda.mean, function(x) mean(x) )) )
#' plot_stat <- do.call(cbind, plot_stat)
#' # Add two columns for the other anatomical regions (just for this example).
#' plot_stat <- cbind(plot_stat, plot_stat*0.75, plot_stat*0.5)
#' colnames(plot_stat) <- c("head", "mesosoma", "metasoma")
#' # Select an arbitrary branch.
#' plot_stat <- plot_stat[5,]
#' # Set scale.
#' scale_lim <- range(plot_stat)
#' # Get color palette.
#' hm.palette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"), space = "Lab")
#' cols_maps <- make_colors_relative_scale(plot_stat, palette = hm.palette(100),
#'                                         lims = scale_lim)
#' # Plot picture.
#' new_pic <- make_pic(picture, anat_layers, cols_maps)
#' grImport::grid.picture(new_pic)
#'
#' @export
make_pic <- function(picture, layers, cols.maps) {

  i <- 2
  for (i in 1:length(layers)) {

    picture@paths[layers[i]]$path@rgb <- cols.maps[names(cols.maps) == names(layers[i])]

  }

  return(picture)

}


######################
#### Plot picture ####
######################


#' @title Plot Picture
#'
#' @description Wrapper function for making a plot of an object of class 'Picture' using the 'make_pic' function.
#'
#' @param picture grImport object. A vector image imported in R using the 'readPicture' function from grImport.
#' @param anat_layers numeric. A named vector with the layer IDs obtained using the 'get_vector_ids_list' function.
#' @param plot_stat numeric. A named vector with values corresponding to the branch statistics to be plotted and names matching the names in the layer IDs.
#' @param color_palette A character vector or function defining a color palette.
#' @param scale_lim numeric. A pair of values defining the lower and upper limits of the scale.
#'
#' @return A plot of the object of class 'Picture' with the assigned colors to different anatomical regions.
#'
#' @author Diego Porto
#'
#' @import dplyr
#' @importFrom grImport grid.picture
#'
#' @examples
#' data("HAO", "hym_graph", "hym_img", "hym_kde")
#' # Get picture.
#' picture <- hym_img
#' # Get picture layers from three anatomical regions.
#' terms_list <- as.list(c("HAO:0000397", "HAO:0000576", "HAO:0000626"))
#' terms_list <- setNames(terms_list, c("head", "mesosoma", "metasoma"))
#' anat_layers <- get_vector_ids_list(terms = terms_list , ONT = HAO, GR = hym_graph)
#' # Get mean rates all branches for the three anatomical regions.
#' plot_stat <- lapply(hym_kde, function(x) unlist(lapply(x$loess.lambda.mean, function(x) mean(x) )) )
#' plot_stat <- do.call(cbind, plot_stat)
#' # Add two columns for the other anatomical regions (just for this example).
#' plot_stat <- cbind(plot_stat, plot_stat*0.75, plot_stat*0.5)
#' colnames(plot_stat) <- c("head", "mesosoma", "metasoma")
#' # Select an arbitrary branch.
#' plot_stat <- plot_stat[5,]
#' # Set scale.
#' scale_lim <- range(plot_stat)
#' # Get color palette.
#' hm.palette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"), space = "Lab")
#'
#' # Plot picture.
#' anat_plot(picture, anat_layers, plot_stat, hm.palette(100), scale_lim)
#' 
#' @export
anat_plot <- function(picture, anat_layers, plot_stat, color_palette, scale_lim) {

  # Set layout to print.
  layout(matrix(c(1:4), ncol = 2, nrow = 2),
         heights = c(0.7, 1), widths = c(0.2, 1))

  # Set color bar.
  color.bar(color_palette, scale_lim[[1]], scale_lim[[2]],
            ticks = round(c(scale_lim[[1]], scale_lim[[2]]/2, scale_lim[[2]]), 2), title = "")

  # Set color scale for pictures.
  cols_maps <- make_colors_relative_scale(plot_stat, palette = color_palette,
                                          lims = scale_lim)
  # Make picture.
  new_pic <- make_pic(picture, anat_layers, cols_maps)

  grid.picture(new_pic)

}
