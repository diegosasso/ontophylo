

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
#' @examples
#' GR <- data(HYM_IMG)$info
#' ontology <- HAO
#' get_vector_ids_per_term(term = 'HAO:0000349', ontology, GR)
#'
#' @export
get_vector_ids_per_term <- function(term = 'HAO:0000349', ONT, GR) {

  GR <- as_tibble(GR)

  GL <- vector("list", nrow(GR))
  names(GL) <- GR$ID

  for (i in 1:nrow(GR)){

    GL[[i]]$name <- GR[[1]][i]
    GL[[i]]$layer <- strsplit(GR[[3]][i], ", ")[[1]]

  }

  all.ids <- names(GL)

  des <- get_descendants(ONT, roots = term, exclude_roots = FALSE)
  pp <- all.ids[(all.ids %in% des)]
  selected.ids <- lapply(GL[pp], function(x) x$layer) %>% unlist %>% as.numeric()
  selected.ids <- selected.ids[!is.na(selected.ids)]

  return(selected.ids)

}


#' @title Wrapper for getting vector layer IDs for multiple terms
#'
#' @description Given an ontology_index object, data.frame with ontology term labels, and data.frame with picture information (see examples),
#' produces a named vector with layer IDs to be used in the 'make_pic' function.
#'
#' @param terms data.frame. A data.frame with ontology terms to get layer IDs for.
#' The first column corresponds to the ontology term labels, the second to the ontology IDs.
#' @param ONT ontology_index object.
#' @param GR data.frame. A data.frame with the picture information. It contains the matches between all ontology term labels and layer IDs in the Picture object.
#' The first column corresponds to the ontology term labels, the second to the ontology IDs, and the third to the layer IDs in the Picture object.
#'
#' @return A named vector with the layer IDs corresponding to or descending from the ontology term label queried.
#'
#' @author Diego S. Porto
#'
#' @examples
#' GR <- data(HYM_IMG)$info
#' ontology <- HAO
#' get_vector_ids_list(term = 'HAO:0000349', ontology, GR)
#'
#' @export
get_vector_ids_list <- function(terms, ONT, GR) {

  layers <- lapply(terms[[2]], function(x) get_vector_ids_per_term(term = x, ONT, GR))

  names(layers) <- terms[[1]]

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
#' @examples
#' #Stat <- setNames(runif(5, 0.1,10), c("cranium", "fore_wing", "hind_wing", "pronotum", "propectus") )
#' #hm.palette <- colorRampPalette(brewer.pal(9, 'Spectral') %>% rev, space = 'Lab')
#' #cols.maps <- make_colors(Stat, palette = hm.palette(100))
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
  odr <- (1+((ncols - 1)/(lims[2] - lims[1]))*(Stat-lims[1])) %>% round(.,0)
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
#' @examples
#' #Stat <- setNames(runif(5, 0.1,10), c("cranium", "fore_wing", "hind_wing", "pronotum", "propectus") )
#' #hm.palette <- colorRampPalette(brewer.pal(9, 'Spectral') %>% rev, space = 'Lab')
#' #cols.maps <- make_colors_relative_scale(Stat, palette = hm.palette(100), lims = c(min(Stat),max(Stat)))
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
  odr <- (1+((ncols - 1)/(lims[2] - lims[1]))*(Stat-lims[1])) %>% round(.,0)
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
#' @param nticks integer. Number of subdivisions of the scale.
#' @param title character. A legend for the scale bar.
#'
#' @return A plot of the color scale bar.
#'
#' @author Sergei Tarasov
#'
#' @examples
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
#' @param picture Picture object.
#' @param layers numeric. A named vector where values indicate the layer IDs in the Picture object and names indicate the anatomy ontology term labels.
#' @param cols.maps character. A named vector where elements correspond to color IDs and names indicate the anatomy ontology term labels.
#'
#' @return An object of class 'Picture' with the assigned colors to different anatomical regions.
#'
#' @author Sergei Tarasov
#'
#' @examples
#'
#' @export
make_pic <- function(picture, layers, cols.maps) {

  i <- 2
  for (i in 1:length(layers)) {

    picture@paths[layers[i]]$path@rgb <- cols.maps[names(cols.maps) == names(layers[i])]

  }

  return(picture)

}
