

##################################################
#### Multidimensional scaling functions (MDS) ####
##################################################


#' @title Multidimensional scaling of character states from one stochastic character map
#'
#' @description Performs multidimensional scaling (MDS) based on hamming distances among character state vectors from one stochastic character map.
#'
#' @param tree.merge simmap object. A stochastic character map after being processed by the discr_Simmap and merge_tree_cat functions.
#' @param add.noise numeric. A vector of length 2 or NULL. Indicates if noise should be added or not. Useful if there are many identical states occupying the same points in the 2D coordinates of the morphospace plot.
#' The noise is calculated as var(V)*add.noise.
#'
#' @return A list of tibbles -- Points, Lines, and Edge.map -- correponding to tree branch information to plot.
#'
#' @author Sergei Tarasov
#'
#' @import dplyr
#' @importFrom tibble as_tibble
#' @importFrom purrr set_names
#' @importFrom stringdist stringdistmatrix
#' @importFrom stats cmdscale runif var
#'
#' @examples
#' data("hym_stm_mds")
#' # Select a few taxa from main lineages of Hymenoptera.
#' tax <- c("Xyela", "Tenthredo", "Orussus", "Pimpla",
#'          "Ceraphron", "Evania", "Pison",
#'          "Ibalia", "Proctotrupes", "Chiloe") 
# Prune original tree.
#' drop_tax <- hym_stm_mds$tip.label[!hym_stm_mds$tip.label %in% tax]
#' hym_stm_mds <- phytools::drop.tip.simmap(hym_stm_mds, drop_tax)
#' # Get a sample of amalgamated stochastic map (phenome).
#' tree <- merge_tree_cat(hym_stm_mds)
#' \donttest{
#'
#'   # Multidimensional scaling for an arbitrary tree.
#'   MD <- suppressWarnings(MultiScale.simmap(tree))
#'
#' }
#'
#' @export
MultiScale.simmap <- function(tree.merge, add.noise = NULL) {

  #cat('MultiScale.simmap')
  # recode states
  Maps <- tree.merge$maps
  lengths <- lapply(Maps, function(x) length(x)) %>% unlist()

  start <- lengths
  start <- cumsum(start)+1
  start <- c(1, start[-length(start)])
  end <- cumsum(lengths)

  states_id <- mapply(function(x,y) c(x:y), start, end)
  # get Map
  Map.new <- mapply(function(x,y) set_names(x,y), Maps, states_id)

  # Map 2 edges
  Map.new.names <- lapply(Map.new, function(x) names(x) %>% as.numeric)
  X <- lapply(Map.new.names, function(y) cbind(y[-length(y)], y[-1]))
  Ed.int <- c()
  for (i in 1:length(X)) {Ed.int <- rbind(Ed.int, X[[i]])}

  Ed.ext.end <- start #lapply(Map.new.names, function(x) x[1]) %>% unlist()
  Ed.ext.start <- end[match(tree.merge$edge[,1], tree.merge$edge[,2])]

  # make root
  root <- Ed.ext.end[is.na(Ed.ext.start)]
  # bind colls without NAs
  Ed.ext <- cbind(Ed.ext.start[!is.na(Ed.ext.start)], Ed.ext.end[!is.na(Ed.ext.start)])
  # add root
  Ed.ext <- rbind(root, Ed.ext)


  Ed.all <- rbind(Ed.int, Ed.ext) %>% as_tibble()

  # get those states which are t-1 for extant taxa which are t
  tax.id.org <- c(1:length(tree.merge$tip.label))
  tax.id.new <- end[match(tax.id.org, tree.merge$edge[,2] )]
  #


  # org states
  #Maps <- tree.merge$maps
  states <- lapply(Maps, function(x) names(x)) %>% unlist
  ham.d <- stringdistmatrix(states, states, method = c("hamming"))
  #colnames(ham.d) <- rownames(ham.d) <- states
  mds1 <- cmdscale(ham.d, k = 2)
  M <- as_tibble(mds1)

  # activated #
  # add noise to separate points with the same coords
  if (!is.null(add.noise)) {
     #add.noise=c(.1,.1)
     noise.x <- var(M$V1)*add.noise[1]
     noise.y <- var(M$V2)*add.noise[2]
     noise <- cbind(runif(nrow(M), -noise.x, noise.x), runif(nrow(M), -noise.y, noise.y))
     M <- (M+noise) %>% as_tibble()
  }


  # get absolute ages of change and join them with M
  H <- phytools::nodeHeights(tree.merge)
  age.loc <- lapply(Maps, function(x) cumsum(x))
  age.glob <- mapply(function(x,y) x+y, age.loc, H[,1])
  Age <- lapply(age.glob, function(x) x) %>% unlist
  M <- bind_cols(M, time = Age)

  # combine extant/extinct species with M
  spp <- rep('no', nrow(M))
  spp[tax.id.new] <- 'yes'
  M <- bind_cols(M, sp_extant = spp)
  # mark root (2 taxa since tree is not rooted)
  isroot <- rep('no', nrow(M))
  isroot[root] <- 'yes'
  M <- bind_cols(M, is_root = isroot)
  #Extant.taxa <- bind_cols(M[tax.id.new,], id=tax.id.new)

  # get start and end coordinates for lines
  L.start <- M[Ed.all$V1,]
  L.start <- rename(L.start, all_of(c(start.V1 = "V1", start.V2 = "V2", start.time = "time", start.sp_extant = "sp_extant")))
  L.end <- M[Ed.all$V2,]
  L.end <- rename(L.end, all_of(c(end.V1 = "V1", end.V2 = "V2", end.time = "time", end.sp_extant = "sp_extant")))
  Lines.coor <- bind_cols(L.start, L.end)

  MD <- list(Points = M, Lines = Lines.coor, Edge.map = Ed.all)

  return(MD)

}


#' @title Adding noise to MDS from one stochastic character map
#'
#' @description Adds noise to the points in the 2D coordinates in the MDS plot.
#' # The noise is calculated as var(V)*add.noise.
#'
#' @param MD tibble. The output of a MD_morpho function.
#' @param add.noise numeric. A vector of length 2 indicating the amount of noise to be added to the x and y coordinates.
#'
#' @return A list of tibbles as in the output of MD_morpho functions, with noise added.
#'
#' @author Sergei Tarasov
#'
#' @import dplyr
#' @importFrom stats runif
#'
#' @examples
#' data("hym_stm_mds")
#' # Select a few taxa from main lineages of Hymenoptera.
#' tax <- c("Xyela", "Tenthredo", "Orussus", "Pimpla",
#'          "Ceraphron", "Evania", "Pison",
#'          "Ibalia", "Proctotrupes", "Chiloe") 
# Prune original tree.
#' drop_tax <- hym_stm_mds$tip.label[!hym_stm_mds$tip.label %in% tax]
#' hym_stm_mds <- phytools::drop.tip.simmap(hym_stm_mds, drop_tax)
#' # Get a sample of amalgamated stochastic map (phenome).
#' tree <- merge_tree_cat(hym_stm_mds)
#'
#' \donttest{
#'
#'   # Multidimensional scaling for an arbitrary tree.
#'   MD <- suppressWarnings(MultiScale.simmap(tree))
#'
#'   # Add noise.
#'   add_noise_MD(MD, c(0.3, 0.3))
#'
#' }
#'
#' @export
add_noise_MD <- function(MD, add.noise) {

  #cat('Noise Simmap')

  # Update Points
  M <- MD$Points
  #noise.x <- var(M$V1)*add.noise[1]
  #noise.y <- var(M$V2)*add.noise[2]
  #noise <- cbind(runif(nrow(M), -add.noise[1], add.noise[1]), runif(nrow(M), -add.noise[2], add.noise[2]))
  M$V1 <- M$V1+runif(nrow(M), -add.noise[1], add.noise[1])
  M$V2 <- M$V2+runif(nrow(M), -add.noise[2], add.noise[2])
  MD$Points <- M

  # Update lines
  L.start <- M[MD$Edge.map$V1,]
  L.start <- rename(L.start, all_of(c(start.V1 = "V1", start.V2 = "V2", start.time = "time", start.sp_extant = "sp_extant")))
  L.end <- M[MD$Edge.map$V2,]
  L.end <- rename(L.end, all_of(c(end.V1 = "V1", end.V2 = "V2", end.time = "time", end.sp_extant = "sp_extant")))
  Lines.coor <- bind_cols(L.start, L.end)
  MD$Lines <- Lines.coor

  return(MD)

}


#' @title Plot morphospace from MDS
#'
#' @description Wrapper function for plotting morphospaces obtained using the MultiScale.simmap' function.
#'
#' @param MD list. A list with the morphospace information obtained using the 'MultiScale.simmap' function.
#' @param Tslice numeric. The value for the temporal slices to be plotted, from root to tip. For example, if Tslice = 25, then all points in the morphospaces from time 0 (root) to 25 will be plotted.
#'
#' @return An object of class ggplot with the morphospace to be plotted.
#'
#' @author Diego Porto
#'
#' @import dplyr
#' @import ggplot2
#'
#' @examples
#' # Select a few taxa from main lineages of Hymenoptera.
#' tax <- c("Xyela", "Tenthredo", "Orussus", "Pimpla",
#'          "Ceraphron", "Evania", "Pison",
#'          "Ibalia", "Proctotrupes", "Chiloe") 
# Prune original tree.
#' drop_tax <- hym_stm_mds$tip.label[!hym_stm_mds$tip.label %in% tax]
#' hym_stm_mds <- phytools::drop.tip.simmap(hym_stm_mds, drop_tax)
#' # Get a sample of amalgamated stochastic map (phenome).
#' tree <- merge_tree_cat(hym_stm_mds)
#' \donttest{
#'  
#'   # Multidimensional scaling for an arbitrary tree.
#'   MD <- suppressWarnings(MultiScale.simmap(tree))
#'
#'   MD_plot <- mds_plot(MD, Tslice = 10)
#'   MD_plot
#'   MD_plot <- mds_plot(MD, Tslice = 50)
#'   MD_plot
#'   MD_plot <- mds_plot(MD, Tslice = 200)
#'   MD_plot
#'   MD_plot <- mds_plot(MD, Tslice = 280)
#'   MD_plot
#'
#' }
#'
#' @export
mds_plot <- function(MD, Tslice = max(MD$Points$time)) {

  # Add tip ids.
  MD$Points <- mutate(MD$Points, tip.id = c(1:nrow(MD$Points)))

  # Get tree height.
  Tmax <- max(MD$Points$time)

  # Filter data according to time slice.
  MD$Points <- filter(MD$Points, MD$Points$time <= Tslice)

  # Generate MDS plot.
  GG <-

    ggplot(MD$Points) +

    geom_segment(data = MD$Lines, aes(x = MD$Lines$start.V1, y = MD$Lines$start.V2, xend = MD$Lines$end.V1, yend = MD$Lines$end.V2),
                 colour = "red", linewidth = 0.3, linetype = 1, alpha = 0.3) +

    geom_point(aes(x = MD$Points$V1, y = MD$Points$V2, color = MD$Points$time, group = MD$Points$tip.id), alpha = 0.8, size = 1.6) +

    scale_colour_gradient(low = "purple",  high = "green", limits = c(0, Tmax),
                          breaks = seq(50, Tmax, 50), labels = seq(50, Tmax, 50) %>% rev) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 16),
          axis.title.y = element_text(size = 18),
          axis.text.y = element_text(size = 16),) +

    xlab("Coor1") + ylab("Coor2") +
    guides(color = guide_colourbar(title = "Time Myr", reverse = F, draw.ulim = T, draw.llim = T))

}
