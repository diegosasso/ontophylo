

################################
#### Path hamming functions ####
################################


#' @title Hamming distances for a list of trees
#'
#' @description Calculates hamming distances for all paths in each discretized tree of a list.
#'
#' @param tree.list multiSimmap object.
#'
#' @return A tibble with information on state changes, time spent on each state, edge IDs, absolute and normalized hamming distances for all edges and all trees in a list.
#'
#' @author Sergei Tarasov
#'
#' @examples
#' Tb.trees <- path_hamming_over_trees_KDE(tree.list)
#'
#' @export
path_hamming_over_trees_KDE <- function(tree.list) {

  Dist.trees <- tibble()
  i <- 1
  for (i in 1:length(tree.list)) {

    tr <- tree.list[[i]]
    tr.i <- path_hamming_over_all_edges(tr)
    tr.i <- mutate(tr.i, tree.id = i)
    tr.i <- mutate(tr.i, tree.tip.id = paste0(tree.id,'-', Focal.Edge.id))
    Dist.trees <- bind_rows(Dist.trees, tr.i)

  }

  return(Dist.trees)

}


#' @title Hamming distances for a tree
#'
#' @description Calculates hamming distances for all paths in a given discretized tree.
#'
#' @param tree.merge simmap object. A discretized simmap using the 'discr_Simmap' function where identical state bins were merged using the 'merge_tree_cat' function.
#'
#' @return A tibble with information on state changes, time spent on each state, edge IDs, absolute and normalized hamming distances for all edges in a tree.
#'
#' @author Sergei Tarasov
#'
#' @examples
#' Tb.trees <- path_hamming_over_trees_KDE(tree.list)
#'
#' @export
path_hamming_over_all_edges <- function(tree.merge) {

  #n.tips <- length(tree.merge$tip.label)
  #Nodes <- tree.merge$edge[,2] %>% unique()
  #n.tips <- nrow(tree.merge$edge)
  Dist.tips <- tibble()
  i <- 2
  #for (i in 1:n.tips){
  for (i in 1:nrow(tree.merge$edge)) {

    Node <- tree.merge$edge[i,2]
    P <- get_states_path(tree.merge, node = Node)
    Dist <- path_hamming(P)

    # add Poissson cumulative changes
    Dist <- mutate(Dist, Pois.count = c(1:nrow(Dist)))

    Dist.tips <- bind_rows(Dist.tips, mutate(Dist, Focal.Edge.id = i))

  }

  return(Dist.tips)

}


# OBS.: Name change: get_states_path_v2 => get_states_path #
# OBS.: get_states_path is based on get_states_path #
# OBS.: Maybe the function needs to be checked again #

#' Get state information about a given path.
#'
#' Internal function. Not exported.
#'
#' @author Sergei Tarasov
#'
get_states_path <- function(tree.merge, node) {

  Maps <- tree.merge$maps
  # get absolute ages
  #H <- phytools:::nodeHeights(tree.merge)
  #age.loc <- lapply(Maps, function(x) cumsum(x))
  #Maps.abs <-mapply(function(x,y) x+y, age.loc, H[,1] )

  # get state path
  tree.merge$edge
  #get_path_edges(tree.merge, node=95)
  E.id <- get_path_edges(tree.merge, node)
  edges.foc <- Maps[E.id]
  state.path <- lapply(edges.foc, function(x) rev(x)) %>% unlist %>% rev

  # make edge ids
  E.as.id <- mapply(function(MM, EE) rep(EE, length(MM)), MM = edges.foc, EE = E.id) %>% unlist %>% rev
  #

  Cum <- c(0, state.path)
  Cum <- cumsum(Cum)
  Int <- cbind(Cum[-length(Cum)], Cum[-1]) %>% as_tibble()
  colnames(Int) <- c('t.start', 't.end')
  Int <- mutate(Int,  States = names(state.path))

  Int <- mutate(Int, Edge.id = E.as.id)

  #### Remove duplicated successive states
  # Int <-tibble(t.start=c(1:11), t.end=c(2:12),
  #              States=c('0001', '0001', '0001', '0002', '0001', '0003', '0003', '0004', '0004', '0004', '0004'),
  #              Edge.id=c(1,1,1,1, 2,2, 3, 3,3,3,4))
  #
  SS <- 1
  while (SS < nrow(Int)) {

    if (Int$States[SS] == Int$States[SS+1]) {

      # merge two identical states
      Int$t.end[SS] <- Int$t.end[SS+1] # reassign time
      Int$Edge.id[SS] <- Int$Edge.id[SS+1] # reassign edge.id; edge.id shows the id of edge when change happens
      Int <- Int[-(SS+1),] # remove duplicate
      SS <- SS

    } else

      SS <- SS+1

  }
  ###

  Int <- mutate(Int,  delta.t = t.end-t.start)

  return(Int)

}


#' Get edges IDs from root to a given node.
#'
#' Internal function. Not exported.
#'
#' @author Sergei Tarasov
#'
get_path_edges <- function(tree.merge, node) {

  E <- tree.merge$edge
  E.des <- c()

  while (length(node) > 0) {

    e.id <- which(E[,2] == node)
    E.des <- c(E.des,e.id)
    node <- E[e.id, 1]

  }

  return(E.des)

}


# OBS.: Note that some states may be the same if they are separated by split #

#' @title Path hamming
#'
#' @description Calculates the hamming distance between states for a given path.
#'
#' @param Path data.frame. A tibble with state information about a given path (from root to a given node).
#' The tibble is the output obtained from the get_states_path function. The columns give information on state changes, time spent on each state, and edge IDs.
#'
#' @return The input tibble with two additional columns giving information on absolute and normalized hamming distances.
#'
#' @author Sergei Tarasov
#'
#' @examples
#' hm <- path_hamming(Path)
#'
#' @export
path_hamming <- function(Path) {

  st <- Path$States
  ham.d <- stringdist::stringdist(st[1], st, method = c("hamming"))
  # normalized hamming
  str.len <- nchar(st[1])
  ham.d.norm <- ham.d/str.len

  return(mutate(Path, Ham.dist = ham.d, Ham.dist.n = ham.d.norm))

}


########################
#### NHPP functions ####
########################


#' @title Get NHPP data for all edges (Markov KDE)
#'
#' @description Gets data on changing times between states for all edges of a given sample of trees for the Markov kernel density estimator (KDE).
#'
#' @param Tb.trees data.frame. A tibble obtained with the 'path_hamming_over_trees_KDE' function.
#'
#' @return A list with changing times between states for all edges of a given sample of trees.
#'
#' @author Sergei Tarasov
#'
#' @examples
#' path.data <- make_data_NHPP_KDE_Markov_kernel(Tb.trees)
#'
#' @export
make_data_NHPP_KDE_Markov_kernel <- function(Tb.trees) {

  Nedges <- Tb.trees$Focal.Edge.id %>% unique() %>% length()

  dt.out <- list()
  i <- 1
  for (i in 1:Nedges) {

    #dt.out[[i]] <- make_data_NHPP_over_tip(Tb.trees, tip_id = i)
    dt.out[[i]] <- make_data_NHPP_over_edge_MarkovKDE(Tb.trees, Focal.Edge = i)

  }

  return(dt.out)

}


#' @title Get NHPP data for a given edge (Markov KDE)
#'
#' @description Gets data on changing times between states for a given edge of a given sample of trees for the Markov kernel density estimator (KDE).
#'
#' @param Tb.trees data.frame. A tibble obtained with the 'path_hamming_over_trees_KDE' function.
#' @param Focal.Edge integer. A value indicating the edge ID.
#'
#' @return A numeric vector with changing times between states for a given edge.
#'
#' @author Sergei Tarasov
#'
#' @examples
#' make_data_NHPP_over_edge_MarkovKDE(Tb.trees, Focal.Edge = 5)
#'
#' @export
make_data_NHPP_over_edge_MarkovKDE <- function(Tb.trees, Focal.Edge) {

  #Focal.Edge <- 1
  f.tip <- filter(Tb.trees, Focal.Edge.id == Focal.Edge)
  change.time <- f.tip$t.start
  change.time <- change.time[-which(change.time == 0)]

  return(change.time)

}


#' @title Add pseudodata
#'
#' @description Adds a vector of pseudodata to the path data obtained from the 'make_data_NHPP_KDE_Markov_kernel' function.
#'
#' @param Edge.groups list. A list with groups of edge IDs.
#' @param Pseudo.data numeric. A vector with values of pdeusodata.
#' @param Path.data numeric. A list of path data obtained from the 'make_data_NHPP_KDE_Markov_kernel' function.
#'
#' @return A list of path data with the pseudodata added.
#'
#' @author Sergei Tarasov
#'
#' @examples
#' add_pseudodata(Edge.groups, Pseudo.data, Path.data)
#'
#' @export
add_pseudodata <- function(Edge.groups, Pseudo.data, Path.data) {

  Pseudo.path.data <- vector(length = length(Path.data), mode = 'list')

  for (i in 1:length(Edge.groups)) {

    j <- 1
    for (j in 1:length(Edge.groups[[i]])) {

      Pseudo.path.data[[Edge.groups[[i]][j]]] <- c(Pseudo.data[[i]], Path.data[[Edge.groups[[i]][j]]])

    }

  }

return(Pseudo.path.data)

}


# OBS.: Uses path root-tips for estimation #

#' @title Estimate bandwidth
#'
#' @description Estimate the bandwidth for the Markov KDE.
#'
#' @param tree.discr simmap or phylo object. A discretized tree using the 'discr_Simmap' function.
#' @param data.path list. A list of path data obtained from the 'make_data_NHPP_KDE_Markov_kernel' function.
#' @param band.width character. Bandwidth selectors for the KDE, as in density.
#'
#' @return A numeric vector.
#'
#' @author Sergei Tarasov
#'
#' @examples
#' estimate_band_W(tree_discr, data.path, band.width = "bw.nrd0")
#'
#' @export
estimate_band_W <- function(tree.discr, data.path, band.width = c('bw.nrd0', 'bw.nrd0', 'bw.ucv', 'bw.bcv', 'bw.SJ')) {

  otus <- c(1:length(tree.discr$tip.label))
  #tree.discr$edge[,2] %in% otus
  Edges <- match(otus, tree.discr$edge[,2])

  Band.width <- c()
  i <- 14
  for (i in Edges) {

    cat('Working on edge ', i, ' \n')
    dt <- data.path[[i]]

    if (band.width == 'bw.nrd')
      h <- bw.nrd(dt)
    if (band.width == 'bw.nrd0')
      h <- bw.nrd0(dt)
    if (band.width == 'bw.ucv')
      h <- bw.ucv(dt, lower = 0.01, upper = 20)
    if (band.width == 'bw.bcv')
      h <- bw.bcv(dt, lower = 0.01, upper = 20)
    if (band.width == 'bw.SJ')
      h <- bw.SJ(dt, lower = 0.01, upper = 20)

    Band.width <- c(Band.width, h)

  }

  return(Band.width)

}


#' @title Estimate the normalized Markov KDE
#'
#' @description Estimated the normalized Markov KDE for each edge averaged across all possible root-tip paths.
#'
#' @param tree.discr simmap or phylo object. A discretized tree using the 'discr_Simmap' function.
#' @param Path.data numeric. A list of path data obtained from the 'make_data_NHPP_KDE_Markov_kernel' function.
#' @param h numeric. A value for the bandwidth calculated using the 'estimate_band_W' function.
#'
#' @return A list with the estimated unnormalized ($Maps.mean) and normalized ($Maps.mean.norm) KDEs for each edge.
#'
#' @author Sergei Tarasov
#'
#' @examples
#' Edge.KDE <- estimate_edge_KDE(tree.discr, Pseudo.path.data, h = 15)
#'
#' @export
estimate_edge_KDE <- function(tree.discr, Path.data, h) {

  #Edge.KDE <- estimate_edge_KDE_unnorm(tree.discr, root.taxon=87, Path.data, band.width)
  Edge.KDE <- estimate_edge_KDE_Markov_kernel_unnorm(tree.discr, Path.data = Path.data , h = h)
  # normalize KDE
  #res <- 1000
  #Tmax <- max(phytools:::nodeHeights(tree.discr))
  #ss <- 0:res/res * max(phytools:::nodeHeights(tree.discr))
  #ss[1]-ss[2]
  #delta <- Tmax/1000

  Maps.dx <- tree.discr$maps
  Maps.mean <- Edge.KDE$Maps.mean # dy

  # estimations
  dy.dx <- mapply(function(x,y) y*x, y = Maps.mean, x = Maps.dx)
  sum.dy.dx <- lapply(dy.dx, sum)
  total.sum <- unlist(sum.dy.dx) %>% sum()

  Maps.mean.norm <- lapply(Maps.mean, function(x) x/total.sum)

  Edge.KDE$Maps.mean.norm <- Maps.mean.norm

  return(Edge.KDE)

}


#' @title Estimate the unnormalized Markov KDE
#'
#' @description Estimated the unnormalized Markov KDE for each edge averaged across all possible root-tip paths.
#'
#' @param tree.discr simmap or phylo object. A discretized tree using the 'discr_Simmap' function.
#' @param Path.data numeric. A list of path data obtained from the 'make_data_NHPP_KDE_Markov_kernel' function.
#' @param h numeric. A value for the bandwidth calculated using the 'estimate_band_W' function.
#'
#' @return A list with the estimated unnormalized KDEs ($Maps.mean) for each edge.
#'
#' @author Sergei Tarasov
#'
#' @examples
#' Edge.KDE <- estimate_edge_KDE_Markov_kernel_unnorm(tree.discr, Pseudo.path.data, h = 15)
#'
#' @export
estimate_edge_KDE_Markov_kernel_unnorm <- function(tree.discr, Path.data, h = 10) {

  #root.taxon=87
  Maps <- tree.discr$maps
  #Maps.mean <- Maps
  #Maps.error <- Maps
  # absolute ages
  #tree.discr$edge
  H <- phytools:::nodeHeights(tree.discr)
  age.loc <- lapply(Maps, function(x) cumsum(x))
  age.glob <- mapply(function(x,y) x+y, age.loc, H[,1])



  Density <- vector(length = length(Maps), mode = 'list')
  Edge <- 1
  for (Edge in 1:nrow(H)) {

    cat('Working on edge ', Edge, ' \n')
    X <- age.glob[[Edge]]
    Y <- Path.data[[Edge]]
    #
    #hist(Y, breaks = 50)

    #Y <- Y[Y>0] !!! UNMASK

    #
    #if (length(Y) == 0){Y <- 0}
    #Y <- 0
    #density(Y)
    #unique(Y)
    #
    if (length(Y) > 0) {

    Density[[Edge]] <- KDE_unnorm_trunc_Markov(x = X, h = h, dat = Y)

    }

    if (length(Y) == 0) {

      Density[[Edge]] <- rep(0, length(X))

      }

  }

  return(list(Maps.mean = Density))

}


#' KDE for unnormalized Markov kernel.
#'
#' Internal function. Not exported.
#'
#' @author Sergei Tarasov
#'
KDE_unnormalized_scalar_Markov_kernel <- function(x, h, dat) {

  y <- truncnorm::dtruncnorm(x, a = dat, b = Inf, mean = dat, sd = h) %>% sum()

  return(y)

}


#' KDE for unnormalized Markov kernel vectorized.
#'
#' Internal function. Not exported.
#'
#' @author Sergei Tarasov
#'
KDE_unnorm_trunc_Markov <- Vectorize(KDE_unnormalized_scalar_Markov_kernel, vectorize.args = 'x')


#' @title Get loess smoothing for the unnormalized Markov KDE
#'
#' @description Calculates loess smoothing for the unnormalized Markov KDE obtained from the 'estimate_edge_KDE_Markov_kernel_unnorm' function.
#'
#' @param tree.discr simmap or phylo object. A discretized tree using the 'discr_Simmap' function.
#' @param Edge.KDE list. A list with the estimated unnormalized KDEs ($Maps.mean) for each edge.
#'
#' @return A list with the loess smoothing calculated for each edge.
#'
#' @author Sergei Tarasov
#'
#' @examples
#' edge.KDE$Maps.mean.loess <- loess_smooting_KDE(tree.discr, edge.KDE)
#'
#' @export
loess_smoothing_KDE <- function(tree.discr, Edge.KDE) {

  # plot edge profiles
  Maps <- tree.discr$maps
  # absolute ages
  H <- phytools:::nodeHeights(tree.discr)
  age.loc <- lapply(Maps, function(x) cumsum(x))
  age.glob <- mapply(function(x,y) x+y, age.loc, H[,1] )

  Edges <- c(1:length(Edge.KDE$Maps.mean))


  Maps.mean.loess <- vector(length = length(Maps), mode = 'list')
  i <- 106
  for (i in Edges) {

    #for (i in 105:105){
    cat('Working on edge ', i, ' \n')
    dt <- tibble(X = age.glob[[i]], Y = Edge.KDE$Maps.mean[[i]])
    #dt <- tibble(X=age.glob[[i]], Y=Edge.KDE[[2]]$Maps.mean[[i]] )

    #loessMod <- loess(Y ~ X, data=dt, span=span)
    #loessMod <- loess.as(dt$X, dt$Y, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F)
    loessMod <- fANCOVA::loess.as(dt$X, dt$Y, degree = 1,
      criterion = c("aicc", "gcv")[2], user.span = NULL, plot = FALSE,
      control = loess.control(surface = "direct"))

    smoothed <- predict(loessMod)
    #str(smoothed)
    Maps.mean.loess[[i]] <- smoothed

  }

  return(Maps.mean.loess)

}


#' @title Normalize loess smoothing
#'
#' @description Normalizes the loess smoothing for the Markov KDE.
#'
#' @param tree.discr simmap or phylo object. A discretized tree using the 'discr_Simmap' function.
#' @param Maps.mean.loess list. A list with the loess smoothing calculated for each edge using the 'loess_smoothing_KDE' function.
#'
#' @return A list with the normalized loess smoothing calculated for each edge.
#'
#' @author Sergei Tarasov
#'
#' @examples
#' edge.KDE$Maps.mean.loess.norm <- normalize_KDE(tree.discr, edge.KDE$Maps.mean.loess)
#'
#' @export
normalize_KDE <- function(tree.discr, Maps.mean.loess) {

  Maps.dx <- tree.discr$maps
  Maps.mean <- Maps.mean.loess # dy

  # estimations
  dy.dx <- mapply(function(x,y) y*x, y = Maps.mean, x = Maps.dx)
  sum.dy.dx <- lapply(dy.dx, sum)
  total.sum <- unlist(sum.dy.dx) %>% sum()

  Maps.mean.norm <- lapply(Maps.mean, function(x) x/total.sum)

  #Edge.KDE$Maps.mean.norm <- Maps.mean.norm

  return(Maps.mean.norm)

}


#' @title Calculate KDE integral over edges
#'
#' @description Checks the integral of normalized Markov KDE or normalized loess smoothing over edges.
#'
#' @param tree.discr simmap or phylo object. A discretized tree using the 'discr_Simmap' function.
#' @param Edge.KDE.list list. A list with the normalized KDEs or loess smoothing for each edge.
#'
#' @return A numeric value for the integral over all edges.
#'
#' @author Sergei Tarasov
#'
#' @examples
#' integrate_edge_KDE(tree.discr, Edge.KDE$Maps.mean.norm)
#' integrate_edge_KDE(tree.discr, Edge.KDE$Maps.mean.loess.norm)
#'
#' @export
integrate_edge_KDE <- function(tree.discr, Edge.KDE.list) {

  #Edge.KDE <- estimate_edge_KDE_unnorm(tree.discr, root.taxon=87, Path.data)

  # normalize KDE
  #res <- 1000
  #Tmax <- max(phytools:::nodeHeights(tree.discr))
  #ss <- 0:res/res * max(phytools:::nodeHeights(tree.discr))
  #ss[1]-ss[2]
  #delta <- Tmax/1000

  Maps.dx <- tree.discr$maps
  #Maps.mean <- Edge.KDE$Maps.mean # dy
  Maps.mean <- Edge.KDE.list #dy

  # estimations
  dy.dx <- mapply(function(x,y) y*x, y = Maps.mean, x = Maps.dx)
  sum.dy.dx <- lapply(dy.dx, sum)
  total.sum <- unlist(sum.dy.dx) %>% sum()

  names(total.sum) <- 'Sum of all edges'

  return(total.sum)

}


#' @title Get analytical posterior
#'
#' @description Calculates the required statitics for the posterior distribution of number of state changes across all branches of all trees.
#'
#' @param tree.list multiSimmap object.
#'
#' @return A list with mean ($Mean), standard deviation ($SD), and 95HPD interval ($Q_2.5 and $Q_97.5) calculated for the posterior distribution.
#'
#' @author Sergei Tarasov
#'
#' @examples
#' tree.list <- Trees.focal.BRs[[1]]
#' posterior_lambda_KDE(tree.list)
#'
#' @export
posterior_lambda_KDE <- function(tree.list) {

  # uniq.trees <- Tb.trees$tree.id %>% unique() %>% length()
  # uniq.tips <- Tb.trees$tip.id %>% unique() %>% length()
  #
  # # n changes
  # n.changes <- c()
  # for(i in 1:uniq.trees){
  #   n.changes[i] <- filter(Tb.trees, tree.id == i) %>% nrow-1
  # }

  #class(tree.list) <- append(class(tree.list),"multiSimmap")
  n.changes <- lapply(tree.list, function(x) countSimmap(x)$N) %>% unlist


  # get required statistics
  Mean <- (0.01 + sum(n.changes))/(0.01 + length(n.changes))
  SD <- sqrt((0.01 + sum(n.changes))/(0.01 + length(n.changes))^2)
  Q_2.5 <- qgamma(p = 0.025, 0.01 + sum(n.changes), 0.01 + length(n.changes))
  Q_97.5 <- qgamma(p = 0.975, 0.01 + sum(n.changes), 0.01 + length(n.changes))

  XX <- list(Mean = Mean, SD = SD, Q_2.5 = Q_2.5, Q_97.5 = Q_97.5)

  return(XX)

}


#' @title Get distributions of analytical posterior
#'
#' @description Simulates a distribution of number of state changes across all branches of all trees
#'
#' @param tree.list multiSimmap object.
#' @param n.sim integer. Number of simulations.
#' @param BR.name character. A label name for the anatomical region.
#'
#' @return A tibble with the simulated distribution.
#'
#' @author Sergei Tarasov
#'
#' @examples
#' tree.list <- Trees.focal.BRs[[1]]
#' posterior_lambda_KDE_Distr(tree.list, n.sim = 10, "br")
#'
#' @export
posterior_lambda_KDE_Distr <- function(tree.list, n.sim = 10, BR.name) {

  # uniq.trees <- Tb.trees$tree.id %>% unique() %>% length()
  # uniq.tips <- Tb.trees$tip.id %>% unique() %>% length()
  #
  # # n changes
  # n.changes <- c()
  # for(i in 1:uniq.trees){
  #   n.changes[i] <- filter(Tb.trees, tree.id == i) %>% nrow-1
  # }

  #class(tree.list) <- append(class(tree.list),"multiSimmap")
  n.changes <- lapply(tree.list, function(x) countSimmap(x)$N) %>% unlist

  #
  sim <- rgamma(n.sim, 0.01 + sum(n.changes), 0.01 + length(n.changes))
  out <- tibble(BR = BR.name, sim)

  return(out)

}

#' @title Make posterior distribution of NHPP
#'
#' @description Produces a posterior distribution from a given list of statistics calculated with the 'posterior_lambda_KDE' function.
#'
#' @param Edge.KDE.stat list. A list with the estimated normalized or loess smoothing KDEs for each edge.
#' @param lambda.post list. A list with the distribution statistics calculated with the 'posterior_lambda_KDE' function.
#' @param lambda.post.stat character. A value with the statistic to be used.
#'
#' @return A list with the distribution of the selected statistic for each edge.
#'
#' @author Sergei Tarasov
#'
#' @examples
#' make_postPois_KDE(Edge.KDE.stat, lambda.post, lambda.post.stat = 'Mean')
#'
#' @export
make_postPois_KDE <- function(Edge.KDE.stat, lambda.post, lambda.post.stat = 'Mean') {

  #Maps <- tree.discr$maps
  # lambda
  lambda <- lambda.post[[lambda.post.stat]] # !!! Perhaps to include lambda right in Edge.KDE!!!
  #XX <- lapply(Edge.KDE$Maps.mean.norm, function(x) x*lambda)
  XX <- lapply(Edge.KDE.stat, function(x) x*lambda)

  return(XX)

}


#' @title Calculate KDE derivative over edges
#'
#' @description Calculates the derivative of the normalized Markov KDE or normalized loess smoothing over edges.
#'
#' @param tree.discr simmap or phylo object. A discretized tree using the 'discr_Simmap' function.
#' @param Edge.KDE.stat list. A list with the estimated normalized or loess smoothing KDEs for each edge.
#'
#' @return A list with the distribution of the derivatives calculated for each edge.
#'
#' @author Sergei Tarasov
#'
#' @examples
#' derivative_KDE(tree.discr, Edge.KDE.stat)
#'
#' @export
derivative_KDE <- function(tree.discr, Edge.KDE.stat) {

  # dx <- tree.discr$maps
  # y <-  Edge.KDE.stat # y
  #
  # # Calculate dy
  # #Maps.mean <- list(c(1,2,3,5), c(8,5,6))
  # dy <- lapply(Maps.mean, function(x) c(x[-1]-x[-length(x)], x[length(x)]-x[(length(x)-1)]) )
  #
  # # estimations
  # dy.dx <- mapply(function(x,y) y/x, y=dy, x=Maps.dx)

  #
  dx.list <- tree.discr$maps
  y.list <- Edge.KDE.stat # y
  E.map <- tree.discr$edge

  Map.deriv <- vector(length = nrow(E.map), mode = 'list')

  # get non-root edges
  roots <- c(1:nrow(E.map))[!(E.map[,1] %in% E.map[,2])]
  not.roots <- c(1:nrow(E.map))[(E.map[,1] %in% E.map[,2])]

  # Working on non-roots
  i <- 7
  for (i in not.roots) {

    y <- y.list[[i]]
    dx <- dx.list[[i]]

    # add first element from combination ancestor + current edge
    ances.edge <- which(E.map[,2] == E.map[i,1])
    #y.anc <- y.list[[ances.edge]][length(y.list[[ances.edge]])]
    #dx.anc <- dx.list[[ances.edge]][length(dx.list[[ances.edge]])]
    y.anc <- y.list[[ances.edge]]
    dy.last.anc <- y.anc[[length(y.anc)]]-y.anc[[length(y.anc)-1]]
    dx.last.anc <- dx.list[[ances.edge]][length(dx.list[[ances.edge]])]
    Der.last.anc <- dy.last.anc/dx.last.anc

    # average between ancestor and first derivatives
    dy <- y[-1]-y[-length(y)]
    dx <- dx[-1]
    Der <- dy/dx
    Der <- c((Der[1]+Der.last.anc)/2, Der)
    #dy.first <- y[1]-y.anc
    #dy.first <- (dy.last.anc+dy[[1]])/2
    #dy <- c(dy.first, dy)

    #dx.first <- dx[1]+dx.anc
    #dx[[1]] <- dx.first
    #

    Map.deriv[[i]] <- Der
    #Map.deriv[[i]] <- dy/dx

  }

  # Working on roots, just duplicating n+1 derivative
  i <- 1
  for (i in roots) {

    y <- y.list[[i]]
    dx <- dx.list[[i]]

    dy <- y[-1]-y[-length(y)]
    dy <- c(dy[[1]], dy)

    Map.deriv[[i]] <- dy/dx

  }

  return(Map.deriv)

}


#' @title Make contMap KDE object
#'
#' @description Produces a contMap object for plotting the NHPP.
#'
#' @param tree.discr phylo object. A discretized tree using the 'discr_Simmap' function.
#' @param Edge.KDE.stat list. A list with the distributions of the estimated parameter of KDEs for each edge.
#'
#' @return A 'contMap' object.
#'
#' @author Sergei Tarasov
#'
#' @examples
#' NHPP.Map <- make_contMap(tree.discr, Edge.KDE$Maps.mean.norm)
#'
#' @export
make_contMap_KDE <- function(tree.discr, Edge.KDE.stat) {

  Maps <- tree.discr$maps

  # lambda
  #lambda <- lambda.post[[lambda.post.stat]] # !!! Perhaps to include lambda right in Edge.KDE!!!
  #Maps.mean <- lapply(Edge.KDE$Maps.mean.norm, function(x) x*lambda)
  #Maps.mean <- Edge.KDE[[lambda.post.stat]]
  Maps.mean <- Edge.KDE.stat

  # absolute ages
  #H <- phytools:::nodeHeights(tree.discr)
  #age.loc <- lapply(Maps, function(x) cumsum(x))
  #age.glob <-mapply(function(x,y) x+y, age.loc, H[,1] )


  # assign colors to Maps.mean
  Maps.col <- Maps
  MM <- Maps.mean %>% unlist
  lims <- c(min(MM), max(MM))

  cols <- rainbow(1001, start = 0, end = 0.7) %>% rev
  names(cols) <- 0:1000
  trans <- 0:1000/1000 * (lims[2] - lims[1]) + lims[1]
  names(trans) <- 0:1000

  i <- 1
  for (i in 1:length(Maps.col)) {

    b <- Maps.mean[[i]]
    d <- sapply(b, phytools:::getState, trans = trans)
    names(Maps.col[[i]]) <- d

  }
  #
  tree.new <- tree.discr
  tree.new$maps <- Maps.col
  #tree.new$maps.sd <- Maps.sd

  xx <- list(tree = tree.new, cols = cols, lims = lims)
  class(xx) <- "contMap"

  return(xx)

}


####################
#### Edge plots ####
####################


#' @title Make edge profiles for plotting
#'
#' @description Gets the information necessary for making an edgeplot, where the tree is plotted in a space where the x axis is the time and y axis the scale of the desired statistics.
#'
#' @param tree.discr phylo object. A discretized tree using the 'discr_Simmap' function.
#' @param Edge.KDE.stat list. A list with the distributions of the estimated parameter of KDEs for each edge.
#'
#' @return A tibble with X and Y coordinates and other information necessary for making an edgeplot.
#'
#' @author Sergei Tarasov
#'
#' @examples
#'
#' @export
edge_profiles4plotting <- function(tree.discr, Edge.KDE.stat) {

  # plot edge profiles
  Maps <- tree.discr$maps
  # absolute ages
  H <- phytools:::nodeHeights(tree.discr)
  age.loc <- lapply(Maps, function(x) cumsum(x))
  age.glob <- mapply(function(x,y) x+y, age.loc, H[,1])

  v.x <- age.glob %>% enframe(name = "edge.id", value = "X") %>% unnest(cols = c(X))
  #v.x$edge.id %>% unique()
  #v.y <- Edge.KDE$lambda.mean  %>% enframe(name = "edge.id", value = "Y") %>% unnest()
  v.y <- Edge.KDE.stat  %>% enframe(name = "edge.id", value = "Y") %>% unnest(cols = c(Y))

  #edge.profs <- bind_cols(v.x,v.y)
  edge.profs <- mutate(v.x,v.y)
  #edge.profs <-select(edge.profs, -edge.id1) %>% mutate(edge.type = 'main')
  edge.profs <- mutate(edge.profs, edge.type = 'main')

  #all(edge.profs$edge.id == edge.profs$edge.id1)
  # add joins for edges
  edge.profs.joint <- join_edges(tree.discr, edge.profs)
  edge.profs.joint <- mutate(edge.profs.joint, edge.type = 'joint')

  edge.profs <- bind_rows(edge.profs, edge.profs.joint)

  #
  # v.x <- age.glob %>% enframe(name = "edge.id", value = "X") %>% unnest()
  # #v.x$edge.id %>% unique()
  # v.y <- Edge.KDE$lambda.mean %>% enframe(name = "edge.id", value = "Y") %>% unnest()
  #
  # edge.profs <- bind_cols(v.x,v.y)
  #
  # # add joins for edges
  # edge.profs.joint <- join_edges(tree.discr, edge.profs)
  #
  # edge.profs <- bind_rows(edge.profs, edge.profs.joint)

  return(edge.profs)

}


#' Join neighboring edges in edge profiles.
#'
#' Internal function. Not exported.
#'
#' @author Sergei Tarasov
#'
join_edges <- function(tree.discr, edge.profs) {

  E.map <- tree.discr$edge

  # getting non-tip edges
  tips <- c(1:Ntip(tree.discr))
  E.map.nontip <- c(1:nrow(E.map))[!(E.map[,2] %in% tips)]

  #Non.tip.edge <- edge.profs$edge.id %>% unique()

  Joints <- tibble()
  i <- 14
  for (i in E.map.nontip ) {

    # focal edge, last X,Y coords
    parent.edge <- filter(edge.profs, edge.id == i)
    Par <- parent.edge[nrow(parent.edge),]

    # descendant eges
    descen.edges <- which(E.map[,1] == E.map[i,2])
    #descen.edges <- E.map[(E.map[,1] %in% E.map[i,2]),]
    Des1 <- filter(edge.profs, edge.id == descen.edges[1])[1,]
    Des2 <- filter(edge.profs, edge.id == descen.edges[2])[1,]

    # joints <- bind_rows(Par, Des1, Par, Des2)
    # joints <- mutate(joints, groups = paste0(i, c('a', 'a', 'b', 'b')))
    joints <- bind_rows(Par, Des1, Par, Des2)
    #joints <- mutate(joints, edge.id = paste0(i, c('a', 'a', 'b', 'b')))
    joints <- mutate(joints, edge.id = (i + c(0.1,0.1, .2,.2)))

    Joints <- bind_rows(Joints, joints)

  }

  return(Joints)

}

