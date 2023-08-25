test_that(
  "Testing paramo()",
  {

    library(ontophylo)
    library(testthat)

    stm <- readRDS("../testdata/stm_test.RDS")

    tree_list <- list("CH1" = stm[1:2], "CH2" = stm[3:4], "CH3" = stm[5:6], "CH4" = stm[7:8])

    tree_list <- lapply(tree_list, function(x) discr_Simmap_all(x, res = 100))

    anat_query <- list("head" = c("CH1", "CH2"), "leg" = c("CH3", "CH4"))

    tree_list_amalg <- paramo(anat_query, tree_list, ntrees = 2)

    # Check number of trees.
    expect_true(all(sapply(tree_list, length) == sapply(tree_list_amalg, length)))

    # Check number of branches.
    n_br1 <- unlist(lapply(tree_list, function(x) lapply(x, function(y) length(y$maps) ) ))
    n_br2 <- unlist(lapply(tree_list_amalg, function(x) lapply(x, function(y) length(y$maps) ) ))
    expect_true(all(n_br1 == n_br1))

    # Check number of amalgamated characters.
    n_chars <- lapply(tree_list_amalg, function(x) lapply(x, function(y) names(unlist(y$maps)) ) )
    n_chars <- sapply(n_chars, function(x) unique(unlist(lapply(x, function(y) nchar(unique(y)) ))) )
    expect_identical(length(tree_list), sum(n_chars))

  }
)
