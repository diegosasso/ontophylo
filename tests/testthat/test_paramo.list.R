test_that(
  "Testing paramo.list()",
  {

    stm <- readRDS("../testdata/stm_test.RDS")

    tree_list <- list("CH1" = stm[1:2], "CH2" = stm[3:4], "CH3" = stm[5:6])

    tree_list <- lapply(tree_list, function(x) discr_Simmap_all(x, res = 100))

    tree_list_amalg <- paramo.list(names(tree_list), tree_list, ntrees = 2)

    # Check number of trees.
    expect_true(all(sapply(tree_list, length) == length(tree_list_amalg)))

    # Check number of branches.
    n_br1 <- unlist(lapply(tree_list, function(x) lapply(x, function(y) length(y$maps) ) ))
    n_br2 <- sapply(tree_list_amalg, function(x) length(x$maps) )
    expect_true(all(n_br1 == n_br2))

    # Check number of amalgamated characters.
    n_chars <- unique(unlist(lapply(tree_list_amalg, function(x) sapply(x$maps, function(y) nchar(names(y)) ) )))
    expect_true(length(tree_list) == n_chars)

  }
)
