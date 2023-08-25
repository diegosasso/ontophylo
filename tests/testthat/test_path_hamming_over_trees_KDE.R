test_that(
  "Testing path_hamming_over_trees_KDE()",
  {

    stm <- readRDS("../testdata/stm_amalg_test.RDS")

    ph <- suppressWarnings(path_hamming_over_trees_KDE(stm))

    ph_st <- unique(ph$States)

    fun_states <- function(x) {unique(unlist(lapply(x$maps, function(y) unique(names(y)) )))}
    states <- unique(unlist(lapply(stm, fun_states)))

    # Check number of trees.
    expect_identical(length(unique(ph$tree.id)), length(stm))

    # Check if all states where sampled.
    expect_identical(ph_st, states)

  }
)
