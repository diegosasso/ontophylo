test_that(
  "Testing path_hamming_over_all_edges()",
  {

    stm <- readRDS("../testdata/stm_amalg_test.RDS")
    stm <- stm[[1]]

    ph <- suppressWarnings(path_hamming_over_all_edges(stm))

    ph_st <- unique(ph$States)

    fun_states <- function(x) {unique(unlist(lapply(x$maps, function(y) unique(names(y)) )))}
    states <- fun_states(stm)

    # Check if all states where sampled.
    expect_identical(ph_st, states)

  }
)
