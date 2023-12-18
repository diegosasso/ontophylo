test_that(
  "Testing path_hamming()",
  {

    stm <- readRDS("../testdata/stm_amalg_test.RDS")
    stm <- stm[[1]]

    node = 3

    st_ph <- suppressWarnings(ontophylo:::get_states_path(stm, node))

    path_hm <- suppressWarnings(ontophylo:::path_hamming(st_ph))

    # Check output columns.
    expect_true(is.numeric(path_hm$Ham.dist))
    expect_true(is.numeric(path_hm$Ham.dist.n))

  }
)
