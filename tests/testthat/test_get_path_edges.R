test_that(
  "Testing get_path_edges()",
  {

    stm <- readRDS("../testdata/stm_amalg_test.RDS")
    stm <- stm[[1]]

    node = 3

    st_ph <- suppressWarnings(ontophylo:::get_path_edges(stm, node))

    # Check path edges.
    expect_identical(head(st_ph, n = 1), which(stm$edge[,2] == node))
    expect_true(tail(st_ph, n = 1) %in% which(stm$edge[,1] == length(stm$tip.label)+1))

  }
)
