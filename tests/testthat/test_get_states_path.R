test_that(
  "Testing get_states_path()",
  {

    stm <- readRDS("../testdata/stm_amalg_test.RDS")
    stm <- stm[[1]]

    node = 40

    st_ph <- suppressWarnings(ontophylo:::get_states_path(stm, node))

    # Check length of paths.
    expect_true(all.equal(sum(stm$edge.length[unique(st_ph$Edge.id)]), sum(st_ph$delta.t)))

  }
)
