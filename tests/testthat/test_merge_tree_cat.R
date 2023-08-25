test_that(
  "Testing merge_tree_cat()",
  {

    stm <- readRDS("../testdata/stm_test.RDS")
    stm <- stm[[1]]

    stm_discr <- discr_Simmap(stm, res = 100)

    stm_merg <- merge_tree_cat(stm_discr)

    edge1 <- stm_discr$maps[[50]]
    edge2 <- stm_merg$maps[[50]]

    # Check output.
    expect_identical(round(sum(edge1), 4), round(unname(edge2), 4))

  }
)
