test_that(
  "Testing merge_tree_cat_list()",
  {

    stm <- readRDS("../testdata/stm_test.RDS")

    stm_discr <- discr_Simmap_all(stm, res = 100)

    stm_merg <- merge_tree_cat_list(stm_discr)

    edge1 <- stm_discr[[1]]$maps[[50]]
    edge2 <- stm_merg[[1]]$maps[[50]]

    # Check output.
    expect_identical(round(sum(edge1), 4), round(unname(edge2), 4))

  }
)
