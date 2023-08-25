test_that(
  "Testing merge_branch_cat()",
  {

    stm <- readRDS("../testdata/stm_test.RDS")
    stm <- stm[[1]]

    stm_discr <- discr_Simmap(stm, res = 100)

    edge1 <- stm_discr$maps[[50]]
    edge2 <- merge_branch_cat(edge1)

    # Check output.
    expect_identical(round(sum(edge1), 4), round(unname(edge2), 4))

  }
)
