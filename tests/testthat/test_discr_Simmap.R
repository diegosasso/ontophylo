test_that(
  "Testing discr_Simmap()",
  {

    stm <- readRDS("../testdata/stm_test.RDS")
    stm <- stm[[1]]

    stm_discr <- discr_Simmap(stm, res = 100)

    edges1 <- sapply(stm$maps, function(x) round(sum(x),4) )
    edges2 <- sapply(stm_discr$maps, function(x) round(sum(x),4) )

    # Check output.
    expect_identical(edges1, edges2)


  }
)
