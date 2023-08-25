test_that(
  "Testing discr_Simmap_all()",
  {

    stm <- readRDS("../testdata/stm_test.RDS")

    stm_discr <- discr_Simmap_all(stm, res = 100)

    edges1 <- sapply(stm[[1]]$maps, function(x) round(sum(x),4) )
    edges2 <- sapply(stm_discr[[1]]$maps, function(x) round(sum(x),4) )

    # Check output.
    expect_identical(edges1, edges2)

  }
)
