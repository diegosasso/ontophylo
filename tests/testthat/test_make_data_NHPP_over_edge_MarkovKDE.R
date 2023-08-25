test_that(
  "Testing make_data_NHPP_over_edge_MarkovKDE()",
  {

    ph <- readRDS("../testdata/ph_test.RDS")

    edge = 10

    nhpp_edge <- ontophylo:::make_data_NHPP_over_edge_MarkovKDE(ph, Focal.Edge = edge)

    x <- ph$t.start[ph$Focal.Edge.id == edge]
    x <- x[-which(x == 0)]

    # Check number of bins.
    expect_identical(length(nhpp_edge), length(x))

  }
)
