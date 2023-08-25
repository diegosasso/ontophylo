test_that(
  "Testing add_pseudodata()",
  {

    ph <- readRDS("../testdata/ph_test.RDS")

    nhpp <- make_data_NHPP_KDE_Markov_kernel(ph)

    psd <- lapply(nhpp, function(x) -x )

    edge_groups <- as.list(unique(ph$Focal.Edge.id))

    nhpp_psd <- add_pseudodata(Edge.groups = edge_groups, Pseudo.data = psd, Path.data = nhpp)

    # Check number of branches.
    expect_identical(length(nhpp), length(nhpp_psd))

    # Check number of bins.
    expect_true(all.equal(sapply(psd, function(x) length(x)*2), sapply(nhpp_psd, function(x) length(x))))

  }
)
