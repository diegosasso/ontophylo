test_that(
  "Testing make_data_NHPP_KDE_Markov_kernel()",
  {

    ph <- readRDS("../testdata/ph_test.RDS")

    nhpp <- make_data_NHPP_KDE_Markov_kernel(ph, add.psd = FALSE)

    x <- setNames(ph$t.start, ph$Focal.Edge.id)
    x <- x[-which(x == 0)]
    x <- split(x, f = names(x))
    x <- sapply(x[paste(1:length(x))], length)

    # Check number of branches.
    expect_identical(length(nhpp), length(unique(ph$Edge.id)))

    # Check number of bins.
    expect_identical(sapply(nhpp, length), unname(x))

  }
)
