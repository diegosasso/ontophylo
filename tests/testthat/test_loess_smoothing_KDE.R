test_that(
  "Testing loess_smoothing_KDE()",
  {

    stm <- readRDS("../testdata/stm_amalg_test.RDS")
    stm <- stm[[1]]

    ph <- readRDS("../testdata/ph_test.RDS")

    nhpp <- make_data_NHPP_KDE_Markov_kernel(ph)

    psd <- lapply(nhpp, function(x) -x )

    edge_groups <- as.list(unique(ph$Focal.Edge.id))

    nhpp_psd <- add_pseudodata(Edge.groups = edge_groups, Pseudo.data = psd, Path.data = nhpp)

    bdw <- estimate_band_W(stm, nhpp_psd, band.width = c('bw.nrd'))

    Edge_KDE <- estimate_edge_KDE_Markov_kernel_unnorm(stm, nhpp_psd, h = bdw)

    Edge_KDE$Maps.mean.loess <- suppressWarnings(loess_smoothing_KDE(stm, Edge_KDE))

    # Check lengths.
    expect_identical(length(Edge_KDE$Maps.mean), length(Edge_KDE$Maps.mean.loess))

  }
)
