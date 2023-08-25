test_that(
  "Testing integrate_edge_KDE()",
  {

    stm <- readRDS("../testdata/stm_amalg_test.RDS")
    stm <- stm[[1]]

    ph <- readRDS("../testdata/ph_test.RDS")

    nhpp <- make_data_NHPP_KDE_Markov_kernel(ph)

    psd <- lapply(nhpp, function(x) -x )

    edge_groups <- as.list(unique(ph$Focal.Edge.id))

    nhpp_psd <- add_pseudodata(Edge.groups = edge_groups, Pseudo.data = psd, Path.data = nhpp)

    bdw <- estimate_band_W(stm, nhpp_psd, band.width = c('bw.nrd'))

    Edge_KDE <- estimate_edge_KDE(stm, nhpp_psd, h = bdw)
    Edge_KDE$Maps.mean.loess <- suppressWarnings(loess_smoothing_KDE(stm, Edge_KDE))
    Edge_KDE$Maps.mean.loess.norm <- normalize_KDE(stm, Edge_KDE$Maps.mean.loess)

    v1 <- integrate_edge_KDE(stm, Edge_KDE$Maps.mean.norm)
    v2 <- integrate_edge_KDE(stm, Edge_KDE$Maps.mean.loess.norm)

    # Check integrals.
    expect_true(all.equal(as.numeric(v1), as.numeric(1)))
    expect_true(all.equal(as.numeric(v2), as.numeric(1)))

  }
)
