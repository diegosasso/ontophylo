test_that(
  "Testing estimate_band_W()",
  {

    stm <- readRDS("../testdata/stm_amalg_test.RDS")
    stm <- stm[[1]]

    ph <- readRDS("../testdata/ph_test.RDS")

    nhpp <- make_data_NHPP_KDE_Markov_kernel(ph)

    psd <- lapply(nhpp, function(x) -x )

    edge_groups <- as.list(unique(ph$Focal.Edge.id))

    nhpp_psd <- add_pseudodata(Edge.groups = edge_groups, Pseudo.data = psd, Path.data = nhpp)

    bdw1 <- estimate_band_W(stm, nhpp_psd, band.width = c('bw.nrd'))
    bdw2 <- estimate_band_W(stm, nhpp_psd, band.width = c('bw.nrd0'))
    bdw3 <- estimate_band_W(stm, nhpp_psd, band.width = c('bw.ucv'))
    bdw4 <- estimate_band_W(stm, nhpp_psd, band.width = c('bw.SJ'))

    # Check output lengths.
    expect_identical(length(bdw1), length(stm$tip.label))
    expect_identical(length(bdw2), length(stm$tip.label))
    expect_identical(length(bdw3), length(stm$tip.label))
    expect_identical(length(bdw4), length(stm$tip.label))

  }
)
