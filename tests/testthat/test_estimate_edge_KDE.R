test_that(
  "Testing estimate_edge_KDE()",
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

    # Check class.
    expect_true(is.list(Edge_KDE))

    # Check names.
    expect_identical(names(Edge_KDE), c("Maps.mean", "Maps.mean.norm"))

    # Check lengths.
    expect_true(all(sapply(Edge_KDE, length) == length(stm$maps)))
    expect_true(all(sapply(Edge_KDE, length) == length(nhpp_psd)))

  }
)
