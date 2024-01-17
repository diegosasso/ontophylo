test_that(
  "Testing make_colors()",
  {

    stat <- setNames(c(0.1, 0.2, 0.3, 0.4, 0.5), c("cranium", "fore_wing", "hind_wing", "pronotum", "propectus") )

    hm.palette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"), space = "Lab")

    cols_maps <- make_colors(stat, palette = hm.palette(100))

    # Check class.
    expect_true(is.character(cols_maps))

    # Check for unique colors.
    expect_true(length(stat) == length(cols_maps))

  }
)
