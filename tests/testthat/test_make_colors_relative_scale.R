test_that(
  "Testing make_colors_relative_scale()",
  {

    stat <- setNames(runif(5, 0.1, 10), c("cranium", "fore_wing", "hind_wing", "pronotum", "propectus") )

    hm.palette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"), space = "Lab")

    cols_maps <- make_colors_relative_scale(stat, palette = hm.palette(100), lims = c(min(stat), max(stat)))

    # Check class.
    expect_true(is.character(cols_maps))

    # Check for unique colors.
    expect_true(anyDuplicated(cols_maps) == 0)

  }
)
