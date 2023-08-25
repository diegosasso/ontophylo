test_that(
  "Testing get_rough_state_cols()",
  {

    stm <- readRDS("../testdata/stm_test.RDS")

    cols <- get_rough_state_cols(stm[[1]])

    # Check class.
    expect_true(is.character(cols))

    # Check number of colors.
    expect_true(length(cols) == length(unique(unlist(lapply(stm[[1]]$maps, function(x) names(x) )))))

  }
)
