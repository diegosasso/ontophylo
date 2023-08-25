test_that(
  "Testing derivative()",
  {

    stm <- readRDS("../testdata/stm_kde_test.RDS")

    Edge_KDE <- readRDS("../testdata/kde_test.RDS")

    Edge_KDE_stat <- Edge_KDE$loess.lambda.mean

    deriv_KDE <- suppressWarnings(derivative_KDE(stm[[1]], Edge_KDE_stat))

    # check class.
    expect_true(is.numeric(unlist(deriv_KDE)))

    # Check for nulls.
    expect_false(any(is.null(unlist(deriv_KDE))))

    # Check lengths.
    expect_identical(length(Edge_KDE_stat), length(deriv_KDE))

  }
)
