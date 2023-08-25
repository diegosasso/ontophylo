test_that(
  "Testing test_make_postPois_KDE()",
  {

    stm <- readRDS("../testdata/stm_kde_test.RDS")

    Edge_KDE <- readRDS("../testdata/kde_test.RDS")

    Edge_KDE_stat <- Edge_KDE$Maps.mean.loess.norm

    lambda_post <- posterior_lambda_KDE(stm)

    lambda_mean <- make_postPois_KDE(Edge_KDE_stat, lambda_post, lambda.post.stat = "Mean")

    # Check lengths.
    expect_identical(length(Edge_KDE_stat), length(lambda_mean))
    expect_identical(sapply(Edge_KDE_stat, length), sapply(lambda_mean, length))

  }
)
