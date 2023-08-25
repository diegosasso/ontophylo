test_that(
  "Testing posterior_lambda_KDE()",
  {

    stm <- readRDS("../testdata/stm_test.RDS")

    lambda_post <- posterior_lambda_KDE(stm)

    # Check class.
    expect_true(is.numeric(lambda_post$Mean))
    expect_true(is.numeric(lambda_post$SD))
    expect_true(is.numeric(lambda_post$Q_2.5))
    expect_true(is.numeric(lambda_post$Q_97.5))

    # Check output.
    expect_true(lambda_post$Mean != 0)
    expect_true(lambda_post$SD != 0)

  }
)
