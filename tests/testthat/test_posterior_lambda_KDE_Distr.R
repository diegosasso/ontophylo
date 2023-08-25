test_that(
  "Testing posterior_lambda_KDE_Distr()",
  {

    stm <- readRDS("../testdata/stm_test.RDS")

    n_sim = 10
    p_lamb_dist <- posterior_lambda_KDE_Distr(stm, n.sim = n_sim, "BR.name")

    # Check outputs.
    expect_true(all(dim(p_lamb_dist) == c(n_sim,2)))
    expect_true(is.numeric(p_lamb_dist$sim))

  }
)
