test_that(
  "Testing add_noise_MD()",
  {

    stm_amalg <- readRDS("../testdata/stm_amalg_test.RDS")

    tree_mds <- stm_amalg[[1]]

    MD <- suppressWarnings(MultiScale.simmap(tree_mds, add.noise = NULL))

    noise = runif(2, 0.01, 1)

    MD_n <- add_noise_MD(MD, noise)

    # Check amount of noise added to points.
    expect_true(all.equal((mean(abs(MD_n$Points$V1 - MD$Points$V1))*2), noise[1], tolerance = 0.1))
    expect_true(all.equal((mean(abs(MD_n$Points$V2 - MD$Points$V2))*2), noise[2], tolerance = 0.1))

  }
)
