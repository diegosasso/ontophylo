test_that(
  "Testing MultiScale.simmap()",
  {

    stm_amalg <- readRDS("../testdata/stm_amalg_test.RDS")

    tree_mds <- stm_amalg[[1]]

    MD <- suppressWarnings(MultiScale.simmap(tree_mds))

    # Check list elements.
    expect_identical(names(MD), c("Points", "Lines", "Edge.map"))

    # Check number of bins.
    expect_true(sum(sapply(tree_mds$maps, length)) == dim(MD$Points)[1])

    # Check tree age.
    expect_true(round(max(MD$Points$time), 4) == round(max(phytools::nodeHeights(tree_mds)), 4))

    # Check number of tips.
    expect_true(sum(MD$Points$sp_extant == "yes") == length(tree_mds$tip.label))

    # Check root.
    expect_true(sum(MD$Points$is_root == "yes") == 2)

  }
)
