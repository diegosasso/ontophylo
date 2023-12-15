test_that(
  "Testing make_contMap_KDE()",
  {

    hym_tree <- readRDS("../testdata/hym_tree.RDS")

    tree_discr <- discr_Simmap(hym_tree, res = 4000)

    Edge_KDE <- readRDS("../testdata/kde_test.RDS")

    Edge_KDE_stat <- Edge_KDE$loess.lambda.mean

    nhpp_map <- make_contMap_KDE(tree_discr, Edge_KDE_stat)

    # Check class.
    expect_true(class(nhpp_map) == "contMap")

    # Check mapping.
    expect_identical(sapply(nhpp_map$tree$maps, length), sapply(Edge_KDE_stat, length))

  }
)
