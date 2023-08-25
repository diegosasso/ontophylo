test_that(
  "Testing edge_profiles4plotting()",
  {

    data("hym_tree")

    tree_discr <- discr_Simmap(hym_tree, res = 4000)

    Edge_KDE <- readRDS("../testdata/kde_test.RDS")

    Edge_KDE_stat <- Edge_KDE$loess.lambda.mean

    edge_prof <- edge_profiles4plotting(tree_discr, Edge_KDE_stat)

    # Check lengths.
    expect_identical(length(unique(edge_prof$edge.id[edge_prof$edge.type == "main"])), length(tree_discr$edge.length))

  }
)
