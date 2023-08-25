test_that(
  "Testing stack_stm()",
  {

    stm <- readRDS("../testdata/stm_test.RDS")
    stm <- stm[1:5]

    stm_discr <- discr_Simmap_all(stm, res = 100)

    stm_amalg <- ontophylo:::stack_stm(stm_discr)

    ind_maps_list <- lapply(stm_discr, function(x) x$maps)
    ind_maps_list <- lapply(ind_maps_list, function(x) lapply(x, function(y) names(y)))

    amalg_maps_list <- lapply(stm_amalg$maps, function(x) names(x))

    # Check number of amalgamated characters.
    expect_true(all(length(ind_maps_list) == unlist(lapply(amalg_maps_list, nchar))))

    # Check number of branches.
    expect_true(all(sapply(ind_maps_list, function(x) length(x)) == length(amalg_maps_list)))

    # Check number of bins on branches.
    expect_identical(sapply(ind_maps_list[[1]], length), sapply(amalg_maps_list, length))

  }
)
