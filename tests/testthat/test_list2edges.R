test_that(
  "Testing list2edges()",
  {

    annot_list <- list("CH1" = c("HAO:0000933", "HAO:0000958"), "CH2" = c("HAO:0000833", "HAO:0000258"))

    res <- list2edges(annot_list)

    # Check class.
    expect_true(is.matrix(res))
    # Check names.
    expect_true(all(names(annot_list) %in% res[,1]))
    # Check terms.
    expect_identical(unname(unlist(annot_list)), res[,2])

  }
)
