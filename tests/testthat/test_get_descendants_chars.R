test_that(
  "Testing get_descendants_chars()",
  {

    data("HAO")

    onto <- HAO

    search_term <- "HAO:0000653"

    annot_list <- list("CH1" = "HAO:0000653", "CH2" = "HAO:0000653")

    onto$terms_selected_id <- annot_list

    desc_chars <- get_descendants_chars(onto, annotations = "manual", "HAO:0000653")

    # Check class.
    expect_true(is.character(desc_chars))
    # Check names.
    expect_identical(names(annot_list), desc_chars)
    # Check lengths.
    expect_identical(length(annot_list), length(desc_chars))

  }
)
