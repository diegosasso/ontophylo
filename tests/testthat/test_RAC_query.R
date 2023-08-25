test_that(
  "Testing RAC_query()",
  {

    data("HAO")

    onto <- HAO

    char_info <- cbind(c("CH1", "CH2", "CH3", "CH4", "CH5", "CH6"),
                       c("HAO:0000234", "HAO:0000234", "HAO:0000853",
                         "HAO:0000854", "HAO:0000053", "HAO:0000626"))

    char_info <- data.frame(char_info)

    query_terms <- c("head", "mesosoma", "metasoma")
    query_wrong <- c("AAAA")

    query <- RAC_query(char_info, onto, query_terms)

    # Check class.
    expect_true(is.list(query))
    # Check names.
    expect_identical(query_terms, names(query))
    # Check query correct output.
    expect_true(all(char_info[[1]] %in% unname(unlist(query))))
    # Check wrong input.
    expect_error(RAC_query(char_info, onto, query_wrong))

  }
)
