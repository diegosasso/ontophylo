test_that(
  "Testing get_vector_ids_per_term()",
  {

    data("HAO")

    graph_table <- readRDS("../testdata/graph_tb_test.RDS")

    term1 <- "HAO:0000397"
    term2 <- "HAO:0000494"

    v1 <- get_vector_ids_per_term(term = term1, ONT = HAO, GR = graph_table)
    v2 <- get_vector_ids_per_term(term = term2, ONT = HAO, GR = graph_table)

    pics1 <- paste0(graph_table$pic_id[1:2], collapse = ", ")
    pics1 <- as.numeric(stringr::str_split(pics1, pattern = ", ")[[1]])
    pics2 <- paste0(graph_table$pic_id[4], collapse = ", ")
    pics2 <- as.numeric(stringr::str_split(pics2, pattern = ", ")[[1]])

    # Check queries of picture layers.
    expect_identical(v1, pics1)
    expect_identical(v2, pics2)

  }
)
