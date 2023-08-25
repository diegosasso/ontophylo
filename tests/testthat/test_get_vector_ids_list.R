test_that(
  "Testing get_vector_ids_list()",
  {

    data("HAO")

    graph_table <- readRDS("../testdata/graph_tb_test.RDS")

    terms_list <- as.list(c("HAO:0000397", "HAO:0000494"))
    terms_list <- setNames(terms_list, c("head", "leg"))

    pic_layers <- get_vector_ids_list(terms = terms_list , ONT = HAO, GR = graph_table)

    pics1 <- paste0(graph_table$pic_id[1:2], collapse = ", ")
    pics1 <- as.numeric(stringr::str_split(pics1, pattern = ", ")[[1]])
    pics2 <- paste0(graph_table$pic_id[4], collapse = ", ")
    pics2 <- as.numeric(stringr::str_split(pics2, pattern = ", ")[[1]])

    v1 <- unname(pic_layers[names(pic_layers) %in% names(terms_list)[1]])
    v2 <- unname(pic_layers[names(pic_layers) %in% names(terms_list)[2]])

    # Check queries of picture layers.
    expect_identical(v1, pics1)
    expect_identical(v2, pics2)

  }
)
