test_that(
  "Testing get_state()",
  {

    vec <- setNames(sort(runif(1001,0,1)), 0:1000)
	
	x1 <- vec[1]
	x2 <- vec[vec == median(vec)]
	x3 <- vec[1001]
	y1 <- ontophylo:::get_state(vec, x1)
	y2 <- ontophylo:::get_state(vec, x2)
	y3 <- ontophylo:::get_state(vec, x3)

    # Check for correct states.
	expect_true(y1 == "0")
	expect_true(y2 == "499")
	expect_true(y3 == "1000")

  }
)
