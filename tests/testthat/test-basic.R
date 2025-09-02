test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("basic check", {
  expect_true(TRUE)
})

test_that("app_server initializes without errors", {
  shiny::testServer(app_server, {
    # inside here, we have a working session, input, output
    expect_true(TRUE)  # if we got here, app_server didn't crash
  })
})
