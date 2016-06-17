library(chords)
context("Utility functions")

test_that("RDS example", {
  expect_equal(length(rdsObjectConstructor()),6)
  expect_equal(length(makeRDSExample()), 6)
})
