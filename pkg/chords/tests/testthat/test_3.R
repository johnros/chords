library(chords)
context("testing simulator")

test_that("simulator", {
  expect_is(makeSparseMatrix(2, 2e3, 1)!=0, 'Matrix')
})
