library(chords)
context("testing estimators")

test_that("various estimators", {
  data(brazil)
  expect_message(rds.object<- initializeRdsObject(brazil))
  expect_equal(length(rds.object$I.t), 303)
  rds.object2 <- Estimate.b.k(rds.object = rds.object )
  expect_equal(length(rds.object2$estimates),10)

  # Integrated
  rds.object3 <- Estimate.b.k(rds.object = rds.object, type = 'integrated')
  expect_equal(length(rds.object3$estimates),10)
  expect_equal(sum(rds.object3$estimates$Nk.estimates), 313)
  
  # jack.control <- makeJackControl(1,1e1)
  # rds.object4 <- Estimate.b.k(rds.object = rds.object, type = 'leave-d-out', jack.control = jack.control)
  # expect_equal(length(rds.object4$estimates),10)
  # expect_equal(sum(rds.object4$estimates$Nk.estimates, na.rm=TRUE), 1330)
  
  
  expect_equal(length(thetaSmoothingNks(rds.object2)),28)
  rds.object5 <- Estimate.b.k(rds.object = rds.object2, type = 'parametric')
  expect_equal(length(rds.object5$estimates),10)
  expect_equal(sum(rds.object5$estimates$Nk.estimates<Inf), 997)
  
  rds.object.6 <- Estimate.b.k(rds.object = rds.object2, type='rescaling')
  rds.object.5$estimates$Nk.estimates
  see(rds.object.5) # will not work. less than 2 converging estimates.
  
  
})
