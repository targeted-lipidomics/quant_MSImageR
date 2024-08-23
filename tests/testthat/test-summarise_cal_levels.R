require(testthat)
require(quantMSImageR)

context("test calibration levels are summarised correctly")

test_that("summarise_cal_levels function", {

  # Create test data
  test_data = readMSIData(file = sprintf("%s/tissue_MRM_data.raw/combined.imzML", system.file('extdata', package = 'quantMSImageR')))
  cal_metadata = read.csv(sprintf("%s/calibration_metadata.csv", system.file('extdata', package = 'quantMSImageR')))

  test_data <- as(test_data, "quant_MSImagingExperiment")

  new_data <- summarise_cal_levels(MSIobject = test_data,
                                   cal_metadata = cal_metadata,
                                   val_slot = "intensity",
                                   cal_label = "Cal",
                                   id = "identifier")


  output = new_data@calibrationInfo@cal_response_data

  expect_equal(nrow(output), 21)
  expect_equal(ncol(output), 7)

  expect_equal(signif(output$response_perpixel[1],6), 38000.6)
  expect_equal(signif(output$response_perpixel[6],6), 5172.83)
  expect_equal(signif(output$response_perpixel[11],6), 40830.7)
  expect_equal(signif(output$response_perpixel[16],6), 9177.55)
  expect_equal(signif(output$response_perpixel[21],6), 4975.11)

})
