require(testthat)
require(quantMSImageR)

context("test calibration curves are created")

test_that("create_cal_curve function", {

  # Create test data
  test_data = as(readMSIData(file = sprintf("%s/tissue_MRM_data.raw/combined.imzML", system.file('extdata', package = 'quantMSImageR'))), "quant_MSImagingExperiment")
  test_data@calibrationInfo@cal_metadata = cal_metadata = read.csv(sprintf("%s/calibration_metadata.csv", system.file('extdata', package = 'quantMSImageR')))
  test_data@calibrationInfo@cal_response_data = read.csv("C:/Users/matsmi/OneDrive - Karolinska Institutet/Dokument/MSI/quantMSImageR/inst/extdata/cal_response_data.csv")

  # Test cal
  new_data <- create_cal_curve(test_data, cal_type = "Cal")

  coeffs = as.numeric(new_data@calibrationInfo@cal_list[[1]]$coefficients)

  expect_equal(signif(coeffs[1], 4), 4414)
  expect_equal(signif(coeffs[2], 6), 339782)
  expect_equal(signif(new_data@calibrationInfo@r2_df$r2[1], 5), 0.92667)

})
