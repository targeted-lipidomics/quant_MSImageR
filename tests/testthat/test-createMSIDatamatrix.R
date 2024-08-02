require(testthat)
require(quantMSImageR)

context("test data matrix is formed from MSI object")

test_that("createMSIDatamatrix function", {

  # Create test data
  test_data = as(readMSIData(file = sprintf("%s/tissue_MRM_data.raw/combined.imzML", system.file('extdata', package = 'quantMSImageR'))),  "quant_MSImagingExperiment")

  new_data = createMSIDatamatrix(test_data, val_slot = "intensity", roi_header = "identifier")

  roi_average_matrix = new_data@tissueInfo@roi_average_matrix

  expect_equal(nrow(roi_average_matrix), 34)
  expect_equal(signif(roi_average_matrix$`12_13-DiHOME`[3], 7), 4006.476)
  expect_equal(signif(roi_average_matrix$`12_13-DiHOME`[21], 7), 41279.94)
  expect_equal(signif(roi_average_matrix$`12_13-DiHOME`[33], 7), 2863.862)

  all_pixel_matrix = new_data@tissueInfo@all_pixel_matrix

  expect_equal(nrow(all_pixel_matrix), 21506)
  expect_equal(signif(all_pixel_matrix$`12_13-DiHOME`[3010], 7), 11717)
  expect_equal(signif(all_pixel_matrix$`12_13-DiHOME`[21500], 7), 2131)
  expect_equal(signif(all_pixel_matrix$`12_13-DiHOME`[10579], 7), 6515)

})
