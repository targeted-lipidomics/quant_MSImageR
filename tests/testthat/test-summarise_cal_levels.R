require(testthat)
require(quantMSImageR)

context("test calibration levels are summarised correctly")

test_that("summarise_cal_levels function", {

  # Create test data
  fdata <- MassDataFrame(mz=c(500, 510, 540, 550))
  pdata <- PositionDataFrame(run="run1",
                             coord=expand.grid(x=1:4, y=1:4),
                             sample_type = "Cal",
                             replicate = "01",
                             ROI = c(rep(NA, 2), rep("01", 2),
                                     rep(NA, 2), rep("02", 2),
                                     rep(NA, 2), rep("03", 2),
                                     rep(NA, 2), rep("04", 2)))
  pdata$sample_ID = sprintf("%s_rep%s_%s", pdata$sample_type, pdata$replicate, pdata$ROI)

  ints <- matrix(nrow=4, ncol=16, byrow = T,
                 data = rep(c(rep(1, 4), rep(10, 4),
                          rep(100, 4), rep(1000, 4)), 4))

  test_data <- MSImagingExperiment(imageData=ints,
                                   featureData=fdata,
                                   pixelData=pdata)

  test_data <- as(test_data, "quant_MSImagingExperiment")
  new_data <- summarise_cal_levels(test_data)

  expect_equal(as.numeric(nrow(new_data)), 4)
  expect_equal(as.numeric(ncol(new_data)), 16)
  expect_equal(new_data@calibrationInfo@pixels_per_level[[1]], 2)
  expect_equal(new_data@calibrationInfo@pixels_per_level[[4]], 2)

  int_vec = c(0.5, 5, 50, 500)

  expect_true(all(new_data@calibrationInfo@response_per_pixel[1,] == int_vec))
  expect_true(all(new_data@calibrationInfo@response_per_pixel[4,] == int_vec))

})
