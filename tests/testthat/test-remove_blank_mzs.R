require(testthat)
require(quantMSImageR)

context("test m/z values where all data is empty are removed")

test_that("remove_blank_mzs function", {

  # Create test data
  fdata <- MassDataFrame(mz=c(500, 510, 540, 550))
  pdata <- PositionDataFrame(run="run1", coord=expand.grid(x=1:2, y=1:2))
  ints <- matrix(nrow=4, ncol=4, byrow=T,
                 c(rep(NA, 4), rep(c(0,10,3, 7), 3)))

  test_data <- MSImagingExperiment(spectraData=ints,
                                   featureData=fdata,
                                   pixelData=pdata)

  test_data <- as(test_data, "quant_MSImagingExperiment")
  new_data <- remove_blank_mzs(test_data)

  expect_equal(as.numeric(nrow(new_data)), 3)
  expect_equal(as.numeric(ncol(new_data)), 4)
})
