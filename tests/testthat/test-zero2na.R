require(testthat)
require(quantMSImageR)

context("test zero values changed to NA")

test_that("zero2na function", {

  # Create test data
  fdata <- MassDataFrame(mz=c(500, 510, 540, 550))
  pdata <- PositionDataFrame(run="run0", coord=expand.grid(x=1:2, y=1:2))
  ints <- matrix(nrow=4, ncol=4, rep(c(0,10,3, 7), 4))

  test_data <- MSImagingExperiment(imageData=ints,
                             featureData=fdata,
                             pixelData=pdata)

  test_data <- as(test_data, "quant_MSImagingExperiment")
  new_data <- zero2na(test_data)

  expect_equal(dim(test_data), dim(new_data))
  expect_true(all(is.na(spectra(new_data)[1,])))
  expect_false(any(is.na(spectra(new_data)[2:4,])))

})
