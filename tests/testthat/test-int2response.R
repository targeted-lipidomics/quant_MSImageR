require(testthat)
require(quantMSImageR)

context("test intensity values normalised to internal standard")

test_that("int2response function", {

  # Create test data
  fdata <- MassDataFrame(mz=c(500, 510, 540, 550), analyte = c("IS", rep("analyte", 3)))
  pdata <- PositionDataFrame(run=c(rep("run1", 2), rep("run2", 2)),
                             coord=expand.grid(x=1:2, y=1:2),
                             sample_ID = c(rep("run1", 2), rep("run2", 2)))

  ints <- matrix(nrow=4, ncol=4,
                 data = c(rep(c(5,10,15, 100), 2),
                          rep(c(10,10,15, 100), 2)))

  test_data <- MSImagingExperiment(imageData=ints,
                                   featureData=fdata,
                                   pixelData=pdata)

  test_data <- as(test_data, "quant_MSImagingExperiment")
  new_data <- int2response(test_data)

  expect_equal(ncol(test_data), ncol(new_data))
  expect_equal(nrow(test_data), (nrow(new_data)+1))
  expect_true(all(spectra(new_data)[1,] == c(2,2,1.0,1.0)))
  expect_true(all(spectra(new_data)[2,] == c(3,3,1.5,1.5)))
  expect_true(all(spectra(new_data)[3,] == c(20,20,10.0,10.0)))

})
