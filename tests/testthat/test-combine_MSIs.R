require(testthat)
require(quantMSImageR)

context("test MSImagingExperiment objects are combined correctly")

test_that("combine_MSIs function", {

  # Create test data
  fdata <- MassDataFrame(mz=1:4,
                         analyte = c("IS", rep("analyte", 3)),
                         precursor_mz = c(500, 510, 540, 550),
                         product_mz = c(500, 510, 540, 550)/2,
                         name = sprintf("lipid_%s", 1:4))

  pdata1 <- PositionDataFrame(run=c(rep("run1", 2), rep("run2", 2)),
                             coord=expand.grid(x=1:2, y=1:2),
                             sample_ID = c(rep("run1", 2), rep("run2", 2)))

  pdata2 <- PositionDataFrame(run=c(rep("run3", 2), rep("run4", 2)),
                              coord=expand.grid(x=1:2, y=1:2),
                              sample_ID = c(rep("run3", 2), rep("run4", 2)))

  ints1 <- matrix(nrow=4, ncol=4,
                 data = c(rep(c(5,10,15, 100), 2),
                          rep(c(10,10,15, 100), 2)))

  ints2 <- matrix(nrow=4, ncol=4,
                  data = c(rep(c(3,6,9, 30), 2),
                           rep(c(6,6,9, 30), 2)))

  test_data1 <- MSImagingExperiment(spectraData=ints1,
                                   featureData=fdata,
                                   pixelData=pdata1)

  test_data2 <- MSImagingExperiment(spectraData=ints2,
                                    featureData=fdata,
                                    pixelData=pdata2)

  new_data <- combine_MSIs(test_data1, test_data2)

  expect_equal((ncol(test_data1) + ncol(test_data1)), ncol(new_data))
  expect_equal(nrow(test_data1), nrow(new_data))

  expect_true(all(spectra(new_data)[1,] == c(5, 5, 10, 10, 3, 3, 6, 6)))
  expect_true(all(spectra(new_data)[2,] == c(10, 10, 10, 10, 6, 6, 6, 6)))
  expect_true(all(spectra(new_data)[3,] == c(15, 15, 15, 15, 9, 9, 9, 9)))
  expect_true(all(spectra(new_data)[4,] == c(100, 100, 100, 100, 30, 30, 30, 30)))

})
