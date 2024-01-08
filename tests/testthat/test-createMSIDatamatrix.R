require(testthat)
require(quantMSImageR)

context("test data matrix is formed from MSI object")

test_that("createMSIDatamatrix function", {

  # Create test data
  fdata <- MassDataFrame(mz=c(500, 510, 540, 550))

  pdata <- PositionDataFrame(run= "run1",
                             coord=expand.grid(x=1:3, y=1:3),
                             replicate = "01",
                             sample_type = "Tissue",
                             ROI = c(rep("roi1", 3), rep("roi2", 3), rep("roi3", 2), NA))
  pdata$sample_ID = sprintf("%s_rep%s_%s", pdata$sample_type, pdata$replicate, pdata$ROI)

  ints <- matrix(nrow=4, ncol=9, byrow=T,
                 data = c(rep(c(0,0,0, 1500), 1),
                          rep(c(0,100,105, 1500), 4),
                          rep(c(60,105,150, 1400), 4)))

  test_data <- MSImagingExperiment(imageData=ints,
                                   featureData=fdata,
                                   pixelData=pdata)

  test_data <- as(test_data, "quant_MSImagingExperiment")

  new_data = createMSIDatamatrix(test_data)

  spectra(new_data)
  conc_matrix = new_data@tissueInfo@conc_matrix

  expect_true(is.na(conc_matrix[1,1]))
  expect_equal(conc_matrix[3,3], 82.5)
  expect_equal(conc_matrix[2,3], 750)
  expect_equal(ncol(conc_matrix), 3)

  expect_equal(as.numeric(nrow(new_data)), nrow(conc_matrix))
  expect_equal(as.numeric(ncol(new_data)), 8)

})
