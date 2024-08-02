require(testthat)
require(quantMSImageR)

context("test intensity values normalised to internal standard")

test_that("int2response function", {

  # Create test data
  fdata <- MassDataFrame(mz=c(500, 510, 540, 550), analyte = c("IS", rep("analyte", 3)))

  pdata <- PositionDataFrame(run=c(rep("run1", 2), rep("run2", 2), rep("Noise", 2)),
                             coord=expand.grid(x=1:3, y=1:2),
                             sample_ID = c(rep("Tissue", 4), rep("Noise", 2)))

  ints <- matrix(nrow=4, ncol=6,
                 data = c(rep(c(5,100,150, 100), 2),
                          rep(c(10,100,105, 300), 2),
                          1, 2, 3, 4,
                          3,6, 9, 12))

  test_data <- MSImagingExperiment(spectraData=ints,
                                   featureData=fdata,
                                   pixelData=pdata)

  test_data <- as(test_data, "quant_MSImagingExperiment")

  new_data <- int2snr(MSIobject = test_data, val_slot = "intensity",
                      noise = "Noise", tissue = "Tissue", snr_thresh = 3, sample_type = "sample_ID")

  expect_true(all(is.na(spectra(new_data, "snr")[, 5:6])))

  expect_equal(spectra(new_data, "snr")[1, ], c(NA, NA, 5, 5, NA, NA))
  expect_equal(spectra(new_data, "snr")[2, ], c(25, 25, 25, 25, NA, NA))
  expect_equal(spectra(new_data, "snr")[3, ], c(25.0, 25.0, 17.5, 17.5, NA, NA))
  expect_equal(spectra(new_data, "snr")[4, ], c(12.5, 12.5, 37.5, 37.5, NA,  NA))

})
