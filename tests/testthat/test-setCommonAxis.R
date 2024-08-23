require(testthat)
require(quantMSImageR)

context("test that list of MSImagingExperiment objects set to a specific feature axis")

test_that("setCommonAxis function", {

  # Create test data
  ref_fdata = MassDataFrame(mz=1:4,
                            analyte = c("IS", rep("analyte", 3)),
                            precursor_mz = c(100, 200, 300, 400),
                            product_mz = c(100, 200, 300, 400)/2,
                            name = c("IS", sprintf("lipid_%s", 1:3)))

  fdata1 <- MassDataFrame(mz=1:5,
                         analyte = c("analyte", "IS", rep("overwrite1", 3)),
                         precursor_mz = c(400, 100, 3000, 4000, 5000),
                         product_mz = c(400, 100, 3000, 4000, 5000)/2,
                         name = c("lipid_3", "IS", sprintf("overwrite_%s", 1:3)))

  fdata2 <- MassDataFrame(mz=6:8,
                          analyte = c(rep("overwrite2", 2), "IS"),
                          precursor_mz = c(6000, 7000, 100),
                          product_mz = c(6000, 7000, 100)/2,
                          name = c(sprintf("overwrite_%s", 4:5), "IS"))

  pdata1 <- PositionDataFrame(run=c(rep("run1", 2), rep("run2", 2)),
                              coord=expand.grid(x=1:2, y=1:2),
                              sample_ID = c(rep("run1", 2), rep("run2", 2)))

  pdata2 <- PositionDataFrame(run=c(rep("run3", 2), rep("run4", 2)),
                              coord=expand.grid(x=1:2, y=1:2),
                              sample_ID = c(rep("run3", 2), rep("run4", 2)))

  ints1 <- matrix(nrow=5, ncol=4,
                  data = c(rep(c(5,10,15, 100), 2),
                           rep(c(10,10,15, 100), 3)))

  ints2 <- matrix(nrow=3, ncol=4,
                  data = c(rep(c(3,6,9, 30), 1),
                           rep(c(6,6,9, 30), 2)))

  MSIobject1 <- MSImagingExperiment(spectraData=ints1,
                                    featureData=fdata1,
                                    pixelData=pdata1)

  MSIobject2 <- MSImagingExperiment(spectraData=ints2,
                                    featureData=fdata2,
                                    pixelData=pdata2)

  test_data = setCommonAxis(MSIobjects = list(MSIobject1, MSIobject2), ref_fdata = ref_fdata)

  expect_equal(length(test_data), 2)



  expect_equal(nrow(fData(test_data[[1]])), 4)

  expect_true(all(is.na(spectra(test_data[[1]])[2, ])))
  expect_true(all(is.na(spectra(test_data[[1]])[3, ])))

  expect_true(all(is.na(spectra(test_data[[2]])[2, ])))
  expect_true(all(is.na(spectra(test_data[[2]])[3, ])))
  expect_true(all(is.na(spectra(test_data[[2]])[4, ])))

  expect_equal(spectra(test_data[[1]])[1, ], c(10, 15, 100, 10))
  expect_equal(spectra(test_data[[1]])[4, ], c(5, 10, 15, 100))

  expect_equal(spectra(test_data[[2]])[1, ], c(9, 6, 6, 30))

})

