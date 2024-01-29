require(testthat)
require(quantMSImageR)

context("test int/response values changed to concentration")

test_that("int2conc function", {

  # Create test data
  fdata <- MassDataFrame(mz=c(500, 510, 540, 550))

  pdata <- PositionDataFrame(run= "run1",
                             coord=expand.grid(x=1:2, y=1:2),
                             sample_type = c(rep("tissue", 3), rep("Cal", 1)))

  ints <- matrix(nrow=4, ncol=4,
                 data = c(rep(c(50,100,105, 1500), 2),
                          rep(c(60,105,150, 1400), 2)))

  test_data <- MSImagingExperiment(imageData=ints,
                                   featureData=fdata,
                                   pixelData=pdata)

  test_data <- as(test_data, "quant_MSImagingExperiment")
  eqn = lm(formula = int ~ ng_per_pixel,
           data = data.frame(int= c(8,10,150,1000),
                             ng_per_pixel = c(1, 20, 200 ,2400)))

  test_data@calibrationInfo@cal_list = list(eqn, eqn, eqn, eqn)

  new_data <- int2conc(test_data, cal_label = "Cal")

  spectra(test_data)
  spectra(new_data)

  expect_equal(as.numeric(nrow(new_data)), 4)
  expect_equal(as.numeric(ncol(new_data)), 3)

  expect_equal(round(spectra(new_data)[,1], 4), c(61.7621, 184.3836, 196.6457, 3617.7846))
  expect_equal(round(spectra(new_data)[,2], 4), round(spectra(new_data)[,1], 4))
  expect_equal(round(spectra(new_data)[,3], 4), c(86.2864, 196.6457, 307.0050, 3372.5417))

})
