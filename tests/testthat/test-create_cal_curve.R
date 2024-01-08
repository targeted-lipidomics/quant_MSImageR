require(testthat)
require(quantMSImageR)

context("test calibration curves are created")

test_that("create_cal_curve function", {


  # Create test data
  fdata <- MassDataFrame(mz=c(500, 510, 540, 550))
  pdata <- PositionDataFrame(run="run1",
                             coord=expand.grid(x=1:4, y=1:5),
                             sample_type = "Cal",
                             replicate = "01",
                             ROI = c(rep(NA, 2), rep("00", 2),
                                     rep(NA, 2), rep("01", 2),
                                     rep(NA, 2), rep("02", 2),
                                     rep(NA, 2), rep("03", 2),
                                     rep(NA, 2), rep("04", 2)))
  pdata$sample_ID = sprintf("%s_rep%s_%s", pdata$sample_type, pdata$replicate, pdata$ROI)

  ints <- matrix(nrow=4, ncol=20, byrow = T,
                 data = rep(c(rep(1, 4), rep(2, 4), rep(10, 4),
                              rep(100, 4), rep(1000, 4)), 4))

  test_data <- MSImagingExperiment(imageData=ints,
                                   featureData=fdata,
                                   pixelData=pdata)

  test_data <- as(test_data, "quant_MSImagingExperiment")
  test_data@calibrationInfo@cal_metadata = data.frame(sample = unique(pdata$sample_ID)[-grep("NA", unique(pdata$sample_ID))],
                                                      amount_ng = c(0, 5, 50, 450, 4800),
                                                      level = c("L0", "L1", "L2", "L3", "L4"))

  # Generate cal ROI summaries for generate curve
  test_data <- summarise_cal_levels(test_data)

  # Test standard addition
  new_data_std_add <- create_cal_curve(test_data, cal_type = "std_addition")

  coeffs = as.numeric(new_data_std_add@calibrationInfo@cal_list[[1]]$coefficients)

  expect_equal(round(coeffs[1], 4), 0.1595)
  expect_equal(round(coeffs[2], 4), 0.208)
  expect_equal(round(new_data_std_add@calibrationInfo@r2_df$r2[1], 5), 0.99996)

  # Test calibration
  new_data_cal <- create_cal_curve(test_data, cal_type = "cal")

  coeffs = as.numeric(new_data_cal@calibrationInfo@cal_list[[1]]$coefficients)

  expect_equal(round(coeffs[1], 4), 0.9439)
  expect_equal(round(coeffs[2], 4), 0.208)
  expect_equal(round(new_data_cal@calibrationInfo@r2_df$r2[1], 5), 0.99996)

})
