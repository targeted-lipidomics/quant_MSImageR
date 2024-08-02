require(testthat)
require(quantMSImageR)

context("test int/response values changed to concentration")

test_that("int2conc function", {

  # Create test data
  test_data = readRDS(file = sprintf("%s/tissue_MRM_data.raw/cal_curve_MSI.RDS", system.file('extdata', package = 'quantMSImageR')))

  new_data <- int2conc(MSIobject = test_data,
                       val_slot = "intensity",
                       pixels = c("Tissue", "Noise"))

  expect_equal(signif(max(spectra(new_data, "conc - pg/pixel")), 4), 0.1025)
  expect_equal(signif(max(spectra(new_data, "conc - pg/mm2")), 4), 41.01)

})
