require(testthat)
require(quantMSImageR)

context("test mrm data read in")

test_that("read_mrm function", {

  # load file
  mrm_folder = system.file('extdata', package = 'quantMSImageR')
  mrm_file = "tissue_MRM_data"

  lib_ion_path = sprintf("%s/ion_library.txt",system.file('extdata', package = 'quantMSImageR'))

  test_data = read_mrm(name = mrm_file, folder = mrm_folder, lib_ion_path = lib_ion_path, polarity="Positive")

  expect_equal(dim(test_data)[[1]], 2)
  expect_equal(dim(test_data)[[2]], 11865)

  expect_equal(as.numeric(spectra(test_data)[1, 856]), 12273448)
  expect_equal(as.numeric(spectra(test_data)[2, 3571]), 13094678.4)

  expect_equal(fData(test_data)$analyte[1], "Analyte")
  expect_equal(fData(test_data)$analyte[2], "IS")
  expect_equal(fData(test_data)$name[1], "SM(d18:1/16:0)")
})

