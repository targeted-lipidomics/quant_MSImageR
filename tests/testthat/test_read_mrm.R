require(testthat)
require(quantMSImageR)

context("test mrm data read in")

test_that("read_mrm function", {

  # load file
  mrm_folder = system.file('extdata', package = 'quantMSImageR')
  mrm_file = "example_MRM_raw"

  lib_ion_path = sprintf("%s/ion_library.txt",system.file('extdata', package = 'quantMSImageR'))

  test_data = read_mrm(name = mrm_file, folder = mrm_folder, lib_ion_path = lib_ion_path, polarity="Positive")

  expect_equal(dim(test_data)[[1]], 2)
  expect_equal(dim(test_data)[[2]], 49686)

  expect_equal(as.numeric(spectra(test_data)[1, 856]), 2600166)
  expect_equal(as.numeric(spectra(test_data)[2, 3571]), 10322)

  expect_equal(fData(test_data)$analyte[1], "Unknown")
  expect_equal(fData(test_data)$name[1], "1")
  expect_equal(fData(test_data)$name[2], "SM(d18:1/16:0)")
})

