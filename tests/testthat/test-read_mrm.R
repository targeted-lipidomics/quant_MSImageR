require(testthat)
require(quantMSImageR)

context("test mrm data read in")

test_that("read_mrm function", {

  # load file
  mrm_folder = system.file('extdata', package = 'quantMSImageR')
  mrm_file = "tissue_MRM_data"

  lib_ion_path = sprintf("%s/ion_library.txt",system.file('extdata', package = 'quantMSImageR'))

  test_data = read_mrm(name = mrm_file, folder = mrm_folder, lib_ion_path = lib_ion_path)

  expect_equal(dim(test_data)[[1]], 4)
  expect_equal(dim(test_data)[[2]], 63070)

  expect_equal(as.numeric(spectra(test_data)[1, 856]), 1129)
  expect_equal(as.numeric(spectra(test_data)[2, 3571]), 852)

  expect_equal(unique(fData(test_data)$analyte), "Analyte")
  expect_equal(fData(test_data)$name[4], "12_13-DiHOME")
})

