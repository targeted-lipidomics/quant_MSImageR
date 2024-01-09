require(testthat)
require(quantMSImageR)

context("test mrm data read in")

test_that("read_mrm function", {

  # load file
  mrm_folder = file.path(system.file(package="quantMSImageR"), "data")
  #mrm_folder = file.path("C:/Users/matsmi/OneDrive - Karolinska Institutet/Dokument/MSI/quantMSImageR", "data")
  mrm_file = "example_MRM_raw"

  test_data = read_mrm(name = mrm_file, folder = mrm_folder)

  expect_equal(dim(test_data)[[1]], 2)
  expect_equal(dim(test_data)[[2]], 49686)

  expect_equal(as.numeric(spectra(test_data)[1, 856]), 2600166)
  expect_equal(as.numeric(spectra(test_data)[2, 3571]), 10322)

})
