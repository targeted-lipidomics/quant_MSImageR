require(reticulate)
require(XML)


#' Function to create data matrix from MSI object
#' @import XML
#' @import reticulate
#'
#' @param czi_fn Name of czi file
#' @param plate_x width (in um) of slide (75000 or 76000)
#' @param type Type of staining
#'
#' @return MSIobject with slots updated for i) matrix of average ng/pixel of m/z (rows = m/z and cols = cal level) in tissue ROIs ii) sample/ROI metadata
#'
#' @export extract_HnE_coords
extract_HnE_coords = function(czi_fn, plate_x = 75000, type = "H&E"){
  #reticulate::py_install("czifile")

  czi_pack = reticulate::import("czifile") # import czifile

  czi = czi_pack$CziFile(czi_fn)

  czi_xml_str = czi$metadata()

  xml_data <- XML::xmlToList(czi_xml_str)

  # Image centre position
  cent_pos = xml_data$Metadata$Experiment$ExperimentBlocks$AcquisitionBlock$SubDimensionSetups$MultiTrackSetup$SubDimensionSetups$RegionsSetup$SampleHolder$TileRegions$TileRegion$CenterPosition
  cent_pos = as.numeric( strsplit(cent_pos, ",")[[1]] )

  # Image size
  im_size = xml_data$Metadata$Experiment$ExperimentBlocks$AcquisitionBlock$SubDimensionSetups$MultiTrackSetup$SubDimensionSetups$RegionsSetup$SampleHolder$TileRegions$TileRegion$ContourSize
  im_size = as.numeric( strsplit(im_size, ",")[[1]] )


  plate_x - cent_pos[1] - (im_size[1] / 2)

  tibble::tibble(
    image = type,
    x_pixels = as.numeric(xml_data[["Metadata"]][["Information"]][["Image"]][["SizeX"]]),
    y_pixels = as.numeric(xml_data[["Metadata"]][["Information"]][["Image"]][["SizeY"]]),
    x_um = im_size[1],
    y_um = im_size[2],
    centre_x = cent_pos[1],
    centre_y = cent_pos[2]) %>%
    dplyr::mutate(x_step = x_um / x_pixels,
                  y_step = y_um / y_pixels,
                  origin_x = plate_x + cent_pos[1] - (im_size[1] / 2),
                  origin_y = cent_pos[2] - (im_size[2] / 2))
}
