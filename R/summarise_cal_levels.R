library(Cardinal)
library(dplyr)

setGeneric("summarise_cal_levels", function(MSIobject, ...) standardGeneric("summarise_cal_levels"))

#' Function to calculate the mean response or intensity per pixel for the ROI at each calibration level across all calibration replicates (ng/pixel).
#' @import Cardinal
#' @import dplyr
#'
#' @param MSIobject MSI object from Cardinal
#' @param cal_label Label in pixel metadata which corresponds to calibration data
#' @return MSIobject with slots updated for i) matrix of average ng/pixel of m/z (rows = m/z and cols = cal level) ii) list of pixel counts per cal level
#'
#' @export
setMethod("summarise_cal_levels", "quant_MSImagingExperiment",
          function(MSIobject, cal_label = "Cal"){

            # create pixel data to associate pixel indices to cal levels
            pixel_data = data.frame(pData(MSIobject)) %>%
              tibble::rownames_to_column("pixel_ind") %>%
              subset(sample_type == cal_label) %>%
              subset(!is.na(ROI))

            # Create empty df to add to
            df = data.frame(matrix(ncol = length(unique(pixel_data$sample_ID)),
                                   nrow = nrow(fData(MSIobject))))
            colnames(df) = unique(pixel_data$sample_ID)
            rownames(df) = fData(MSIobject)@mz

            pixel_count = list()

            for(sample in unique(pixel_data$sample_ID)){
              inds = which(pixel_data$sample_ID == sample)

              pixels = as.numeric(pixel_data$pixel_ind[inds])

              pixel_count[[sample]] = length(pixels)

              mean_int_vec = c()

              for(mz_ind in 1:nrow(fData(MSIobject))){

                ints = spectra(MSIobject)[mz_ind, pixels]
                ints = replace(ints, ints ==0, NA)

                mean_int_per_pixel = mean(ints, na.rm=T) / length(pixels)

                mean_int_vec = c(mean_int_vec, mean_int_per_pixel)
              }

              df[[sample]] = mean_int_vec

            }

            MSIobject@calibrationInfo@response_per_pixel = df
            MSIobject@calibrationInfo@pixels_per_level = pixel_count


            return(MSIobject)
          })
