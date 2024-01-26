library(Cardinal)
library(dplyr)

setGeneric("createMSIDatamatrix", function(MSIobject, ...) standardGeneric("createMSIDatamatrix"))

#' Function to create data matrix from MSI object
#' @import Cardinal
#' @import dplyr
#' @include setClasses.R
#'
#' @param MSIobject MSI object from Cardinal, pData to include..... sample_ID..
#' @param inputNA Whether to convert 0's in matrix to NA (default = T)
#' @param roi_header Header in pData pertaining to ROIs to average. Set to NA to skip generating average df
#' @return MSIobject with slots updated for i) matrix of average ng/pixel of m/z (rows = m/z and cols = cal level) in tissue ROIs ii) sample/ROI metadata
#'
#' @export
setMethod("createMSIDatamatrix", "quant_MSImagingExperiment",
          function(MSIobject, inputNA = T, roi_header = NA){

            # Subset pixels in ROIs only (based on roi_header)
            if(is.na(roi_header)){
              pData(MSIobject)$ROI = 1:nrow(pData(MSIobject))
            } else{
              pData(MSIobject)$ROI = pData(MSIobject)[[roi_header]]
              MSIobject = MSIobject[, - which(is.na(pData(MSIobject)$ROI))]
            }

            # Update pixel data
            pixel_df = data.frame(pData(MSIobject)) %>%
              subset(!is.na(ROI)) %>%
              tibble::rownames_to_column("pixel_ind")


            # All pixel df
            all_pixel_df = data.frame(spectra(MSIobject)[,])
            colnames(all_pixel_df ) = sprintf("pixel_%s", pixel_df$pixel_ind)
            rownames(all_pixel_df ) = fData(MSIobject)$name

            if(inputNA){
              all_pixel_df <- replace(all_pixel_df, all_pixel_df==0, NA)
            }


            # Create empty average df
            if(!is.na(roi_header)){
              ave_df = data.frame(matrix(ncol = length(unique(pData(MSIobject)$sample_ID)),
                                         nrow = nrow(fData(MSIobject))))

              colnames(ave_df) = unique(pData(MSIobject)$sample_ID)
              rownames(ave_df) = fData(MSIobject)$name


              for(roi in unique(pixel_df$sample_ID)){

                pixel_inds = as.numeric(pixel_df$pixel_ind[which(pixel_df$sample_ID == roi)])

                temp_df = data.frame(spectra(MSIobject)[, pixel_inds])
                rowMeans(temp_df, na.rm=T)

                ave_df[[roi]] = rowMeans(temp_df, na.rm=T)

              }

              if(inputNA){
                ave_df <- replace(ave_df, ave_df==0, NA)
              }

              MSIobject@tissueInfo@roi_average_matrix = ave_df
            }

            MSIobject@tissueInfo@all_pixel_matrix = all_pixel_df
            MSIobject@tissueInfo@sample_metadata = pixel_df

            return(MSIobject)
          })
