library(Cardinal)
library(dplyr)

setGeneric("createMSIDatamatrix", function(MSIobject, ...) standardGeneric("createMSIDatamatrix"))

#' Function to create data matrix from MSI object
#' @import Cardinal
#' @import dplyr
#' @include setClasses.R
#'
#' @param MSIobject MSI object from Cardinal
#' @param inputNA Whether to convert 0's in matrix to NA (default = T)
#' @return MSIobject with slots updated for i) matrix of average ng/pixel of m/z (rows = m/z and cols = cal level) in tissue ROIs ii) sample/ROI metadata
#'
#' @export
setMethod("createMSIDatamatrix", "quant_MSImagingExperiment",
          function(MSIobject, inputNA = T){


            MSIobject = MSIobject[, - which(is.na(pData(MSIobject)$ROI))]

            pixel_df = data.frame(pData(MSIobject)) %>%
              subset(!is.na(ROI)) %>%
              tibble::rownames_to_column("pixel_ind")

            df = data.frame(matrix(ncol = length(unique(pData(MSIobject)$sample_ID)),
                                   nrow = nrow(fData(MSIobject))))

            colnames(df) = unique(pData(MSIobject)$sample_ID)
            rownames(df) = fData(MSIobject)@mz

            for(roi in unique(pixel_df$sample_ID)){

              pixel_inds = as.numeric(pixel_df$pixel_ind[which(pixel_df$sample_ID == roi)])

              temp_df = data.frame(spectra(MSIobject)[, pixel_inds])
              rowMeans(temp_df, na.rm=T)

              df[[roi]] = rowMeans(temp_df, na.rm=T)

            }

            if(inputNA){
              df <- replace(df, df==0, NA)
            }

            MSIobject@tissueInfo@conc_matrix = df
            MSIobject@tissueInfo@sample_metadata = pixel_df

            return(MSIobject)
          })
