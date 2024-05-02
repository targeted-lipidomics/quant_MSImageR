library(Cardinal)

setGeneric("int2snr", function(MSIobject, ...) standardGeneric("int2snr"))

#' Function to convert the intensity values to SNR per pixel based on same transitions in noise/background pixels.Run after IS normalization.
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject MSI MSIobject from Cardinal - including pData(MSIobject)$sample_type == Noise
#' @param noise character in pData(MSIobject)$sample_type which indicates background / noise pixels
#' @param tissue character in pData(MSIobject)$sample_type which indicates tissue pixels to calculate SNR for
#' @param snr_thresh value indicating the minimum SNR value to accept (below this value SNR = 0)
#' @return MSIobject with intensity values replaced with SNR values
#'
#' @export
setMethod("int2snr", "quant_MSImagingExperiment",
          function(MSIobject, noise = "Noise", tissue = "Tissue", snr_thresh = 3, ..){ #MSIobject = msi_combined_response

            if(!any(pData(MSIobject)$sample_type == "Noise")){
              print("No noise pixels. Return same values")
              return(MSIobject)
            }

            #Set noise and tissue pixels
            noise_pixels = which(pData(MSIobject)$sample_type == noise)
            tissue_pixels = which(pData(MSIobject)$sample_type == tissue)

            # Iterate over features in study
            for(mz_ind in 1:nrow(fData(MSIobject))){

              #print(sprintf("mz - %s", mz_ind))

              # Save noise response vector
              noise_vec = spectra(MSIobject)[mz_ind, noise_pixels]

              #Deal with 0 values in noise vector - 10% of lowest
              noise_vec[which(noise_vec ==0)] = NA
              mv_impute = min(noise_vec, na.rm = T) / 10
              noise_vec[which(is.na(noise_vec))] = mv_impute

              noise_level = mean(noise_vec)

              # Calculate S/N of tissue pixels
              tissue_vec = spectra(MSIobject)[mz_ind, tissue_pixels]
              snr =  tissue_vec / noise_level
              snr =  ifelse(snr < snr_thresh, NA, snr)
              spectra(MSIobject)[mz_ind, tissue_pixels] = snr

            }

            return(MSIobject)
          })
