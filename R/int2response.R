library(Cardinal)

setGeneric("int2response", function(MSIobject, ...) standardGeneric("int2response"))

#' Function to normalise the intensity values to response per pixel if internal standard is present.
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject MSI MSIobject from Cardinal
#' @param val_slot character defining slot name to normalise - takes "intensity" as default
#' @param IS_name character defining the IS to use in fData(MSIobject) under the analyte header,
#' @return MSIobject with intensity values replaced with response
#'
#' @export
setMethod("int2response", "quant_MSImagingExperiment",
          function(MSIobject, val_slot = "intensity", IS_name = "IS", ...){

            if(!any(fData(MSIobject)$analyte == IS_name)){
              print("No IS in this study so normalisation skipped")

              spectra(MSIobject, "response") = spectraData(MSIobject)[[val_slot]]

              return(MSIobject)
            }

            IS_ind = which(fData(MSIobject)$analyte == IS_name)

            spectra(MSIobject, "response") = matrix(nrow = nrow(MSIobject), ncol = ncol(MSIobject))

            # Iterate over samples in study
            for(sample in unique(pData(MSIobject)$sample_ID)){

              # Select IS mz and pixels for specific sample
              pixel_ind = which(pData(MSIobject)$sample_ID == sample)

              # Save IS intensity vector
              IS_vec = spectraData(MSIobject)[[val_slot]][IS_ind, pixel_ind]

              for(mz_ind in 1:nrow(fData(MSIobject))){

                #print(sprintf("mz - %s", mz_ind))

                ints = spectraData(MSIobject)[[val_slot]][mz_ind, pixel_ind]
                response = ints / IS_vec
                spectra(MSIobject, "response")[mz_ind, pixel_ind] = response

              }
            }

            # Remove IS m/z
            MSIobject = MSIobject[-IS_ind, ]

            return(MSIobject)
          })
