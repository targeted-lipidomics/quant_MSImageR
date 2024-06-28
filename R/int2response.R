library(Cardinal)

setGeneric("int2response", function(MSIobject) standardGeneric("int2response"))

#' Function to normalise the intensity values to response per pixel if internal standard is present.
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject MSI MSIobject from Cardinal
#' @return MSIobject with intensity values replaced with response
#'
#' @export
setMethod("int2response", "quant_MSImagingExperiment",
          function(MSIobject){

            if(!any(fData(MSIobject)$analyte == "IS")){
              print("No IS in this study so normalisation skipped")
              return(MSIobject)
            }

            IS_ind = which(fData(MSIobject)$analyte == "IS")

            # Iterate over samples in study
            for(sample in unique(pData(MSIobject)$sample_ID)){

              #print(sample)

              # Select IS mz and pixels for specific sample
              pixel_ind = which(pData(MSIobject)$sample_ID == sample)

              # Save IS intensity vector
              IS_vec = spectra(MSIobject)[IS_ind, pixel_ind]

              for(mz_ind in 1:nrow(fData(MSIobject))){

                #print(sprintf("mz - %s", mz_ind))

                ints = spectra(MSIobject)[mz_ind, pixel_ind]
                response = ints / IS_vec
                spectra(MSIobject)[mz_ind, pixel_ind] = response

              }
            }

            # Remove IS m/z
            MSIobject = MSIobject[-IS_ind, ]

            return(MSIobject)
          })
