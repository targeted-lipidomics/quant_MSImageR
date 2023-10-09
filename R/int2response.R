#' Function to normalise the intensity values to response per pixel if internal standard is present.
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

            # Iterate over samples in study
            for(sample in unique(pData(MSIobject)$sample_ID)){

              #print(sample)

              # Select IS mz and pixels for specific sample
              pixel_ind = which(pData(MSIobject)$sample_ID == sample)
              IS_ind = which(fData(MSIobject)$analyte == "IS")

              # Save IS intensity vector
              IS_vec = spectra(MSIobject)[IS_ind, pixel_ind]

              #Deal with 0 values in IS vector - 10% of lowest
              IS_vec[which(IS_vec ==0)] = NA
              mv_impute = min(IS_vec, na.rm = T) / 10
              IS_vec[which(is.na(IS_vec))] = mv_impute

              for(mz_ind in 1:nrow(fData(MSIobject))){

                #print(sprintf("mz - %s", mz_ind))

                ints = spectra(MSIobject)[mz_ind, pixel_ind]
                response = ints / IS_vec
                spectra(MSIobject)[mz_ind, pixel_ind] = response

              }
            }
            return(MSIobject)
          })
