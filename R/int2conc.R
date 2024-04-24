library(Cardinal)

setGeneric("int2conc", function(MSIobject, ...) standardGeneric("int2conc"))

#' Function to update intensity with concentration values
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject MSI object from Cardinal
#' @param pixels Label in pixel metadata the pixels to quantify
#' @return MSIobject with intensity values replaced with concentration values (ng/pixel)
#'
#' @export
setMethod("int2conc", "quant_MSImagingExperiment",
          function(MSIobject, pixels = "Tissue"){

            cal_list = MSIobject@calibrationInfo@cal_list
            pixel_inds = which(pData(MSIobject)$sample_type == pixels)
            MSIobject = MSIobject[, pixel_inds]

            no_cal_indices = c()
            for(i in 1:nrow(fData(MSIobject))){

              image(MSIobject, mz = mz(MSIobject)[i], contrast.enhance="histogram",
                    superpose = FALSE)

              eqn = cal_list[[i]]

              if(typeof(eqn) != "list"){

                print(sprintf("No calibration curve fo m/z %s", i))
                no_cal_indices = c(no_cal_indices, i)

              } else{

                pixel_ints = spectra(MSIobject)[i,]
                pixel_concs = as.vector(sapply(pixel_ints, function(y) inverse.predict(eqn, y)$Prediction))

                spectra(MSIobject)[i,] = pixel_concs

              }
            }

            if(length(no_cal_indices) > 0){
              MSIobject = MSIobject[-no_cal_indices, ]
            }

            return(MSIobject)
          })
