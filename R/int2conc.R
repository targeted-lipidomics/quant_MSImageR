#' Function to update intensity with concentration values
#'
#' @param MSIobject MSI object from Cardinal
#' @param cal_label Label in pixel metadata which corresponds to calibration data
#' @return MSIobject with intensity values replaced with concentration values (ng/pixel)
#'
#' @export
setMethod("int2conc", "quant_MSImagingExperiment",
          function(MSIobject, cal_label = "Cal"){

            cal_list = MSIobject@calibrationInfo@cal_list
            pixel_inds = which(pData(MSIobject)$sample_type != cal_label)
            MSIobject = MSIobject[, pixel_inds]

            no_cal_indices = c()
            for(i in 1:length(mz(MSIobject))){

              #print(i)

              #if(round(as.numeric(mz(MSIobject)[i], 3)) != round(as.numeric(names(cal_list)[i]), 3)){
              #  return("mz in cal not match mz from MSIobject")
              #}

              image(MSIobject, mz = mz(MSIobject)[i], contrast.enhance="histogram",
                    superpose = FALSE)

              eqn = cal_list[[i]]

              if(typeof(eqn) != "list"){

                print(sprintf("No calibration curve fo m/z %s", i))
                no_cal_indices = c(no_cal_indices, i)

              } else{

                tissue_ints = spectra(MSIobject)[i,]
                tissue_concs = as.vector(sapply(tissue_ints, function(y) inverse.predict(eqn, y)$Prediction))

                spectra(MSIobject)[i,] = tissue_concs

              }
            }

            MSIobject = MSIobject[-no_cal_indices, ]
          })
