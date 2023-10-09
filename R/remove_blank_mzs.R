#' Function to remove m/z values with no data.
#'
#' @param MSIobject MSI object from Cardinal
#' @return MSIobject with m/z values from experiment with no data removed.
#'
#' @export
setMethod("remove_blank_mzs", "quant_MSImagingExperiment",
          function(MSIobject){

            remove_inds = c()
            for(mz_ind in 1:nrow(fData(MSIobject))){
              if(all(is.na(spectra(MSIobject)[mz_ind, ]))){
                remove_inds = c(remove_inds, mz_ind)
              }
            }
            MSIobject = MSIobject[-remove_inds, ]

            return(MSIobject)
          })
