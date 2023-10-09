#' Function to convert intensity values from 0 to NA in MSI dataset.
#'
#' @param MSIobject MSI object from Cardinal
#' @return MSIobject with intensity values replaced with response
#'
#' @export
setMethod("zero2na", "quant_MSImagingExperiment",
          function(MSIobject){

            for(i in 1:length(mz(MSIobject))){

              #print(i)

              ints = spectra(MSIobject)[i,]
              if(all(is.na(ints))){
                print(sprintf("all intensities are NA for m/z %s. Doing nothing.", i))
              }
              else if(sum(ints, na.rm=T) == 0){
                print(sprintf("all intensities are 0 for m/z %s. Making NA.", i))
                spectra(MSIobject)[i, ] = NA
              }
            }

            return(MSIobject)

          })
