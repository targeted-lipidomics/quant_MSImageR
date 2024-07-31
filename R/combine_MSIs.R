library(Cardinal)

setGeneric("combine_MSIs", function(MSIobject, ...) standardGeneric("combine_MSIs"))

#' Function to convert intensity values from 0 to NA in MSI dataset.
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject MSI object from Cardinal
#' @return MSIobject with intensity values replaced with response
#' @param ... additional MSI objects to combine
#'
#' @export
setMethod("combine_MSIs", "MSImagingExperiment",
          function(MSIobject, ...){

            objects <- c(as.list(environment()), list(...))

            f_data = fData(MSIobject)

            for(ind in 2:length(objects)){

              MSIobject = as( cbind(MSIobject, objects[[ind]]), 'MSImagingExperiment')

            }

            MSIobject = as(MSIobject, "quant_MSImagingExperiment")


            fData(MSIobject) = f_data

            return(MSIobject)

          })
