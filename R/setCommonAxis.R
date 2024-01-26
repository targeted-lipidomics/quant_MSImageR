library(Cardinal)

setGeneric("setCommonAxis", function(MSIobjects, ...) standardGeneric("setCommonAxis"))

#' Function to update intensity with concentration values
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobjects List of MSIobjects
#' @param ref_object MSIobject with all transitions of interets included
#' @return List of MSIobjects with common fData and matching spectra
#'
#' @export
setMethod("setCommonAxis", "list",
          function(MSIobjects, ref_object){

            features = lapply( 1:length(MSIobjects), FUN=function(x){
              data.frame(fData(MSIobjects[[x]])) })

            ref_features = data.frame(fData(ref_object))

            # Common axis
            for(i in 1:length(MSIobjects)){

              if(all(dim(ref_features) == dim(features[[i]]))){
                if(all(ref_features == features[[i]])) next
              }

              feat = features[[i]]

              # Correct fData
              feat = merge(x=ref_features, y=feat, by = c("precursor_mz", "product_mz"),
                           all.x = T, all.y = F, suffixes = c("","_old")) %>%
                mutate(name = name_old, analyte = analyte_old) %>%
                arrange(mz)

              # Set empty iData
              idata = matrix(nrow=nrow(feat), ncol=ncol(MSIobjects[[i]]))

              # Update iData to spectral info corresponding to fData channels
              for(mz in (feat %>% subset(!is.na(analyte)))$mz){

                old_mz = feat$mz_old[which(feat$mz == mz)]

                idata[mz, ] = spectra(MSIobjects[[i]][old_mz,])

              }

              MSIobjects[[i]] <- MSImagingExperiment(imageData= idata,
                                         featureData= MassDataFrame(mz=as.numeric(feat$mz),
                                                                   analyte = feat$analyte,
                                                                   precursor_mz = as.numeric(feat$precursor_mz),
                                                                   product_mz = as.numeric(feat$product_mz),
                                                                   name = feat$name),
                                         pixelData=pData(MSIobjects[[i]]))
            }

            return(MSIobjects)

          })
