library(Cardinal)
library(dplyr)
library(chemCal)
library(viridis)

setGeneric("imageR", function(MSIobject, ...) standardGeneric("imageR"))

#' Function to create ion images (using ggplot)
#'
#' @import Cardinal
#' @import dplyr
#' @import chemCal
#' @include setClasses.R
#'
#' @param value character stating what values pertain to
#' @param scale suppress, histogram, sqrt
#' @param percentile percentile to suppress colour scale (for suppress.)if scale == "suppress")
#' @param threshold
#' @param sample_lab character header from pData(MSIobject) to label images (defaults to "sample_ID")
#' @param pixels character from pData(MSIobject)$sample_type to take pixels to plot
#' @return ggplot object
#'
#' @export
setMethod("imageR", "quant_MSImagingExperiment",
          function(MSIobject, value = "response %", scale = "suppress", threshold = 10,
                   sample_lab = "sample_ID", pixels = "Tissue", percentile=99.0){   #MSIobject = msi_combined_snr

            MSIobject = MSIobject[, which(pData(MSIobject)$sample_type == pixels)]

            image_df = tibble(x=numeric(), y=numeric(), response = numeric(), sample=character())

            # Iterate over features in study
            for(mz_ind in 1:nrow(fData(MSIobject))){

              # Generate image data matrix
              temp_df = tibble(x = pData(MSIobject)@coord@listData[["x"]],
                                y = pData(MSIobject)@coord@listData[["y"]],
                                response = spectra(MSIobject)[mz_ind, ],
                                sample = pData(MSIobject)[[sample_lab]],
                                feature = fData(MSIobject)$name[mz_ind])


              if(scale == "sqrt"){
                temp_df = mutate(temp_df, response = sqrt(response))
              }
              if(scale == "suppress"){
                vals = temp_df$response
                cutoff <- quantile(vals, percentile /100, na.rm=TRUE)
                if (cutoff > min(vals, na.rm=TRUE)){
                  vals[vals > cutoff] <- cutoff
                }

                temp_df = mutate(temp_df, response = vals)
              }
              if(scale == "histogram"){
                print("not implemented")
              }


              # Remove x% of values - threshold
              vals = temp_df$response
              cutoff <- quantile(vals, threshold /100, na.rm=TRUE)
              if (cutoff < max(vals, na.rm=TRUE)){
                vals[vals < cutoff] <- 0
              }
              temp_df = mutate(temp_df, response = vals)


              # Set NA to 0 for plotting!
              temp_df = mutate(temp_df, response = ifelse(is.na(response), 0, response))

              image_df = rbind(image_df, temp_df)

            }

            p = ggplot(data=temp_df,aes(x=x,y=-y,fill=response))+
              geom_tile() +
              theme_minimal() +
              theme(axis.title = element_blank(),
              axis.text = element_blank(),
                axis.line = element_blank(),
                panel.grid = element_blank(),
                plot.title = element_text(hjust = 0.5, face="bold", size = 15)) +
              scale_fill_viridis(na.value = "white") +
              labs(fill=value) +
              facet_grid(sample~feature)

            return(p)

          })
