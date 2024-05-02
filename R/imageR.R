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
#' @param threshold percintile to remove from colour scale (i.e background noise)
#' @param sample_lab character header from pData(MSIobject) to label images (defaults to "sample_ID")
#' @param pixels character from pData(MSIobject)$sample_type to take pixels to plot (defaults to NA)
#' @param overlay whether to overlay features (not yet implemented!)
#' @param feat_ind Index of feature from fData(MSIobject) to image
#' @param perc_scale Logical to normalise scale to % (TRUE) or displar raw values (FALSE)
#' @return ggplot object
#'
#' @export
setMethod("imageR", "quant_MSImagingExperiment",
          function(MSIobject, value = "response %", scale = "suppress", threshold = 1,
                   sample_lab = "sample_ID", pixels = NA, percentile=99.0, overlay = F,
                   feat_ind = 1, perc_scale = F){

            MSIobject = as(MSIobject[feat_ind, ], "quant_MSImagingExperiment")

            if(!is.na(pixels)){
              MSIobject = MSIobject[, which(pData(MSIobject)$sample_type == pixels)]
            }

            # Generate image data matrix
            image_df = tibble(x = pData(MSIobject)@coord@listData[["x"]],
                              y = pData(MSIobject)@coord@listData[["y"]],
                              response = as.numeric(spectra(MSIobject)),
                              sample = pData(MSIobject)[[sample_lab]],
                              feature = fData(MSIobject)$name)


            if(scale == "sqrt"){
              image_df = mutate(image_df, response = sqrt(response))
            }
            if(scale == "suppress"){
              vals = image_df$response
              max_val = quantile(vals, percentile /100, na.rm=TRUE)
              if (max_val > min(vals, na.rm=TRUE)){
                vals[vals > max_val] = max_val
              }

              image_df$response = vals
            }
            if(scale == "histogram"){

              # Cardinal implementation	https://rdrr.io/bioc/Cardinal/src/R/DIP.R

              vals = image_df$response

              unique_percentiles = unique(quantile(vals, seq(from=0, to=1, length.out=100), na.rm=TRUE))
              if (length(unique_percentiles) < 5){
                return("Range of intensity values too narrow for histogram normalization.") }

              vals_interval = cut(vals, unique_percentiles, include.lowest=TRUE)
              vals_new = as.numeric(vals_interval) / length(levels(vals_interval))

              scale = mean(vals, na.rm=TRUE) / mean(vals_new, na.rm=TRUE) # So mean value same as prior to scaling
              vals_new = scale * vals_new

              image_df$response = vals_new

            }

            if(!is.na(threshold)){
              print("threshold")

              # Remove x% of values - threshold
              vals = image_df$response
              min_val = quantile(vals, threshold /100, na.rm=TRUE)
              if (min_val < max(vals, na.rm=TRUE)){
                vals[vals < min_val] = 0
              }
              image_df$response = vals
            }

            if(overlay ==T){
              return("Overlay functionality not added yet. Set to FALSE and rerun.")
            }

            # Set NA and negative values to 0 for plotting!
            image_df = image_df %>%
              mutate(response = ifelse(is.na(response), 0, response)) %>%
              mutate(response = ifelse(response < 0, 0, response))

            if(perc_scale == T){
              vals = image_df$response
              if(max(vals) > 0){
                perc_vals = 100 * (vals / max(vals))

                image_df$response = perc_vals

              }
            }

            p = ggplot(data=image_df,aes(x=x,y=-y,fill=response))+
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

            p

            return(p)

        })
