library(Cardinal)
library(dplyr)
library(chemCal)

setGeneric("create_cal_curve", function(MSIobject, ...) standardGeneric("create_cal_curve"))

#' Function to create calibration curves (response v concentration, where concentration is ng/pixel)
#'
#' @import Cardinal
#' @import dplyr
#' @import chemCal
#' @include setClasses.R
#'
#' @param response_matrix matrix of average ng/pixel of m/z (rows = m/z and cols = cal level)
#' @param cal_type string of approach to generate claibration curve - 'std_addition' is default
#' @return MSIobject with slots updated for i) cal_list - List of linear models for each m/z (response v concentration, where concentration is ng/pixel) and ii) r2 values for each calibration iii) calibration metadata
#'
#' @export
setMethod("create_cal_curve", "quant_MSImagingExperiment",
          function(MSIobject, cal_type = "std_addition"){

            pixel_count = MSIobject@calibrationInfo@pixels_per_level

            cal_metadata = MSIobject@calibrationInfo@cal_metadata %>%
              mutate(pixel_count = sapply(sample,
                                          function(sample_name) pixel_count[[sample_name]]),
                     ng_per_pixel = amount_ng / pixel_count)

            response_matrix = MSIobject@calibrationInfo@response_per_pixel

            cal_list = list()

            # Data frame of r2 values for each calibration
            r2_df = data.frame(mz=rep(NA,nrow(response_matrix)),
                               r2 = rep(NA,nrow(response_matrix)))

            for(i in 1:nrow(response_matrix)){

              mz = as.character(rownames(response_matrix)[i])
              #print(mz)

              ints = response_matrix[i,]

              if(sum(!is.na(ints)) > 2){

                temp_df = data.frame(int = as.numeric(ints),
                                     sample_name = colnames(response_matrix))
                temp_df$ng_per_pixel = sapply(temp_df$sample_name,
                                              function(sample_name) cal_metadata$ng_per_pixel[which(cal_metadata$sample == sample_name)])

                eqn = lm(int~ng_per_pixel, data = temp_df, na.action = na.exclude)

                if(cal_type == "std_addition"){

                  # Calculate background conc
                  background_conc = inverse.predict(eqn, 0)$Prediction

                  # Update conc values (original conc - background)
                  temp_df$ng_per_pixel = temp_df$ng_per_pixel - background_conc

                  # Update equation
                  eqn = lm(int~ng_per_pixel, data = temp_df, na.action = na.exclude)

                  r2_df$mz[i] = mz
                  r2_df$r2[i] = summary(eqn)[["r.squared"]]
                }

                cal_list[[mz]] = eqn
              } else{
                cal_list[[mz]] = "NO DATA"
              }
            }

            MSIobject@calibrationInfo@cal_list = cal_list
            MSIobject@calibrationInfo@r2_df = r2_df
            MSIobject@calibrationInfo@cal_metadata = cal_metadata

            return(MSIobject)
          })
