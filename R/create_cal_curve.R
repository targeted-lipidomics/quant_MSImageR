library(Cardinal)
library(dplyr)
library(chemCal)

setGeneric("create_cal_curve", function(MSIobject, ...) standardGeneric("create_cal_curve"))

#' Function to create calibration curves (response v concentration, where concentration is pg/pixel)
#'
#' @import Cardinal
#' @import dplyr
#' @import chemCal
#' @include setClasses.R
#'
#' @param response_matrix matrix of average pg/pixel of m/z (rows = m/z and cols = cal level)
#' @param cal_type string of approach to generate claibration curve - 'std_addition' if standards are on tissue and 'cal' if direct onto glass slide.
#' @return MSIobject with slots updated for i) cal_list - List of linear models for each m/z (response v concentration, where concentration is pg/pixel) and ii) r2 values for each calibration iii) calibration metadata
#'
#' @export
setMethod("create_cal_curve", "quant_MSImagingExperiment",
          function(MSIobject, cal_type = "std_addition"){ #MSIobject = msi_combined_sumCal

            pixel_count = MSIobject@calibrationInfo@pixels_per_level

            cal_metadata = MSIobject@calibrationInfo@cal_metadata %>%
              mutate(pixel_count = sapply(sample,
                                          function(sample_name) pixel_count[[sample_name]])) %>%
              mutate(pg_per_pixel = amount_pg / pixel_count)

            background_sample = cal_metadata$sample[which(cal_metadata$pg_per_pixel == 0)]

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

                temp_df = data.frame(response = as.numeric(ints),
                                     sample_name = colnames(response_matrix)) %>%
                  mutate(lev = sapply(sample_name,
                                      function(sn) strsplit(sn, "_")[[1]][length(strsplit(sn, "_")[[1]])] ))

                temp_df$pg_per_pixel = sapply(temp_df$lev,
                                              function(l) cal_metadata$pg_per_pixel[which(cal_metadata$level == l)])

                eqn = lm(response~pg_per_pixel, data = temp_df, na.action = na.exclude)

                p = ggplot(temp_df, aes(x=pg_per_pixel, y = response)) +
                  stat_poly_line(method = "lm", col = "red", se=F, size = 2, linetype = "dashed") +
                  stat_poly_eq(use_label(c("eq", "R2"))) +
                  geom_point(size = 2) +
                  theme_Publication()

                if(cal_type == "std_addition"){

                  # Calculate background conc
                  background_conc = inverse.predict(eqn, 0)$Prediction

                  # Update conc values (original conc - background)
                  temp_df$pg_per_pixel = temp_df$pg_per_pixel - background_conc


                  temp_df = subset(temp_df, lev != "background")

                  # Update equation
                  eqn = lm(response~pg_per_pixel, data = temp_df, na.action = na.exclude)
                }

                r2_df$mz[i] = mz
                r2_df$r2[i] = summary(eqn)[["r.squared"]]
                cal_list[[mz]] = eqn

              } else{
                cal_list[[mz]] = "NO DATA"
              }
            }

            cal_metadata = dplyr::left_join(temp_df, cal_metadata , join_by(lev == level)) %>%
              rename(pg_per_pixel.std_add = pg_per_pixel.x,
                     pg_per_pixel = pg_per_pixel.y) %>%
              select(-c("sample"))

            MSIobject@calibrationInfo@cal_list = cal_list
            MSIobject@calibrationInfo@r2_df = r2_df
            MSIobject@calibrationInfo@cal_metadata = cal_metadata

            return(MSIobject)
          })
