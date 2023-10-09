#### Class to store information about the calibration spots
## --------------------------------------------------
calibrationInfo = setClass("calibrationInfo",
         slots = c(
           cal_metadata = "data.frame",
           cal_list = "list",
           response_per_pixel = "data.frame",
           pixels_per_level = "list",
           r2_df = "data.frame"
         )
)

#### Class to store information about the calibration spots
## --------------------------------------------------
tissueInfo = setClass("tissueInfo",
                           slots = c(
                             conc_matrix = "data.frame",
                             sample_metadata = "data.frame",
                             feature_metadata = "data.frame"
                           )
)

#### Class containting calibration metadata and MSImaging experiment
## --------------------------------------------------
quant_MSImagingExperiment = setClass("quant_MSImagingExperiment",
         contains = 'MSContinuousImagingExperiment',
         slots = c(
           calibrationInfo = "calibrationInfo",
           tissueInfo = "tissueInfo"
         )
)
