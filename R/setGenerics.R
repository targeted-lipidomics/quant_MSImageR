#### Set generic functions ####
## -----------------------------------------------
setGeneric("zero2na", function(MSIobject) standardGeneric("zero2na"))

setGeneric("remove_blank_mzs", function(MSIobject) standardGeneric("remove_blank_mzs"))

setGeneric("int2response", function(MSIobject) standardGeneric("int2response"))

setGeneric("summarise_cal_levels", function(MSIobject, ...) standardGeneric("summarise_cal_levels"))

setGeneric("create_cal_curve", function(MSIobject, ...) standardGeneric("create_cal_curve"))

setGeneric("int2conc", function(MSIobject, ...) standardGeneric("int2conc"))

setGeneric("createMSIDatamatrix", function(MSIobject, ...) standardGeneric("createMSIDatamatrix"))
