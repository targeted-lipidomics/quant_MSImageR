createDatamatrix <- function (MSIobject, inputNA = T) {
  MSIobject = as(MSIobject, "MSContinuousImagingExperiment")
  df <- data.frame(spectra(MSIobject))
  if (rownames(df)[1] == "1") {
    df <- df %>% tibble::add_column(mz = as.numeric(mz(MSIobject))) %>% 
      tibble::column_to_rownames(var = "mz") %>%
      dplyr::rename_all(~sprintf("pixel_%s", 1:ncol(MSIobject)))
  }
  if (inputNA) {
    df <- replace(df, df == 0, NA)
  }
  return(df)
}

zero2na = function(MSIobject){
  
  MSIobject = as(MSIobject, 'MSContinuousImagingExperiment')
  
  for(i in 1:length(mz(MSIobject))){
    
    print(i)
    
    ints = spectra(MSIobject)[i,]
    if(all(is.na(ints))){
      print("all intensities are NA. Doing nothing.")
    }
    else if(all(ints == 0)){
        print("all intensities are 0")
        spectra(MSIobject)[i, ] = NA
    }
  }
  
  return(MSIobject)
  
}

int2response = function(MSIobject){
  
  pData(MSIobject)$sample_ID = sprintf("%s_%s", pData(MSIobject)$sample_type, pData(MSIobject)$replicate)
  
  # Iterate over samples in study
  for(sample in unique(pData(MSIobject)$sample_ID)){
    
    print(sample)
    
    # Select IS mz and pixels for specific sample
    pixel_ind = which(pData(MSIobject)$sample_ID == sample)
    IS_ind = which(fData(MSIobject)$analyte == "IS")
    
    # Save IS intensity vector
    IS_vec = spectra(MSIobject)[IS_ind, pixel_ind]
    
    #Deal with 0 values in IS vector - 10% of lowest
    IS_vec[which(IS_vec ==0)] = NA
    mv_impute = min(IS_vec, na.rm = T) / 10
    IS_vec[which(is.na(IS_vec))] = mv_impute
    
    for(mz_ind in 1:nrow(fData(MSIobject))){
      
      print(sprintf("mz - %s", mz_ind))
      
      ints = spectra(MSIobject)[mz_ind, pixel_ind]
      response = ints / IS_vec
      spectra(MSIobject)[mz_ind, pixel_ind] = response
      
    }
  }
  return(MSIobject)
}

MSIobject = msi_response

summarise_cal_levels <- function(MSIobject,
                                 cal_label = "Calibration"){
  
  # Create empty df to add to
  df = data.frame(matrix(ncol = length(unique(pixel_data$sample_ID)),
                  nrow = nrow(fData(MSIobject))))
  colnames(df) = unique(pixel_data$sample_ID)
  rownames(df) = fData(MSIobject)@mz
  
  # create pixel data to associate pixel indices to cal levels
  pixel_data = data.frame(pData(MSIobject)) %>%
    tibble::rownames_to_column("pixel_ind") %>%
    subset(sample_type == cal_label) %>%
    mutate(sample_ID = sprintf("%s_%s", sample_ID, ROI))
     
  for(sample in unique(pixel_data$sample_ID)){
    inds = which(pixel_data$sample_ID == sample)
    
    pixels = as.numeric(pixel_data$pixel_ind[inds])
    
    mean_int_vec = c()
    
    for(mz_ind in 1:nrow(fData(MSIobject))){
      
      ints = spectra(MSIobject)[mz_ind, pixels]
      ints = replace(ints, ints ==0, NA)
      
      mean_int_vec = c(mean_int_vec, mean(ints, na.rm=T))
    }
    
    df[[sample]] = mean_int_vec
    
  }
  
  return(df)
}


int2conc = function(MSIobject,
                    cal_label = "Calibration",
                    cal_list){
  
  pixel_inds = which(pData(MSIobject)$sample_type != cal_label)
  MSIobject = MSIobject[, pixel_inds]
  
  no_cal_indices = c()
  for(i in 1:length(mz(MSIobject))){
    
    print(i)
    
    #if(round(as.numeric(mz(MSIobject)[i], 3)) != round(as.numeric(names(cal_list)[i]), 3)){
    #  return("mz in cal not match mz from MSIobject")
    #}
    
    image(MSIobject, mz = mz(MSIobject)[i], contrast.enhance="histogram",
          superpose = FALSE)
    
    eqn = cal_list[[i]]
    
    if(typeof(eqn) != "list"){
      
      print("no calibration curve")
      no_cal_indices = c(no_cal_indices, i)
      
    } else{
      
      tissue_ints = spectra(MSIobject)[i,]
      tissue_concs = as.vector(sapply(tissue_ints, function(y) inverse.predict(eqn, y)$Prediction))
      
      spectra(MSIobject)[i,] = tissue_concs
      
    }
  }
  
  MSIobject = MSIobject[-no_cal_indices, ]
}

#' Function to create data matrix from MSI object
#'
#' @param MSIobject MSI object from Cardinal
#' @param inputNA Whether to convert 0's in matrix to NA (default = T)
#' @return dataframe with rows = ? cols = ?
#'
#' @export

createDatamatrix <- function(MSIobject, inputNA = T){
  
  df = data.frame(matrix(ncol = length(unique(pixel_df$sample_ID)),
                         nrow = nrow(fData(MSIobject))))
  colnames(df) = unique(pixel_df$sample_ID)
  rownames(df) = fData(MSIobject)@mz
  
  MSIobject = as(MSIobject, "MSContinuousImagingExperiment")
  
  MSIobject = MSIobject[, - which(is.na(pData(MSIobject)$ROI))]
  
  pixel_df = data.frame(pData(MSIobject)) %>%
    mutate(sample_ID = sprintf("%s_%s", ROI, replicate)) %>%
    tibble::rownames_to_column("pixel_ind")
  
  for(roi in unique(pixel_df$sample_ID)){
    
    pixel_inds = as.numeric(pixel_df$pixel_ind[which(pixel_df$sample_ID == roi)])
    
    temp_df = data.frame(spectra(MSIobject)[, pixel_inds])
    rowMeans(temp_df, na.rm=T)
    
    df[[roi]] = rowMeans(temp_df, na.rm=T)
    
  }
  
  if(inputNA){
    df <- replace(df, df==0, NA)
  }
  
  return(list(df, pixel_df, data.frame(fData(MSIobject))))
}
