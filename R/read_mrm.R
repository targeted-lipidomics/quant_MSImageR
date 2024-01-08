read_mrm = function(name,
                    folder){

  imaging_folder = sprintf("%s/%s.raw/imaging", folder, name)

  analyte_fn = list.files(imaging_folder, full.names = T)

  analyte_df = read.table(analyte_fn, skip=2, check.names = F, fill = TRUE, sep="\t", header=F)

  transitions = t(analyte_df[1:3,]) %>%
    `colnames<-`(c("transition_id", "precursor_mz", "product_mz")) %>%
    na.omit() %>% data.frame()

  col_heads = c("ind", "x_loci", "y_loci", sprintf("transition_%s", transitions$transition), "sample", "pixel")

  analyte_df = analyte_df %>%
    `colnames<-`(col_heads) %>%
    filter(!row_number() %in% 1:3) %>%
    mutate(x = NA,
           y = NA)

  for(i in 1:length(unique(analyte_df$y_loci))){
    y_val = unique(analyte_df$y_loci)[i]
    ind = which(analyte_df$y_loci == y_val)
    analyte_df$y[ind] = i
  }

  for(i in 1:length(unique(analyte_df$x_loci))){
    x_val = unique(analyte_df$x_loci)[i]
    ind = which(analyte_df$x_loci == x_val)
    analyte_df$x[ind] = i
  }

  # pixel metadata
  coord <- analyte_df%>%select(x,y)
  run <- factor(rep(name, nrow(coord)))
  pdata <- PositionDataFrame(run=run, coord=coord)

  # feature metadata
  fdata <- MassDataFrame(mz=transitions$transition_id,
                         analyte = "analyte",
                         precursor_mz = transitions$precursor_mz,
                         product_mz = transitions$product_mz)

  # intensity data
  idata = t(analyte_df[, grep("transition", colnames(analyte_df))])


  out <- MSImagingExperiment(imageData= idata,
                             featureData=fdata,
                             pixelData=pdata)
  return(out)
}




