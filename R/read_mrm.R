require(Cardinal)

#' Function to create data matrix from MSI object
#' @import Cardinal
#' @import dplyr
#' @include setClasses.R
#'
#' @param name Name of Waters raw imaging file (including the imaging/Analyte 1.txt file after processing)
#' @param folder Location of folder with raw imaging data in
#' @param lib_ion_path Full path to library of transitions for targeted MSI experiments. Must include headers: 'transition_id',	'precursor_mz',	'product_mz',	'collision_eV',	'cone_V',	'Polarity',	'Type'. Where 'Type' is "Analyte" for analytes in tissue and "IS" for internal standards.
#' @param polarity String indicating whether data was collected in 'Positive' or 'Negative' polarity.
#'
#' @return MSIobject with slots updated for i) matrix of average ng/pixel of m/z (rows = m/z and cols = cal level) in tissue ROIs ii) sample/ROI metadata
#'
#' @export read_mrm
read_mrm = function(name, folder, lib_ion_path, overwrite=T){

  # set Imaging folder
  imaging_folder = sprintf("%s/%s.raw/imaging", folder, name)

  # Set rds filename
  rds_fn = sprintf("%s/%s.raw/MSImagingExperiment.rds", folder, name)

  if(file.exists(rds_fn) & overwrite==F){
    # Read file
    out = readRDS(rds_fn)

  } else{

    # Read in ion library - from database
    ion_lib = read.table(file = lib_ion_path, sep = "\t", header = T)

    # Find experimental parameters
    inf_file = sprintf("%s/%s.raw/_extern.inf", folder, name)
    lines <- readLines(inf_file)
    ystep = strsplit(lines[ grep("DesiYStep", lines) ], "\t")[[1]]
    ystep = as.numeric( ystep [length(ystep )] ) * 1000
    polarity = strsplit(lines[ grep("MS1 DC Polarity", lines) ], "\t")[[1]]
    polarity = polarity  [length(polarity)]

    # Experimental metadata
    exp_metadata <- CardinalIO::ImzMeta()
    exp_metadata$pixelSize = ystep

    # Read MRM imaging data
    analyte_fn = list.files(imaging_folder, full.names = T)
    analyte_df = read.table(analyte_fn, fill = TRUE, sep="\t", header=F, blank.lines.skip = T)[-1,]

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

    # Gather info of MRM transitions from the ion library
    ion_lib = ion_lib %>%
      subset(Polarity == polarity) %>%
      dplyr::right_join(y=transitions,  by = c('precursor_mz', 'product_mz'), suffix = c("_name", "_int")) %>%
      mutate(transition_id_name = ifelse(is.na(transition_id_name), transition_id_int, transition_id_name),
             Polarity = ifelse(is.na(Polarity), polarity, Polarity),
             Type = ifelse(is.na(Type), "Unknown", Type)) %>%
      arrange(transition_id_int)

    # feature metadata
    fdata <- MassDataFrame(mz=ion_lib$transition_id_int,
                           analyte = ion_lib$Type,
                           precursor_mz = ion_lib$precursor_mz,
                           product_mz = transitions$product_mz,
                           name = ion_lib$transition_id_name)


    # intensity data
    idata = t(analyte_df[, grep("transition", colnames(analyte_df))])


    out <- MSImagingExperiment(spectraData=idata,
                               featureData=fdata,
                               pixelData=pdata,
                               experimentData = exp_metadata)

    # Save file
    saveRDS(out, file = rds_fn)

  }

  return(out)
}




