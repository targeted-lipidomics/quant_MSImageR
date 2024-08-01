## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = FALSE)

# set default options
knitr::opts_chunk$set(echo = T,
                      message = FALSE,
                      warning = FALSE,
                      fig.align="center",
                      fig.width = 5,
                      fig.height = 5,
                      dpi = 120,
                      fig.retina = 3)


## ----load-packages, message=FALSE---------------------------------------------

library(quantMSImageR)
library(Cardinal)
library(dplyr)
library(tibble)
library(tidyr)
library(chemCal)
library(ggplot2)
library(ggpmisc)
library(viridis)
library(DT)


## ----set-names, message=FALSE, results = 'hide'-------------------------------

data_path = system.file('extdata', package = 'quantMSImageR') # Path to test data
#data_path = "C:/Users/matsmi/OneDrive - Karolinska Institutet/Dokument/MSI/quantMSImageR/data"

# Set filenames
ion_lib_fn = sprintf("%s/ion_library.txt", data_path)
#ion_lib_fn = "C:/Users/matsmi/OneDrive - Karolinska Institutet/Dokument/MSI/quantMSImageR/data/ion_library.txt"

tissue_fn = "tissue_MRM_data"
cal_fn = "cal_MRM_data"


sample_meta_fn =  sprintf("%s/sample_metadata.csv", data_path)
cal_meta_fn = sprintf("%s/calibration_metadata.csv", data_path)

cal_roi_fn = sprintf("%s/%s.raw/cal_ROIs.csv", data_path, cal_fn)
tissue_pixels_fn = sprintf("%s/%s.raw/tissue_ROIs.csv", data_path, tissue_fn)
tissue_anat_fn = sprintf("%s/%s.raw/anat_ROIs.csv", data_path, tissue_fn)


## ----load-metafiles, message=FALSE--------------------------------------------

sample_metadata = read.csv(sample_meta_fn)

cal_metadata = read.csv(cal_meta_fn)
# cal_metadata requires 'identifier' header which matches the labels from the calibration sample

cal_roi_df = read.csv(cal_roi_fn)
tissue_pixels_df = read.csv(tissue_pixels_fn)
tissue_anat_ROIs = read.csv(tissue_anat_fn)


## ----load-DESI-files, message=FALSE-------------------------------------------

tissue = read_mrm(name = tissue_fn, folder = data_path, lib_ion_path = ion_lib_fn)

# Read calibration mix MRM files
cal_mix =  read_mrm(name = cal_fn, folder = data_path, lib_ion_path = ion_lib_fn)


## ----set-MRM-axes, message=FALSE----------------------------------------------

tissue = setCommonAxis(MSIobjects = list(tissue, cal_mix), ref_fdata = fData(cal_mix))[[1]]

# Tissue feature data
DT::datatable(data.frame(fData(tissue)), width = '75%')

# Calibration miz feature data
DT::datatable(data.frame(fData(cal_mix)), width = '75%')


## ----manual-ROI-select, message=FALSE-----------------------------------------
# Image cal mix
image(cal_mix, enhance= "histogram")

# Example
#cal1_rep1 <- selectROI(cal_mix, contrast.enhance="suppression")


## ----cal-ROI-data, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

# Check ROI info
DT::datatable( head(cal_roi_df) )


## ----cal_pdata, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

cal_level = makeFactor(L1_r1_lipid01 = cal_roi_df$cal1_rep1, L2_r1_lipid01 = cal_roi_df$cal2_rep1,
                       L3_r1_lipid01 = cal_roi_df$cal3_rep1, L4_r1_lipid01 = cal_roi_df$cal3p5_rep1,
                       L5_r1_lipid01 = cal_roi_df$cal4_rep1, L6_r1_lipid01 = cal_roi_df$cal4p5_rep1,
                       L7_r1_lipid01 = cal_roi_df$cal5_rep1,
                       L1_r2_lipid01 = cal_roi_df$cal1_rep2, L2_r2_lipid01 = cal_roi_df$cal2_rep2,
                       L3_r2_lipid01 = cal_roi_df$cal3_rep2, L4_r2_lipid01 = cal_roi_df$cal3p5_rep2,
                       L5_r2_lipid01 = cal_roi_df$cal4_rep2, L6_r2_lipid01 = cal_roi_df$cal4p5_rep2,
                       L7_r2_lipid01 = cal_roi_df$cal5_rep2,
                       L1_r3_lipid01 = cal_roi_df$cal1_rep3, L2_r3_lipid01 = cal_roi_df$cal2_rep3,
                       L3_r3_lipid01 = cal_roi_df$cal3_rep3, L4_r3_lipid01 = cal_roi_df$cal3p5_rep3,
                       L5_r3_lipid01 = cal_roi_df$cal4_rep3, L6_r3_lipid01 = cal_roi_df$cal4p5_rep3,
                       L7_r3_lipid01 = cal_roi_df$cal5_rep3)

#Subset metadata file for calibration DESI-MRM info
sample_metadata_subset = subset(sample_metadata, Sample == cal_fn)

## Update pixel metadata
pData(cal_mix)$sample_type = sample_metadata_subset$Sample_type
pData(cal_mix)$replicate = sample_metadata_subset$Replicate
pData(cal_mix)$identifier = cal_level # match identifier header in cal_metadata
pData(cal_mix)$treatment = sample_metadata_subset$Treatment
pData(cal_mix)$ROI = sapply(as.character(cal_level), FUN = function(x) strsplit(x, "_r")[[1]][1])
pData(cal_mix)$sample_ID = sample_metadata_subset$Sample_ID
pData(cal_mix)$pixel_size = as.numeric(experimentData(cal_mix)$pixelSize)


## ----image-cal-ROIs, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

# Image cal mix labelled
image(cal_mix, "ROI", zlab = "intensity")


## ----cal-pdata-data2, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

# View updated pixel metadata
DT::datatable(head(data.frame(pData(cal_mix))), width = '90%')


## ----write-cal, message=FALSE-------------------------------------------------

# Save cal_mix imzML for long term storage
cal_imzfile <- tempfile(fileext="_cal.imzML")
writeMSIData(cal_mix, file = cal_imzfile)
list.files(cal_imzfile)

# Read imzML if needed
#cal_mix = readMSIData(cal_imzfile)


## ----set-tissue-pixels, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

## Image tissue
image(tissue, enhance= "histogram", zlab = "intensity")

tissue_pixels = tissue_pixels_df$tissue_pixels
noise_pixels = tissue_pixels_df$noise_pixels

noise = tissue[, noise_pixels]
tissue = tissue[, tissue_pixels]


## ----image-tissue-pixels, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----
# Image tissue subset
image(tissue, enhance= "histogram", zlab = "intensity")

## ----set-tissue-ROIs, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

tissue_rois = makeFactor(airway_01 = tissue_anat_ROIs$airway_01,
                         airway_02 = tissue_anat_ROIs$airway_02,
                         airway_03 = tissue_anat_ROIs$airway_03,
                         airway_04 = tissue_anat_ROIs$airway_04,
                         airway_05 = tissue_anat_ROIs$airway_05,
                         airway_06 = tissue_anat_ROIs$airway_06,
                         parenchyma_01 = tissue_anat_ROIs$parenchyma_01,
                         parenchyma_02 = tissue_anat_ROIs$parenchyma_02,
                         parenchyma_03 = tissue_anat_ROIs$parenchyma_03,
                         parenchyma_04 = tissue_anat_ROIs$parenchyma_04,
                         parenchyma_05 = tissue_anat_ROIs$parenchyma_05,
                         parenchyma_06 = tissue_anat_ROIs$parenchyma_06)


## ----tissue_pdata, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

#Subset metadata file for tissue DESI-MRM info
sample_metadata_subset = subset(sample_metadata, Sample == tissue_fn)

## Update pixel metadata
pData(tissue)$sample_type = sample_metadata_subset$Sample_type
pData(tissue)$replicate = sample_metadata_subset$Replicate
pData(tissue)$identifier = tissue_rois 
pData(tissue)$treatment = sample_metadata_subset$Treatment
pData(tissue)$ROI = sapply(as.character(tissue_rois), FUN = function(x) strsplit(x, "_")[[1]][1])
pData(tissue)$sample_ID = sample_metadata_subset$Sample_ID
pData(tissue)$pixel_size = as.numeric(experimentData(tissue)$pixelSize)


## ----image-tissue-ROIs, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

# Image cal mix labelled
image(tissue, "ROI", zlab = "intensity")


## ----tissue-pdata-data2, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

# View updated pixel metadata
DT::datatable(head(data.frame(pData(tissue))))


## ----noise-pdata, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

sample_metadata_subset = subset(sample_metadata, Sample == tissue_fn)

pData(noise)$sample_type = "Noise"
pData(noise)$replicate = sample_metadata_subset$Replicate
pData(noise)$identifier = "noise"
pData(noise)$treatment = sample_metadata_subset$Treatment
pData(noise)$ROI = "noise"
pData(noise)$sample_ID = sample_metadata_subset$Sample_ID
pData(noise)$pixel_size = as.numeric(experimentData(noise)$pixelSize)


## ----image-noise-ROI, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

# Image tissue labelled
image(noise, enhance= "histogram", zlab = "intensity")


## ----noise-pdata2, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

# View updated pixel metadata
DT::datatable(head(data.frame(pData(noise))))


## ----write-tissue, message=FALSE----------------------------------------------

# Save tissue imzML for long term storage
writeMSIData(tissue, file = sprintf("%s/%s.raw/tissue.imzML", data_path, tissue_fn))
#tissue = readMSIData(file = sprintf("%s/%s.raw/tissue.imzML", data_path, tissue_fn))

# Save noise imzML for long term storage
writeMSIData(noise, file = sprintf("%s/%s.raw/noise.imzML", data_path, tissue_fn))
#noise = readMSIData(file = sprintf("%s/%s.raw/noise.imzML", data_path, tissue_fn))


## ----merge-data, message=FALSE, fig.align='center', fig.show='hold', out.width="85%",out.height="49%"----

msi_combined = combine_MSIs(cal_mix, tissue, noise)
msi_combined

image(msi_combined, enhance= "histogram", zlab = "intensity")

# Save combined imzML for long term storage
writeMSIData(msi_combined, file = sprintf("%s/%s.raw/combined.imzML", data_path, tissue_fn))
#tissue = readMSIData(file = sprintf("%s/%s.raw/combined.imzML", data_path, tissue_fn))


## ----set-NA, message=FALSE, fig.align='center', fig.show='hold', out.width="90%",out.height="49%"----

msi_combined = zero2na(MSIobject = msi_combined, val_slot = "intensity")
msi_combined

#image(msi_combined, enhance= "histogram", zlab = "intensity")


## ----remove-MRMs, message=FALSE, fig.align='center', fig.show='hold', out.width="90%",out.height="49%"----

msi_combined = remove_blank_mzs(MSIobject = msi_combined)
msi_combined

#image(msi_combined, enhance= "histogram", zlab = "intensity")


## ----normalise, message=FALSE, fig.align='center', fig.show='hold', out.width="90%",out.height="49%"----

msi_combined = int2response(MSIobject = msi_combined, val_slot = "intensity")
msi_combined

#image(msi_combined, response ~ x * y, enhance= "histogram", zlab = "response")


## ----calc-SNR, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

msi_combined = int2snr(MSIobject = msi_combined, val_slot = "response",
                       noise = "Noise", tissue = "Tissue", snr_thresh = 3)
msi_combined

#image(msi_combined, snr ~ x * y, enhance= "histogram", zlab = "snr")


## ----summarise-cal, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

msi_combined = summarise_cal_levels(MSIobject = msi_combined,
                                    cal_metadata = cal_metadata,
                                    val_slot = "response",
                                    cal_label = "Cal",
                                    id = "identifier")

# Number of pixels per calibration spot
DT::datatable(data.frame(msi_combined@calibrationInfo@cal_response_data))


## ----cal-curves, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

msi_combined = create_cal_curve(MSIobject = msi_combined,
                                cal_type = "Cal")


## ----check-cal, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

# Plot linear model
d = msi_combined@calibrationInfo@cal_response_data %>%
  mutate(labs = sapply(cal_spot, FUN=function(x){
    return( paste(strsplit(x, "_")[[1]][3:2], collapse="_") )})) %>%
  subset(!labs %in% c("L7_rep1", "L7_rep2"))

p = ggplot(d, aes(x=pg_perpixel, y = response_perpixel, label =labs)) +
  geom_smooth(method='lm', formula= y~x, col = "red", se=F, size = 1, linetype = "dashed") +
  geom_point(size = 2) +
  theme_Publication() +
  geom_text(hjust=-0.1, vjust=-0.1) +
  labs(x = "pg / pixel")
p


## ----check-cal2, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----
# linear model equation
msi_combined@calibrationInfo@cal_list

# r2 values for each calibration
DT::datatable(data.frame(msi_combined@calibrationInfo@r2_df))

## ----quant-tissue, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

msi_combined = int2conc(MSIobject = msi_combined,
                        val_slot = "response",
                        pixels = c("Tissue", "Noise"))


## ----save-tissue, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----
# Save the combined DESI-MRM object as RDS - enables multiple spectra data channels to be stored in single object
rdsfile <- tempfile(fileext=".RDS")

saveRDS(msi_combined, file = sprintf("%s/%s.raw/combined_data.RDS", data_path, tissue_fn))

list.files(cal_imzfile)

#readRDS(sprintf("%s/%s.raw/combined_data.RDS", data_path, tissue_fn))


## ----image-int, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

msi_tissue = msi_combined[, which(pData(msi_combined)$sample_type == "Tissue")]

# image intensity of 12,13-DiHOME across tissue
image(msi_tissue, intensity ~ x * y, enhance= "histogram", zlab = "intensity")

# Save imzML for long term storage
imzfile <- tempfile(fileext="_tissue_intensity.imzML")
writeMSIData(msi_tissue, intensity=spectra(msi_tissue, "intensity"), file = imzfile)
list.files(imzfile)


## ----image-response, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----
# image response of 12,13-DiHOME across tissue
image(msi_tissue, response ~ x * y, enhance= "histogram", zlab = "response")

# Save imzML for long term storage
imzfile <- tempfile(fileext="_tissue_response.imzML")
writeMSIData(msi_tissue, intensity=spectra(msi_tissue, "response"), file = imzfile)
list.files(imzfile)

## ----image-SNR, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

# image signla-noise-ratio of 12,13-DiHOME across tissue
image(msi_tissue, snr ~ x * y, enhance= "histogram", zlab = "SNR")

# Save imzML for long term storage
imzfile <- tempfile(fileext="_tissue_snr.imzML")
writeMSIData(msi_tissue, intensity=spectra(msi_tissue, "snr"), file = imzfile)
list.files(imzfile)

## ----image-conc, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

msi_tissue_conc = msi_tissue[, which(spectra(msi_combined,"snr")[1,] >= 5)]

# image conc (pg/pixel) of 12,13-DiHOME across tissue
image(msi_tissue_conc, `conc - pg/pixel` ~ x * y, enhance= "histogram", zlab = "conc - pg/pixel")

# Save imzML for long term storage
imzfile <- tempfile(fileext="_tissue__pg_pixel.imzML")
writeMSIData(msi_tissue, intensity=spectra(msi_tissue, "conc - pg/pixel"), file = imzfile)
list.files(imzfile)


## ----image-conc2, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----
# image conc (pg/pixel) of 12,13-DiHOME across tissue
image(msi_tissue_conc, `conc - pg/mm2` ~ x * y, enhance= "histogram", zlab = "conc - pg/mm2")

# Save imzML for long term storage
imzfile <- tempfile(fileext="_tissue_pg_mm2.imzML")
writeMSIData(msi_tissue, intensity=spectra(msi_tissue, "conc - pg/mm2"), file = imzfile)
list.files(imzfile)


## ----image-conc-quant, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

imageR(msi_tissue,
           val_slot = "conc - pg/pixel",
           value = "pg / pixel",
           scale = "suppress", # "suppress" "histogram"
           threshold = 5,
           sample_lab = "sample_ID",
           pixels = "Tissue",
           percentile=99.5,
           overlay = F,
           feat_ind = 1,
           blank_back = F)


## ----image-conc-quant2, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----
imageR(msi_tissue,
           val_slot = "conc - pg/mm2",
           value = "pg / mm2",
           scale = "suppress", # "suppress" "histogram"
           threshold = 5,
           sample_lab = "sample_ID",
           pixels = "Tissue",
           percentile=99.5,
           overlay = F,
           feat_ind = 1,
           blank_back = F)


## ----write-rds, message=FALSE-------------------------------------------------

# Save cal_mix imzML for long term storage
rdsfile <- tempfile(fileext=".RDS")

saveRDS(msi_combined, file = sprintf("%s/%s.raw/combined_data.RDS", data_path, tissue_fn))

list.files(cal_imzfile)

# Read RDS if needed
#readRDS(sprintf("%s/%s.raw/combined_data.RDS", data_path, tissue_fn))


## ----msi-matrices, message=FALSE----------------------------------------------
# Extract the average amount per pixel at each ROI in tissue (stored in the S4 object)
msi_tissue <- createMSIDatamatrix(MSIobject = msi_tissue, val_slot = "response", roi_header = "identifier")

# Data matrix with average concentration over ROI for each lipid
roi_average_matrix = msi_tissue@tissueInfo@roi_average_matrix
DT::datatable(roi_average_matrix, width = '50%')

# Data matrix with concentration per pixel for entire sample for each lipid
all_pixel_matrix = msi_tissue@tissueInfo@all_pixel_matrix
DT::datatable(head(all_pixel_matrix), width = '50%')

# Data pertaining to sample and pixel (including ROI) metadata
sample_metadata = msi_tissue@tissueInfo@sample_metadata
DT::datatable(head(sample_metadata))

# Data pertaining to feature metadata
feature_metadata = data.frame(fData(msi_tissue))
DT::datatable(feature_metadata, width = '75%')


## ----uva-ave, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

image(tissue, "identifier", zlab = "intensity")

uva_df_ave = roi_average_matrix %>%
  tibble::rownames_to_column("roi_lab") %>%
  tidyr::pivot_longer(cols = -roi_lab, names_to = "lipid", values_to = "int") %>%
  dplyr::left_join(y= distinct(sample_metadata %>% select(any_of(c("identifier", "ROI")))), by = c("roi_lab" = "identifier"), keep=F) %>%
  rename("tissue_type" = "ROI")

p_val_ave = t.test(uva_df_ave$int[uva_df_ave$tissue_type == "airway"], uva_df_ave$int[uva_df_ave$tissue_type == "parenchyma"])[["p.value"]]


## ----boxplot-ave, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

ggplot(uva_df_ave, aes(x=tissue_type, col=tissue_type, y = int)) + 
  theme_classic(base_size = 18) +
  scale_colour_Publication() +
  geom_boxplot(alpha = 0.1, size = 1) +
  geom_point(size = 2.5, alpha = 0.2) +
  labs(x="Tissue type", y = "DESI-MRM response") + 
  theme(legend.position="bottom",
        axis.text.x = element_blank(), #element_text(size=11, face="bold"),
        axis.text.y = element_text(size=11, face="bold"),
        axis.title.x =  element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y =  element_text(size=11, face="bold"),
        strip.text = element_text(size=15, face="bold")) +
  facet_wrap(~lipid, scale = "free", nrow=1)



## ----uva-all, message=FALSE, fig.align='center', fig.show='hold', out.width="49%",out.height="49%"----

image(tissue, "ROI", zlab = "intensity")

uva_df_all = all_pixel_matrix %>%
  tibble::rownames_to_column("pixel_ind") %>%
  tidyr::pivot_longer(cols = -pixel_ind, names_to = "lipid", values_to = "int") %>%
  dplyr::left_join(y= distinct(sample_metadata %>% select(any_of(c("pixel_ind", "ROI")))), by = "pixel_ind", keep=F) %>%
  rename("tissue_type" = "ROI") %>%
  group_by(tissue_type) %>%
  mutate(num = n()) %>%
  slice_sample(n = 300)

p_val_all = t.test(uva_df_all$int[uva_df_all$tissue_type == "airway"], uva_df_all$int[uva_df_all$tissue_type == "parenchyma"])[["p.value"]]


## ----boxplot-all, message=FALSE-----------------------------------------------

ggplot(uva_df_all, aes(x=tissue_type, col=tissue_type, y = int)) + 
  theme_classic(base_size = 18) +
  scale_colour_Publication() +
  geom_boxplot(alpha = 0.1, size = 1, outliers=FALSE) +
  labs(x="Tissue type", y = "DESI-MRM response") + 
  theme(legend.position="bottom",
        axis.text.x = element_blank(), #element_text(size=11, face="bold"),
        axis.text.y = element_text(size=11, face="bold"),
        axis.title.x =  element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y =  element_text(size=11, face="bold"),
        strip.text = element_text(size=15, face="bold")) +
  facet_wrap(~lipid, scale = "free", nrow=1)


## ----message=FALSE------------------------------------------------------------

sessionInfo()


