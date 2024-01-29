
library(quantMSImageR)
library(Cardinal)
library(dplyr)
library(chemCal)
library(openxlsx)

# Set paths
base_path = "C:/Users/matsmi/OneDrive - Karolinska Institutet/Dokument/MSI/quantMSImageR"
data_path = "C:/Users/matsmi/OneDrive - Karolinska Institutet/Dokument/MSI/quantMSImageR/data"
out_path = "C:/Users/matsmi/OneDrive - Karolinska Institutet/Dokument/MSI/quantMSImageR/results"

lib_ion_path = sprintf("%s/ion_library.txt", data_path)

fns = c("23Jan_lung_01-LHS-Acq01", "23Jan_lung_01-RHS_Acq02", "23Jan_lung_02-LHS_Acq03", "23Jan_lung_02-RHS_Acq04")

fn = fns[4]


#tissue = read_mrm(name = fn, folder = data_path, lib_ion_path = lib_ion_path , polarity = "Negative")

#tissue_pixels = selectROI(tissue, contrast.enhance="histogram", feature= 1)
#tissue = tissue[, tissue_pixels]

# Save and load Cardinal::image
#save.image(file = sprintf("%s/rdata/%s.RData", out_path, fn))
load(file = sprintf("%s/rdata/%s.RData", out_path, fn))


DT::datatable(data.frame(fData(tissue)))
write.csv(data.frame(fData(tissue)), sprintf("%s/fdata/%s.csv", out_path, fn), row.names=F)


#gc()
#image(tissue, contrast.enhance="histogram", smooth.image = 'gaussian', superpose = FALSE, normalize.image = "linear", feature = 1)
#gc()
#image(tissue, contrast.enhance="histogram", smooth.image = 'gaussian', superpose = FALSE, normalize.image = "linear", feature = 2)
#gc()
#image(tissue, contrast.enhance="histogram", smooth.image = 'gaussian', superpose = FALSE, normalize.image = "linear", feature = 3)
#gc()
#image(tissue, contrast.enhance="histogram", smooth.image = 'gaussian', superpose = FALSE, normalize.image = "linear", feature = 4)
#gc()
#image(tissue, contrast.enhance="histogram", smooth.image = 'gaussian', superpose = FALSE, normalize.image = "linear", feature = 5)
#gc()
#image(tissue, contrast.enhance="histogram", smooth.image = 'gaussian', superpose = FALSE, normalize.image = "linear", feature = 6)
#gc()

############ Plot (histogram enhance)
gc()
png(sprintf("%s/ion_images/%s_histogram.png", out_path, fn), width = 9, height = 7, units = "in", res = 800)
plot.new()
Cardinal::image(tissue, contrast.enhance="histogram", smooth.image = 'gaussian', superpose = FALSE, normalize.image = "linear", feature =fData(tissue)@mz[1:nrow(fData(tissue))])
dev.off()


############ Plot (suppression enhance)
gc()
png(sprintf("%s/ion_images/%s_suppression.png", out_path, fn), width = 9, height = 7, units = "in", res = 800)
plot.new()
Cardinal::image(tissue, contrast.enhance="suppression", smooth.image = 'gaussian', superpose = FALSE, normalize.image = "linear", feature = fData(tissue)@mz[1:nrow(fData(tissue))])
dev.off()
