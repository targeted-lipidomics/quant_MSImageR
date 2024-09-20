
desi_coord_gen = function(fn, type = "DESI", plate_x = 75000, plate_y = 25000){

  desi_anal_fn = sprintf("%s/imaging/Analyte 1.txt", fn)

  desi_anal_df = read.table(desi_anal_fn, fill = TRUE, sep="\t", header=F, blank.lines.skip = T)[-1,1:3] %>%
    `colnames<-`(c("pixel_ind", "x_loci", "y_loci")) %>%
    filter(!row_number() %in% 1:3)



  desi_meta_fn = sprintf("%s/_extern.inf", fn)

  meta_lines <- readLines(desi_meta_fn)

  metadata <- list()

  for (line in meta_lines) {
    # Skip empty lines and comments
    if (line == "" || grepl("^;", line)) next

    # Split the line by the delimiter

    split_line <- strsplit(line, "\t\t")[[1]]

    # Trim whitespace
    key <- trimws(split_line[1])
    value <- trimws(split_line[2])

    # Store the key-value pair in the list
    metadata[[key]] <- value
  }

  x_offset = 1000 * as.numeric( metadata[["DesiXStart"]] )
  y_offset = 1000 * as.numeric( metadata[["DesiYStart"]] )

  pixel_size = as.numeric(metadata[["DesiYStep"]] )

  tibble::tibble(
    image = type,
    desi_position = metadata[["DesiSlot"]],
    origin_x = (plate_x / 2)  + x_offset,
    origin_y = (plate_y / 2) + y_offset,
    x_um = 1000 * as.numeric( metadata[["DesiXLength"]] ),
    y_um = 1000 * as.numeric( metadata[["DesiYLength"]] ),
    x_step = 1000 * pixel_size,
    y_step = 1000 * pixel_size,
    x_pixels = length(unique(desi_anal_df$x_loci)),
    y_pixels = length(unique(desi_anal_df$y_loci)))

}

desi_pixelInfo = desi_coord_gen(fn = sprintf("%s/tissue_MRM_data.raw", system.file('extdata', package = 'quantMSImageR')),
                                type = "DESI", plate_x = 75000, plate_y = 25000)
