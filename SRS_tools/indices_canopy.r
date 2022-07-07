#Rscript

###########################################
##    Mapping alpha and beta diversity   ##
###########################################

#####Packages :  expint,
#                pracma,
#                R.utils,
#                raster,
#                sp,
#                matrixStats,
#                ggplot2,
#                expandFunctions,
#                stringr,
#                XML,
#                rgdal,
#                stars,
#####Load arguments

args <- commandArgs(trailingOnly = TRUE)

#####Import the S2 data

if (length(args) < 1) {
    stop("This tool needs at least 1 argument")
}else{
    data_raster <- args[1]
    rasterheader <- args[2]
    source(args[3])
    source(args[4])
    source(args[5])
    indice_choice <- as.character(args[6])

}

########################################################################
##                  COMPUTE SPECTRAL INDEX : NDVI                     ##
########################################################################
# Read raster
nv_data <- raster::raster(data_raster)
nv_data <- raster::aggregate(nv_data, fact = 50)
# Convert raster to SpatialPointsDataFrame
r_pts <- raster::rasterToPoints(nv_data, spatial = T)

# reproject sp object
geo.prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" 
r_pts <- sp::spTransform(r_pts, sp::CRS(geo.prj)) 


# Assign coordinates to @data slot, display first 6 rows of data.frame
r_pts@data <- data.frame(r_pts@data, longitude = sp::coordinates(r_pts)[, 1],
                         latitude = sp::coordinates(r_pts)[, 2])                         


Refl <- raster::brick(data_raster)
# get raster band name and clean format. Expecting band name and wavelength to be documented in image
HDR_Refl <- read_ENVI_header(get_HDR_name(data_raster))
SensorBands <- HDR_Refl$wavelength
# compute a set of spectral indices defined by IndexList from S2 data
IndexList <- c(indice_choice)
# ReflFactor = 10000 when reflectance is coded as INT16
Refl <- raster::aggregate(Refl, fact = 50)
Indices <- ComputeSpectralIndices_Raster(Refl = Refl, SensorBands = SensorBands,
                                                  Sel_Indices = IndexList,
                                                  ReflFactor = 10000, StackOut=F)

# create directory for Spectral indices
results_site_path <- "RESULTS"
SI_path <- file.path(results_site_path, 'SpectralIndices')
dir.create(path = SI_path, showWarnings = FALSE, recursive = TRUE)
# Save spectral indices
for (SpIndx in names(Indices$SpectralIndices)) {
  Index_Path <- file.path(SI_path, paste(basename(data_raster), '_', SpIndx, sep = ''))
  spec_indices <- stars::write_stars(stars::st_as_stars(Indices$SpectralIndices[[SpIndx]]), dsn = Index_Path, driver = "ENVI", type = 'Float32')
  # write band name in HDR
  HDR <- read_ENVI_header(get_HDR_name(Index_Path))
  HDR$`band names` <- SpIndx
  HDR_name <- write_ENVI_header(HDR = HDR, HDRpath = get_HDR_name(Index_Path))
}

spec_indices <- as.data.frame(spec_indices)
r_pts@data[, indice_choice] <- spec_indices[, 3]
write.table(r_pts@data, file = "Spec_Index.tabular", sep = "\t", dec = ".", na = " ", row.names = F, col.names = T, quote = FALSE)

# Update Cloud mask based on radiometric filtering
# eliminate pixels with NDVI < NDVI_Thresh because not enough vegetation
##NDVI_Thresh <- 0.5
##Elim <- which(values(Indices$SpectralIndices[['NDVI']]) < NDVI_Thresh)
##CloudInit <- stars::read_stars(cloudmasks$BinaryMask)
##CloudInit$CloudMask_Binary[Elim] <- 0
# save updated cloud mask
##Cloud_File <- file.path(Cloud_path, 'CloudMask_Binary_Update')
##stars::write_stars(CloudInit, dsn = Cloud_File, driver = "ENVI", type = 'Byte')

spectrale_indices <- function(data, indice_choice) {
  graph_indices <- ggplot2::ggplot(data) +
  ggplot2::geom_point(ggplot2::aes_string(x = data[, 2], y = data[, 3], color = data[, 4]), size = 1, shape = "square") + ggplot2::scale_colour_gradient2(low = "white", high = "#087543", na.value = "grey50") +
  ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), plot.title = ggplot2::element_text(color = "black", size = 12, face = "bold")) + ggplot2::ggtitle(indice_choice)
  
ggplot2::ggsave(paste0(indice_choice, ".png"), graph_indices, width = 12, height = 10, units = "cm")
}

spectrale_indices(data = r_pts@data, indice_choice = indice_choice)

