#Rscript

###########################################
##    Mapping alpha and beta diversity   ##
###########################################

#####Packages :  expint,
#                pracma,
#                utils,
#                raster,
#                dplyr,
#                liquidSVM,
#                progress,
#                matrixStats,
#                ggplot2,
#                prospect,
#                expandFunctions,
#                stringr,
#                XML,
#                rgdal,
#                stars,
#                simsalapar
#                devtools
#install.packages('terra', repos='https://rspatial.r-universe.dev', dependencies = T)
#install.packages("liquidSVM")
#install.packages("expandFunctions")


#devtools::install_gitlab('jbferet/prospect')
#devtools::install_gitlab('jbferet/prosail')

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

}

########################################################################
##                  COMPUTE SPECTRAL INDEX : NDVI                     ##
########################################################################
# Read raster

Refl <- raster::brick(data_raster)
# get raster band name and clean format. Expecting band name and wavelength to be documented in image
HDR_Refl <- read_ENVI_header(get_HDR_name(data_raster))
SensorBands <- HDR_Refl$wavelength
# compute a set of spectral indices defined by IndexList from S2 data
IndexList <- c('NDVI')
# ReflFactor = 10000 when reflectance is coded as INT16
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

write.table(spec_indices, file = "BiodivIndex.tabular", sep = "\t", dec = ".", na = " ", row.names = F, col.names = T, quote = FALSE)

# Update Cloud mask based on radiometric filtering
# eliminate pixels with NDVI < NDVI_Thresh because not enough vegetation
##NDVI_Thresh <- 0.5
##Elim <- which(values(Indices$SpectralIndices[['NDVI']]) < NDVI_Thresh)
##CloudInit <- stars::read_stars(cloudmasks$BinaryMask)
##CloudInit$CloudMask_Binary[Elim] <- 0
# save updated cloud mask
##Cloud_File <- file.path(Cloud_path, 'CloudMask_Binary_Update')
##stars::write_stars(CloudInit, dsn = Cloud_File, driver = "ENVI", type = 'Byte')

ndvi_raster <- Indices$SpectralIndices[[SpIndx]]
ndvi_plot <- rasterVis::levelplot(ndvi_raster, layout = c(0,1,1), main = "NDVI")
ndvi_plot


########################################################################
##                   COMPUTE BIODIVERSITY INDICES                     ##
########################################################################
copNDVI <- raster::raster(data_raster)
copNDVIlr <- raster::reclassify(copNDVI, cbind(252, 255, NA), right=TRUE)

#Resample using raster::aggregate and a linear factor of 10
copNDVIlr <- raster::aggregate(copNDVIlr, fact = 20)
#Set float numbers as integers to further speed up the calculation
storage.mode(copNDVIlr[]) = "integer"

#levelplot(copNDVI,layout=c(0,1,1), main="NDVI 21st of June 1999-2017 - ~8km pixel resolution")
#levelplot(copNDVIlr,layout=c(0,1,1), main="NDVI 21st of June 1999-2017 - ~150km pixel resolution")

#Shannon's Diversity
sha <- rasterdiv::Shannon(copNDVIlr, window = 9, na.tolerance = 0.1, np = 1)
sha_plot <- rasterVis::levelplot(sha, layout = c(0, 1, 1), main = "sha")
sha_plot

#Pielou's Evenness
pie <- rasterdiv::Pielou(copNDVIlr, window = 9, na.tolerance = 0.1, np = 1)
pie_plot <- rasterVis::levelplot(pie,layout = c(0, 1, 1), main = "pie")
pie_plot

#Berger-Parker's Index
ber <- rasterdiv::BergerParker(copNDVIlr, window = 9, na.tolerance = 0.1, np = 1)
ber_plot <- rasterVis::levelplot(ber, layout = c(0, 1, 1), main = "ber")
ber_plot

#Rao's quadratic Entropy
#rao <- rasterdiv::Rao(copNDVIlr,window=2,na.tolerance=0.1,dist_m="euclidean",shannon=FALSE,np=1)
#rao_plot <- rasterVis::levelplot(rao,layout=c(0,1,1), main="rao")

#Parametric Rao's quadratic entropy with alpha ranging from 1 to 5
prao <- rasterdiv::paRao(ndvi_raster, window = 19, alpha = 1:5, na.tolerance = 0.1, dist_m = "euclidean", np = 1)
prao_plot <- rasterVis::levelplot(prao$window.19[[1]], layout = c(0, 1, 1), main = "prao")
prao_plot

#Cumulative Residual Entropy
cre <- rasterdiv::CRE(ndvi_raster, window = 9, na.tolerance = 0.1, np = 1)
cre_plot <- rasterVis::levelplot(cre, layout = c(0, 1, 1), main = "cre")
cre_plot

#Hill's numbers
hil <- rasterdiv::Hill(ndvi_raster, window = 9, alpha = seq(0, 2, 0.5), na.tolerance = 0.1, np = 1)
hil_plot <- rasterVis::levelplot(hil[[3]], layout = c(0, 1, 1), main = "hil")
hil_plot

#Renyi's Index
ren <- rasterdiv::Renyi(ndvi_raster, window = 9, alpha = seq(0, 2, 0.5), na.tolerance = 0.1, np = 1)
ren_plot <- rasterVis::levelplot(ren[[3]], layout = c(0, 1, 1), main = "ren")
ren_plot












