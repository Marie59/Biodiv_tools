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
install.packages('terra', repos='https://rspatial.r-universe.dev', dependencies = T)
install.packages("liquidSVM")
#install.packages("expandFunctions")
#install.packages("git2r")
devtools::install_gitlab('jbferet/prospect')
devtools::install_gitlab('jbferet/prosail')

#####Load arguments

args <- commandArgs(trailingOnly = TRUE)

#####Import the S2 data

if (length(args) < 1) {
    stop("This tool needs at least 1 argument")
}else{
    data_raster <- args[1]
    rasterheader <- args[2]
    xmlfile <- args[3]
    data_S2 <- args[4]
    source(args[5])  
    source(args[6])
    source(args[7])  
    source(args[8])  
    source(args[9])  
    source(args[10])
    source(args[11])  
    source(args[12])
    source(args[13])
    source(args[14])  
    source(args[15])
    source(args[16])
    source(args[17])  


}

########################################################################
##                      COMPUTE SPECTRAL INDEX                        ##
########################################################################
# Read raster

Refl <- raster::brick(data_raster)
# get raster band name and clean format. Expecting band name and wavelength to be documented in image
HDR_Refl <- prosail::read_ENVI_header(prosail::get_HDR_name(data_raster))
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

write.table(spec_indices, file = "BiodivIndex.tabular", sep = "\t", dec = ".", na = " ", row.names = F, col.names = F, quote = FALSE)
# Update Cloud mask based on radiometric filtering
# eliminate pixels with NDVI < NDVI_Thresh because not enough vegetation
##NDVI_Thresh <- 0.5
##Elim <- which(values(Indices$SpectralIndices[['NDVI']]) < NDVI_Thresh)
##CloudInit <- stars::read_stars(cloudmasks$BinaryMask)
##CloudInit$CloudMask_Binary[Elim] <- 0
# save updated cloud mask
##Cloud_File <- file.path(Cloud_path, 'CloudMask_Binary_Update')
##stars::write_stars(CloudInit, dsn = Cloud_File, driver = "ENVI", type = 'Byte')



########################################################################
##      COMPUTE BIOPHYSICAL VARIABLES BASED ON PROSAIL INVERSION      ##
########################################################################
# get S2 geometry
# read metadata file from S2 image
S2Geom <- get_S2geometry(MTD_TL_xml = xmlfile)

# Train PROSAIL inversion
minval <- data.frame('CHL'=10,'CAR'=0,'EWT' = 0.005,'ANT' = 0,'LMA' = 0.005,'N' = 1.0,'psoil' = 0.0, 'BROWN'=0.0,
                     'LIDFa' = 30, 'lai' = 0.5,'q'=0.1,'tto' = 0,'tts' = min(S2Geom$SZA), 'psi' = 5)
maxval <- data.frame('CHL'=90,'CAR'=20,'EWT' = 0.04,'ANT' = 3,'LMA' = 0.04,'N' = 2.0, 'psoil' = 1.0, 'BROWN'=0.5,
                     'LIDFa' = 70, 'lai' = 7,'q'=0.25,'tto' = 7,'tts' = max(S2Geom$SZA), 'psi' = 355)

# get sensor response for Sentinel-2
SensorName <- HDR_Refl$`sensor type`
SRF <- GetRadiometry(SensorName,Path_SensorResponse = NULL)
# adjust optical constants from 1nm sampling into spectral S2 spectral sampling
wvl <- SpecPROSPECT$lambda
SpecSensor <- PrepareSensorSimulation(SpecPROSPECT,SpecSOIL,SpecATM,SRF)
SpecPROSPECT_Sensor <- SpecSensor[[1]]$SpecPROSPECT_Sensor
SpecSOIL_Sensor <- SpecSensor[[1]]$SpecSOIL_Sensor
SpecATM_Sensor <- SpecSensor[[1]]$SpecATM_Sensor

# define spectral bands required to train SVR model for each variable
S2BandSelect <- list()
S2BandSelect$CHL <- S2BandSelect$lai <- S2BandSelect$EWT <- S2BandSelect$LMA <- c('B03','B04','B05','B06','B07','B08','B11','B12')
ImgBandNames <- strsplit(HDR_Refl$`band names`,split = ',')[[1]]
# get variable ID for train_prosail_inversion
Bands2Select <- list()
for (bpvar in names(S2BandSelect)){
  Bands2Select[[bpvar]] <- match(S2BandSelect[[bpvar]],ImgBandNames)
}

# define noise level for each variable
NoiseLevel <- list()
NoiseLevel$EWT <- 0.025
NoiseLevel$CHL <- 0.01
NoiseLevel$LMA <- NoiseLevel$lai <- 0.05

# where results will be stored
PROSAIL_ResPath <- file.path(results_site_path,'PRO4SAIL_INVERSION')
dir.create(path = PROSAIL_ResPath,showWarnings = FALSE,recursive = TRUE)

modelSVR <- train_prosail_inversion(minval=minval,maxval=maxval,Parms2Estimate=c('CHL','EWT','LMA','lai'),
                                    Bands2Select=Bands2Select,NoiseLevel=NoiseLevel, SAILversion = '4SAIL',
                                    SpecPROSPECT = SpecPROSPECT_Sensor, SpecSOIL = SpecSOIL_Sensor, SpecATM = SpecATM_Sensor,
                                    Path_Results=PROSAIL_ResPath,nbModels = 10,nbSamples = 1000,FigPlot = FALSE)

# Apply SVR model on Sentinel-2 data
Apply_prosail_inversion(raster_path = Refl_path, HybridModel = modelSVR, PathOut = PROSAIL_ResPath,
                        SelectedBands = S2BandSelect,bandname = ImgBandNames,
                        MaskRaster = Cloud_File, MultiplyingFactor = 10000)


##################Gros doute test
