# load biodivMapR and useful libraries 
library(biodivMapR)
library(utils)
library(stars)
library(sen2r)


# url for the S2 subset
#url <- 'https://gitlab.com/jbferet/myshareddata/-/raw/master/biodivMapR_S2_Sample/RASTER/S2A_T33NUD_20180104_Subset'
# create a temporary directory (choose your own data directory)
tmpdir <- "/home/pndb/Scriptr_MJ/data/L1A/L1A/parametres"
# name your binary raster with the same name as the online file
NameRaster <- '03022021_L1A_Cloud_SatNdx_SNR_SRTM_Clipped.csv'
destfile <- file.path(tmpdir,NameRaster,fsep = '/')
#download.file(url = url, destfile = destfile, method = 'auto', quiet = FALSE, mode = "wb")

# url for the S2 subset header
#urlhdr <-  'https://gitlab.com/jbferet/myshareddata/-/raw/master/biodivMapR_S2_Sample/RASTER/S2A_T33NUD_20180104_Subset.hdr'
# name your raster HDR with the same name as the binary raster, with .hdr extension
destfile_HDR <- get_HDR_name(destfile,showWarnings = FALSE)
#download.file(url = urlhdr, destfile = destfile_HDR, method = 'auto', quiet = FALSE, mode = "w")

# read ENVI file with stars
#Stars_S2 <- read_stars(destfile, along = 'band',proxy = FALSE)
# write it as a tiff image
# create a specific directory for the tiff image and name your raster

r <- file.path(tmpdir,'03022021_L1A_Cloud_SatNdx_SNR_SRTM_Clipped.tif',fsep = '/')
#dir.create(desttiff,showWarnings = FALSE)
#destfiletiff <- file.path(desttiff,'S2_Subset.tif',fsep = '/')
#r <- write_stars(Stars_S2, dsn=destfiletiff, driver =  'GTiff', type='Int16')

# read ENVI file with stars
create_hdr(ImPath = r, Sensor = 'SENTINEL_2A', 
           SpectralBands = NULL, BandName = NULL, WLunits = NULL)

# read ENVI file with stars
BandName <- c('band_02', 'band 03', 'band_04', 'band_05', 'band_06', 
              'band_07', 'band_08', 'band_08A', 'band_11', 'band_12')
SpectralBands <- c(496.6, 560.0, 664.5, 703.9, 740.2, 
                   782.5, 835.1, 864.8, 1613.7, 2202.4)
WLunits <- 'Nanometers'
create_hdr(ImPath = r, Sensor = 'MyOwnSensor', 
           SpectralBands = SpectralBands,BandName = BandName, WLunits = WLunits)

# library
#library(zip)
# name zip file including plots located on the tile
#destzip <- file.path(tmpdir,'S2A_T33NUD_Plots.zip',fsep = '/')
# url for the zip file
#url <- 'https://gitlab.com/jbferet/myshareddata/-/raw/master/biodivMapR_S2_Sample/VECTOR/S2A_T33NUD_Plots.zip'
#download.file(url = url, destfile = destzip)
#destunz <- file.path(tmpdir,'S2A_T33NUD_Plots',fsep = '/')
#unzip(zipfile = destzip,exdir = destunz)
