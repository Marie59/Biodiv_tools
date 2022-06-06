#Rscript

###########################################
##    Mapping alpha and beta diversity   ##
###########################################

#####Packages : sen2r,
#    raster,
#    sf,
#    rgeos,
#    sp,
#    raster,
#    tools,
#    stars,
#    stringr,
#    progress,
#    methods,
#    utils,
#    rgdal,
#    R.utils,
#    gdalUtils,
#    fasterize,
#    XML

#####Load arguments

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    stop("This tool needs at least 1 argument")
}else{
    data <- args[1]
    source(args[2])
    data_source <- as.character(args[3])
    id <- as.character(args[4])
    pwd <- as.character(args[5])
}   

###Identification to acces scihub with sen2r package
sen2r::write_scihub_login(id, pwd)

##____________________________________________________________________##
##        Define where data is stored and where to write results      ##
##--------------------------------------------------------------------##

#Create a directory where to unzip your folder of data
dir.create("data_dir")
unzip(data, exdir = "data_dir")

# define raster path
Path_S2 <- file.path("data_dir",list.files("data_dir",pattern = '.SAFE'))

# Result directory
result_path <- "results"
dir.create(path = result_path,showWarnings = FALSE,recursive = TRUE)

##____________________________________________________________________##
##                  Extract, resample & stack data                    ##
##--------------------------------------------------------------------##
# define resolution
resolution <- 10
# define source of data
S2source <- data_source

S2obj <- extract_from_S2_L2A(Path_dir_S2 = Path_S2,
                                        path_vector = NULL,
                                        S2source = S2source,
                                        resolution = resolution)

# update shapefile if needed (reprojection)
###path_vector <- S2obj$path_vector

# create specific result directory corresponding to granule name
#results_site_path <- file.path(result_path,basename(S2obj$S2_Bands$GRANULE))
#dir.create(path = results_site_path, showWarnings = FALSE, recursive = TRUE)
##____________________________________________________________________##
##                        Write CLOUD MASK                            ##
##--------------------------------------------------------------------##

# directory for cloud mask
Cloud_path <- file.path(result_path,'CloudMask')
dir.create(path = Cloud_path, showWarnings = FALSE, recursive = TRUE)
# Filename for cloud mask
cloudmasks <- save_cloud_s2(S2_stars = S2obj$S2_Stack,
                                       Cloud_path = Cloud_path,
                                       S2source = S2source, SaveRaw = T)

##____________________________________________________________________##
##                        Write REFLECTANCE                           ##
##--------------------------------------------------------------------##
 
# directory for Reflectance
Refl_dir <- file.path(result_path,'Reflectance')
dir.create(path = Refl_dir, showWarnings = FALSE, recursive = TRUE)
# filename for Reflectance
Refl_path <- file.path(Refl_dir, paste(basename(S2obj$S2_Bands$GRANULE), '_Refl', sep = ''))

# Save Reflectance file as ENVI image with BIL interleaves
tile_S2 <- substring(strsplit(basename(S2obj$S2_Bands$GRANULE),'_')[[1]][2],2)
dateAcq_S2 <- as.Date(substring(strsplit(basename(S2obj$S2_Bands$GRANULE),'_')[[1]][4],1,8),format="%Y%m%d")

save_data <- save_reflectance_s2(S2_stars = S2obj$S2_Stack, Refl_path = Refl_path,
                               S2Sat = NULL, tile_S2 = tile_S2, dateAcq_S2 = dateAcq_S2,
                               Format = 'ENVI', datatype = 'Int16', MTD = S2obj$S2_Bands$metadata, MTD_MSI = S2obj$S2_Bands$metadata_MSI)



