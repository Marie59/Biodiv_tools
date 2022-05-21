#Rscript

###########################################
##    Mapping alpha and beta diversity   ##
###########################################

#####Packages : stars
#               utils
#               biodivmapr
#               raster
#               sf
#               mapview
#               leafpop
#               RColorBrewer
#               labdsv
#               rgdal
#               ggplot2
#               gridExtra
remotes::install_github('jbferet/biodivMapR')
#####Load arguments

args <- commandArgs(trailingOnly = TRUE)


# url for the S2 subset

if (length(args) < 1) {
    stop("This tool needs at least 1 argument")
}else{
    data_raster <- args[1]
    rasterheader <- args[2]
    source(args[3])
    source(args[4])
    source(args[5])
    source(args[6])
    source(args[7])
    source(args[8])
    source(args[9])
    source(args[10])
    source(args[11])
    source(args[12])
}   

#####Import data

destfile <- data_raster
tmpdir <- "/home/pndb-cr/Marie/Scriptr_MJ"
# name your binary raster with the same name as the online file
NameRaster <- 'S2A_T33NUD_20180104_Subset'
destfile <- file.path(tmpdir,NameRaster,fsep = '/')

# name your binary raster with the same name as the online file
# name your raster HDR with the same name as the binary raster, with .hdr extension
destfile_HDR <- get_HDR_name(destfile,showWarnings = FALSE)

url <- 'https://gitlab.com/jbferet/myshareddata/-/raw/master/biodivMapR_S2_Sample/RASTER/S2A_T33NUD_20180104_Subset'
download.file(url = url, destfile = destfile, method = 'auto', quiet = FALSE, mode = "wb")

urlhdr <- 'https://gitlab.com/jbferet/myshareddata/-/raw/master/biodivMapR_S2_Sample/RASTER/S2A_T33NUD_20180104_Subset.hdr'
download.file(url = urlhdr, destfile = destfile_HDR, method = 'auto', quiet = FALSE, mode = "w")

# read ENVI file with starss
#Stars_S2 <- stars::read_stars(destfile, along = 'band',proxy = FALSE)
# write it as a tiff image
#destfiletiff <- file.path(tiff_file)
#r <- stars::write_stars(Stars_S2, dsn=destfiletiff, driver =  'GTiff', type='Int16')

# read ENVI file with stars
#biodivmapr::create_hdr(ImPath = destfiletiff, Sensor = 'SENTINEL_2A', 
#           SpectralBands = NULL, BandName = NULL, WLunits = NULL)

# read ENVI file with stars
BandName <- c('band_02', 'band 03', 'band_04', 'band_05', 'band_06', 
              'band_07', 'band_08', 'band_08A', 'band_11', 'band_12')
SpectralBands <- c(496.6, 560.0, 664.5, 703.9, 740.2, 
                   782.5, 835.1, 864.8, 1613.7, 2202.4)
WLunits <- 'Nanometers'
#create_hdr(ImPath = destfiletiff, Sensor = 'MyOwnSensor', 
#           SpectralBands = SpectralBands,BandName = BandName, WLunits = WLunits)



################################################################################
##              DEFINE PARAMETERS FOR DATASET TO BE PROCESSED                 ##
################################################################################
# expected to be in ENVI HDR  

Input_Image_File <- destfile

# path for the Mask raster corresponding to image to process
# expected to be in ENVI HDR format, 1 band, integer 8bits
# expected values in the raster: 0 = masked, 1 = selected
# set to FALSE if no mask available
Input_Mask_File <- FALSE

# relative or absolute path for the Directory where results will be stored
# For each image processed, a subdirectory will be created after its name
Output_Dir  = 'RESULTS'

# SPATIAL RESOLUTION
# resolution of spatial units for alpha and beta diversity maps (in pixels), relative to original image
# if Res.Map = 10 for images with 10 m spatial resolution, then spatial units will be 10 pixels x 10m = 100m x 100m surfaces
# rule of thumb: spatial units between 0.25 and 4 ha usually match with ground data
# too small window_size results in low number of pixels per spatial unit, hence limited range of variation of diversity in the image
window_size <- 10

# PCA FILTERING:        Set to TRUE if you want second filtering based on PCA outliers to be processed. Slower
FilterPCA <- TRUE

# type of PCA:
# PCA: no rescaling of the data
# SPCA: rescaling of the data
TypePCA <-'SPCA'


################################################################################
##                    DEFINE PARAMETERS FOR METHOD                            ##
################################################################################
nbCPU <- 4
MaxRAM <- 0.5
nbclusters <- 50

################################################################################
##                              PROCESS IMAGE                                 ##
################################################################################
# 1- Filter data in order to discard non vegetated / shaded / cloudy pixels
NDVI_Thresh <- 0.5
Blue_Thresh <- 500
NIR_Thresh  <- 1500
Continuum_Removal <- TRUE

print("PERFORM RADIOMETRIC FILTERING")

ImPathShade <- perform_radiometric_filtering(
  Image_Path = Input_Image_File, Mask_Path = Input_Mask_File, Output_Dir = Output_Dir, NDVI_Thresh = NDVI_Thresh, 
  Blue_Thresh = Blue_Thresh, NIR_Thresh = NIR_Thresh)

print("PERFORM PCA ON RASTER")
PCA_Output <- perform_PCA(Input_Image_File = Input_Image_File, Input_Mask_File = Input_Mask_File,
                          Output_Dir = Output_Dir, TypePCA = TypePCA, FilterPCA = FilterPCA, nbCPU = nbCPU, MaxRAM = MaxRAM)
PCA_Files <- PCA_Output$PCA_Files
Pix_Per_Partition <- PCA_Output$Pix_Per_Partition
nb_partitions <- PCA_Output$nb_partitions
Input_Mask_File <- PCA_Output$ImPathShade
PCA_model <- PCA_Output$PCA_model
SpectralFilter <- PCA_Output$SpectralFilter
# path for the updated mask
Input_Mask_File <- PCA_Output$MaskPath

# 3- Select principal components from the PCA raster
#select_PCA_components(Input_Image_File,Output_Dir,PCA.Files)
# Select components from the PCA/SPCA/MNF raster
# Sel_PC = path of the file where selected components are stored
stop('lo')
Sel_PC <- select_PCA_components(Input_Image_File = Input_Image_File, Output_Dir = Output_Dir, PCA_Files = PCA_Files, TypePCA = TypePCA, File_Open = TRUE)

################################################################################
##                      MAP ALPHA AND BETA DIVERSITY                          ##
################################################################################
print("MAP SPECTRAL SPECIES")

Kmeans_info <- map_spectral_species(Input_Image_File = Input_Image_File, Output_Dir = Output_Dir, PCA_Files = PCA_Files, Input_Mask_File = Input_Mask_File, Pix_Per_Partition = Pix_Per_Partition, nb_partitions = nb_partitions, nbCPU = nbCPU, MaxRAM = MaxRAM, nbclusters = nbclusters, TypePCA = TypePCA)


print("MAP ALPHA DIVERSITY")
Index_Alpha <- c('Shannon')
alpha_div <- map_alpha_div(Input_Image_File = Input_Image_File, Output_Dir = Output_Dir, TypePCA = TypePCA, window_size = window_size, nbCPU = nbCPU, MaxRAM = MaxRAM, Index_Alpha = Index_Alpha, nbclusters = nbclusters)


#alpha_file <- raster::raster(alpha_div)
#alpha_map <- mapview::mapview(alpha_file, layer.name = "Estimated Shannon index", col.regions = brewer.pal(7, "Dark2")) + mapview::mapview(alphamean_file, legend = FALSE, col.regions = brewer.pal(7, "Dark2"))
#return(alpha_map)
##### Trouver le moyen de mapview le truc la 


print("MAP BETA DIVERSITY")
beta_div <- map_beta_div(Input_Image_File = Input_Image_File, Output_Dir = Output_Dir, TypePCA = TypePCA, window_size = window_size, nb_partitions=nb_partitions, nbCPU = nbCPU, MaxRAM = MaxRAM, nbclusters = nbclusters)

#beta_file <- raster::raster(beta_div)
#beta_map <- mapview::mapview(beta_file, layer.name = "PCoA", col.regions = c("red", "blue", "green")) + mapview::mapview(betafull_file, legend = FALSE, col.regions = c("red", "blue", "green"))
#return(beta_map)
##### Trouver le moyen de mapview le truc la 


################################################################################
##          COMPUTE ALPHA AND BETA DIVERSITY FROM FIELD PLOTS                 ##
################################################################################
## read selected features from dimensionality reduction 
Selected_Features <- read.table(Sel_PC)[[1]]
## path for selected components
mapper <- map_functional_div(Original_Image_File = Input_Image_File, Functional_File = PCA_Files,  Selected_Features = Selected_Features, Output_Dir = Output_Dir, window_size = window_size, nbCPU = nbCPU, MaxRAM = MaxRAM,TypePCA = TypePCA)

# location of the directory where shapefiles used for validation are saved
VectorDir <- destunz
# list vector data
Path_Vector <- list_shp(VectorDir)
Name_Vector <- tools::file_path_sans_ext(basename(Path_Vector))
# location of the spectral species raster needed for validation
Path_SpectralSpecies <- Kmeans_info$SpectralSpecies
# get diversity indicators corresponding to shapefiles (no partitioning of spectral dibversity based on field plots so far...)
Biodiv_Indicators <- diversity_from_plots(Raster_SpectralSpecies = Path_SpectralSpecies, Plots = Path_Vector, nbclusters = nbclusters, Raster_Functional = PCA_Files, Selected_Features = Selected_Features)

Shannon_RS <- c(Biodiv_Indicators$Shannon)[[1]]
FRic <- c(Biodiv_Indicators$FunctionalDiversity$FRic)
FEve <- c(Biodiv_Indicators$FunctionalDiversity$FEve)
FDiv <- c(Biodiv_Indicators$FunctionalDiversity$FDiv)
# if no name for plots
Biodiv_Indicators$Name_Plot = seq(1, length(Biodiv_Indicators$Shannon[[1]]), by = 1)


####################################################
# write RS indicators
####################################################
# write a table for Shannon index

write.table(Shannon_RS, file = "ShannonIndex.csv", sep="\t", dec=".", na=" ", row.names = Biodiv_Indicators$Name_Plot, col.names= F, quote=FALSE)

            
# write a table for all spectral diversity indices corresponding to alpha diversity
Results <- data.frame(Name_Vector, Biodiv_Indicators$Richness, Biodiv_Indicators$Fisher,
                      Biodiv_Indicators$Shannon, Biodiv_Indicators$Simpson,
                      Biodiv_Indicators$FunctionalDiversity$FRic,
                      Biodiv_Indicators$FunctionalDiversity$FEve,
                      Biodiv_Indicators$FunctionalDiversity$FDiv)

names(Results)  = c("ID_Plot", "Species_Richness", "Fisher", "Shannon", "Simpson", "FRic", "FEve", "FDiv")
write.table(Results, file = "AlphaDiversity.csv", sep="\t", dec=".", na=" ", row.names = F, col.names= T, quote=FALSE)
  

# write a table for Bray Curtis dissimilarity
BC_mean <- Biodiv_Indicators$BCdiss
colnames(BC_mean) <- rownames(BC_mean) <- Biodiv_Indicators$Name_Plot
write.table(BC_mean, file = "BrayCurtis.csv", sep="\t", dec=".", na=" ", row.names = F, col.names= T, quote=FALSE)



####################################################
# illustrate results
####################################################
# apply ordination using PCoA (same as done for map_beta_div)

MatBCdist <- as.dist(BC_mean, diag = FALSE, upper = FALSE)
BetaPCO <- labdsv::pco(MatBCdist, k = 3)
#plot(BetaPCO)
# assign a type of vegetation to each plot, assuming that the type of vegetation 
# is defined by the name of the shapefile         
 
nbSamples <- shpName <- c()
for (i in 1:length(Path_Vector)){
  shp <- Path_Vector[i]
  nbSamples[i] <- length(rgdal::readOGR(shp,verbose = FALSE))
  shpName[i] <- file_path_sans_ext(basename(shp))
}

Type_Vegetation = c()
for (i in 1: length(nbSamples)){
  for (j in 1:nbSamples[i]){
    Type_Vegetation = c(Type_Vegetation,shpName[i])
  }
}

#  data frame including a selection of alpha diversity metrics and beta diversity expressed as coordinates in the PCoA space
Results <- data.frame('vgtype' = Type_Vegetation, 'pco1' = BetaPCO$points[,1], 'pco2' = BetaPCO$points[,2], 'pco3' = BetaPCO$points[,3], 'shannon' = Shannon_RS, 'FRic' = FRic, 'FEve' = FEve, 'FDiv' = FDiv)
                      
# plot field data in the PCoA space, with size corresponding to shannon index
g1 <- ggplot2::ggplot(Results, aes(x = pco1, y = pco2, color = vgtype, size = shannon)) + geom_point(alpha = 0.6) + scale_color_manual(values = c("#e6140a", "#e6d214", "#e68214", "#145ae6"))

g2 <- ggplot2::ggplot(Results, aes(x = pco1, y = pco3, color = vgtype, size = shannon)) + geom_point(alpha = 0.6) + scale_color_manual(values = c("#e6140a", "#e6d214", "#e68214", "#145ae6"))

g3 <- ggplot2::ggplot(Results, aes(x = pco2, y = pco3, color = vgtype, size = shannon)) + geom_point(alpha = 0.6) + scale_color_manual(values = c("#e6140a", "#e6d214", "#e68214", "#145ae6"))                     

#extract legend                 
get_legend <- function(a.gplot) {
    tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

legend <- get_legend(g3)
gAll <- gridExtra::grid.arrange(arrangeGrob(g1 + theme(legend.position = "none"), g2 + theme(legend.position = "none"), g3 + theme(legend.position = "none"), nrow=1), legend, nrow = 2, heights = c(5, 4)) 


filename <- ggplot2::ggsave("BetaDiversity_PcoA1_vs_PcoA2_vs_PcoA3.png", gAll, scale = 1, width = 12, height = 7, units = "in", dpi = 600, limitsize = TRUE)

filename


